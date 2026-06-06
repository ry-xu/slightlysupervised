#!/usr/bin/env python3
"""Rank countries by how much they look like a lowercase 'w' (Liberation
Serif, the open Times-New-Roman-compatible font). Score = best IoU between the
country and the glyph under a similarity transform + reflection (translation,
rotation, uniform scale, mirror). 1.0 = identical to the letter."""
import json, math
import numpy as np
from shapely.geometry import shape, Polygon
from shapely.ops import unary_union, transform
from pyproj import Transformer
from scipy.optimize import minimize

# ---- target glyph, normalized to centroid origin & unit area --------------
GP = np.array(json.load(open('w_glyph.json'))[0])
def normalize(P):
    poly = Polygon(P)
    P = P - np.array(poly.centroid.coords)[0]
    P = P / math.sqrt(abs(Polygon(P).area))
    return P
GP = normalize(GP)
GLYPH_AREA = 1.0

def place(P, tx, ty, theta, s, mirror):
    Q = P.copy()
    if mirror:
        Q = Q * np.array([-1.0, 1.0])
    c, sn = math.cos(theta), math.sin(theta)
    R = np.array([[c, -sn], [sn, c]])
    Q = (Q @ R.T) * s
    Q[:, 0] += tx; Q[:, 1] += ty
    return Q

# ---- countries ------------------------------------------------------------
data = json.load(open('countries_50m.geojson'))
SKIP = {'Antarctica'}
MERGE = {'Somaliland': 'Somalia'}
geoms = {}
for feat in data['features']:
    p = feat['properties']
    name = p.get('NAME') or p.get('ADMIN')
    if name in SKIP or p.get('TYPE') not in ('Sovereign country', 'Country', 'Dependency'):
        continue
    geoms.setdefault(MERGE.get(name, name), []).append(shape(feat['geometry']))

def best_w(poly):
    A = poly.area
    def neg_iou(p):
        Q = place(GP, p[0], p[1], p[2], math.exp(p[3]), p[4] > 0)
        g = Polygon(Q)
        if not g.is_valid:
            g = g.buffer(0)
        inter = poly.intersection(g).area
        return -inter / (A + g.area - inter)
    best = None
    for mirror in (0, 1):
        for k in range(18):                         # coarse rotation scan
            theta = 2*math.pi*k/18
            for ls in (-0.15, 0.0, 0.18):           # coarse scale scan
                f = neg_iou([0, 0, theta, ls, mirror])
                if best is None or f < best[0]:
                    best = (f, [0, 0, theta, ls, mirror])
    # local refine (keep reflection fixed)
    x0 = best[1]; m = x0[4]
    res = minimize(lambda q: neg_iou([q[0], q[1], q[2], q[3], m]), x0[:4],
                   method='Nelder-Mead',
                   options={'xatol': 0.01, 'fatol': 1e-4, 'maxiter': 4000})
    if -res.fun > -best[0]:
        best = (res.fun, list(res.x) + [m])
    return best[1], -best[0]

rows = []
for name, gs in geoms.items():
    g = unary_union(gs)
    lon, lat = g.centroid.x, g.centroid.y
    tf = Transformer.from_crs('EPSG:4326',
        f'+proj=laea +lat_0={lat} +lon_0={lon} +ellps=WGS84', always_xy=True)
    gp = transform(lambda x, y, z=None: tf.transform(x, y), g)
    poly = gp if gp.geom_type == 'Polygon' else max(gp.geoms, key=lambda x: x.area)
    poly = poly.simplify(poly.length / 600)         # speed; keep shape
    # normalize country to centroid origin & unit area to match glyph scale
    cen = np.array(poly.centroid.coords)[0]
    norm = transform(lambda x, y, z=None:
                     ((np.asarray(x)-cen[0])/math.sqrt(poly.area),
                      (np.asarray(y)-cen[1])/math.sqrt(poly.area)), poly)
    params, iou = best_w(norm)
    rows.append({'name': name, 'poly': norm, 'params': params, 'iou': iou,
                 'area': poly.area/1e6})

rows.sort(key=lambda d: -d['iou'])

print("=== MOST LIKE A LOWERCASE 'w' (best IoU; sim transform + mirror) ===")
for i, d in enumerate(rows[:20], 1):
    print(f'{i:>2}. {d["iou"]:.3f}  {d["name"]:<22} ({d["area"]:,.0f} km^2)')

big = [d for d in rows if d['area'] >= 100_000]
print('\n=== substantial only (>= 100,000 km^2) ===')
for i, d in enumerate(big[:20], 1):
    print(f'{i:>2}. {d["iou"]:.3f}  {d["name"]:<22} ({d["area"]:,.0f} km^2)')

import csv
with open('w_likeness_iou.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['rank', 'iou', 'country', 'area_km2'])
    for i, d in enumerate(rows, 1):
        w.writerow([i, f'{d["iou"]:.4f}', d['name'], f'{d["area"]:.0f}'])

# ---- visualization --------------------------------------------------------
import cairosvg
def make_svg(items, title, fname):
    cols, cell = 5, 230
    rn = (len(items)+cols-1)//cols
    W, H = cols*cell, rn*cell + 40
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="sans-serif" style="background:#fff">']
    svg.append(f'<text x="10" y="24" font-size="16" font-weight="bold">{title}</text>')
    for k, d in enumerate(items):
        r, col = divmod(k, cols)
        cx0, cy0 = col*cell+cell/2, r*cell+cell/2+40
        P = np.array(d['poly'].exterior.coords)
        G = place(GP, *d['params'][:3], math.exp(d['params'][3]), d['params'][4] > 0)
        allp = np.vstack([P, G])
        cen = allp.mean(0); rad = np.sqrt(((allp-cen)**2).sum(1)).max(); sc = cell*0.34/rad
        def tp(px, py): return f'{cx0+sc*(px-cen[0]):.1f},{cy0-sc*(py-cen[1]):.1f}'
        path = 'M'+'L'.join(tp(px, py) for px, py in P)+'Z'
        gpath = 'M'+'L'.join(tp(px, py) for px, py in G)+'Z'
        svg.append(f'<path d="{path}" fill="#dc2626" fill-opacity="0.34" stroke="#7f1d1d" stroke-width="1.6"/>')
        svg.append(f'<path d="{gpath}" fill="#1d4ed8" fill-opacity="0.20" stroke="#1d4ed8" stroke-width="1.8"/>')
        svg.append(f'<text x="{cx0}" y="{r*cell+40+16}" font-size="13" text-anchor="middle" font-weight="bold" fill="#7f1d1d">{k+1}. {d["name"]}</text>')
        svg.append(f'<text x="{cx0}" y="{(r+1)*cell+40-8}" font-size="12" text-anchor="middle" fill="#555">IoU={d["iou"]:.3f}</text>')
    svg.append('</svg>')
    open(fname, 'w').write('\n'.join(svg))
    cairosvg.svg2png(url=fname, write_to=fname.replace('.svg', '.png'), scale=2.0)
    print('written', fname.replace('.svg', '.png'))

make_svg(rows[:10], "Top 10 countries most like a lowercase 'w' (red = country, blue = best-fit 'w')", 'top_w.svg')
make_svg(big[:10], "Top 10 substantial countries most like 'w' (>= 100,000 km^2)", 'top_w_substantial.svg')
