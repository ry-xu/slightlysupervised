#!/usr/bin/env python3
"""Rank countries by how triangular they are: overlap (IoU) with their
best-fit EQUILATERAL triangle (center, rotation and size optimized).
1.0 = perfect equilateral triangle."""
import json, math
import numpy as np
from shapely.geometry import shape, Polygon
from shapely.ops import unary_union, transform
from pyproj import Transformer
from scipy.optimize import minimize

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

def equilateral(cx, cy, phi, R):
    return [(cx + R*math.cos(phi + math.radians(a)),
             cy + R*math.sin(phi + math.radians(a))) for a in (90, 210, 330)]

def best_triangle(poly, restarts=6):
    """Best-fit EQUILATERAL triangle: params = (cx, cy, phi, R)."""
    A = poly.area
    c = (poly.centroid.x, poly.centroid.y)
    R0 = math.sqrt(4 * A / (3 * math.sqrt(3)))     # circumradius giving area ~ A
    best = None
    def neg_iou(p):
        cx, cy, phi, R = p
        if R <= 0:
            return 1.0
        tri = Polygon(equilateral(cx, cy, phi, R))
        inter = poly.intersection(tri).area
        return -inter / (A + tri.area - inter)
    for k in range(restarts):
        phi0 = (math.pi/3) * k / restarts          # only need to scan 0..60 deg
        res = minimize(neg_iou, [c[0], c[1], phi0, R0], method='Nelder-Mead',
                       options={'xatol': R0*0.003, 'fatol': 1e-4, 'maxiter': 3000})
        if best is None or res.fun < best[0]:
            best = (res.fun, res.x)
    cx, cy, phi, R = best[1]
    return np.array(equilateral(cx, cy, phi, R)).ravel(), -best[0]

rows = []
for name, gs in geoms.items():
    g = unary_union(gs)
    lon, lat = g.centroid.x, g.centroid.y
    tf = Transformer.from_crs('EPSG:4326',
        f'+proj=laea +lat_0={lat} +lon_0={lon} +ellps=WGS84', always_xy=True)
    gp = transform(lambda x, y, z=None: tf.transform(x, y), g)
    poly = gp if gp.geom_type == 'Polygon' else max(gp.geoms, key=lambda x: x.area)
    verts, iou = best_triangle(poly)
    rows.append({'name': name, 'poly': poly, 'verts': verts, 'iou': iou,
                 'area': poly.area/1e6})

rows.sort(key=lambda d: -d['iou'])

print('=== MOST TRIANGULAR by best-fit-triangle IoU (mainland) ===')
for i, d in enumerate(rows[:20], 1):
    print(f'{i:>2}. {d["iou"]:.3f}  {d["name"]:<22} ({d["area"]:,.0f} km^2)')

big = [d for d in rows if d['area'] >= 100_000]
print('\n=== substantial only (>= 100,000 km^2) ===')
for i, d in enumerate(big[:20], 1):
    print(f'{i:>2}. {d["iou"]:.3f}  {d["name"]:<22} ({d["area"]:,.0f} km^2)')

import csv
with open('triangularity_iou.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['rank', 'iou', 'country', 'area_km2'])
    for i, d in enumerate(rows, 1):
        w.writerow([i, f'{d["iou"]:.4f}', d['name'], f'{d["area"]:.0f}'])

# ---- visualization ---------------------------------------------------------
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
        V = np.array(d['verts']).reshape(3, 2)
        cen = P.mean(0)
        allp = np.vstack([P, V]) - cen
        rad = np.sqrt((allp**2).sum(1)).max()
        sc = cell*0.33/rad
        def tp(px, py): return f'{cx0+sc*(px-cen[0]):.1f},{cy0-sc*(py-cen[1]):.1f}'
        path = 'M'+'L'.join(tp(px, py) for px, py in P)+'Z'
        tri = 'M'+'L'.join(tp(px, py) for px, py in V)+'Z'
        svg.append(f'<path d="{tri}" fill="#16a34a" fill-opacity="0.14" stroke="#15803d" stroke-width="1.8" stroke-dasharray="5 3"/>')
        svg.append(f'<path d="{path}" fill="#dc2626" fill-opacity="0.32" stroke="#7f1d1d" stroke-width="1.6"/>')
        svg.append(f'<text x="{cx0}" y="{r*cell+40+16}" font-size="13" text-anchor="middle" font-weight="bold" fill="#7f1d1d">{k+1}. {d["name"]}</text>')
        svg.append(f'<text x="{cx0}" y="{(r+1)*cell+40-8}" font-size="12" text-anchor="middle" fill="#555">IoU={d["iou"]:.3f}</text>')
    svg.append('</svg>')
    open(fname, 'w').write('\n'.join(svg))
    cairosvg.svg2png(url=fname, write_to=fname.replace('.svg', '.png'), scale=2.0)
    print('written', fname.replace('.svg', '.png'))

make_svg(rows[:10], 'Top 10 most triangular countries (red = country, green = best-fit triangle; IoU 1.0 = perfect)', 'top_triangle.svg')
make_svg(big[:10], 'Top 10 most triangular substantial countries (>= 100,000 km^2)', 'top_triangle_substantial.svg')
