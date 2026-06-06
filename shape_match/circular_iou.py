#!/usr/bin/env python3
"""Rank countries by how circular they are, measured as the overlap (IoU)
with their BEST-FIT circle (center + radius optimized to maximize
intersection-over-union). 1.0 = perfect circle."""
import json, math
import numpy as np
from shapely.geometry import shape, Point
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

def best_circle(poly):
    A = poly.area
    P = np.array(poly.exterior.coords)
    x, y = P[:, 0], P[:, 1]
    M = np.c_[2*x, 2*y, np.ones(len(x))]
    cx0, cy0, c = np.linalg.lstsq(M, x**2 + y**2, rcond=None)[0]
    r0 = math.sqrt(max(c + cx0**2 + cy0**2, 1.0))
    def neg_iou(pp, qs=48):
        cx, cy, r = pp
        if r <= 0:
            return 1.0
        disk = Point(cx, cy).buffer(r, quad_segs=qs)
        inter = poly.intersection(disk).area
        return -inter / (A + disk.area - inter)
    res = minimize(neg_iou, [cx0, cy0, r0], method='Nelder-Mead',
                   options={'xatol': r0*0.002, 'fatol': 1e-4, 'maxiter': 2000})
    cx, cy, r = res.x
    iou = -neg_iou((cx, cy, r), qs=256)   # high-precision final score
    return cx, cy, r, iou

rows = []
for name, gs in geoms.items():
    g = unary_union(gs)
    lon, lat = g.centroid.x, g.centroid.y
    tf = Transformer.from_crs('EPSG:4326',
        f'+proj=laea +lat_0={lat} +lon_0={lon} +ellps=WGS84', always_xy=True)
    gp = transform(lambda x, y, z=None: tf.transform(x, y), g)
    poly = gp if gp.geom_type == 'Polygon' else max(gp.geoms, key=lambda x: x.area)
    cx, cy, r, iou = best_circle(poly)
    rows.append({'name': name, 'poly': poly, 'cx': cx, 'cy': cy, 'r': r,
                 'iou': iou, 'area': poly.area/1e6})

rows.sort(key=lambda d: -d['iou'])

print('=== MOST CIRCULAR by best-fit-circle IoU (mainland) ===')
for i, d in enumerate(rows[:20], 1):
    print(f'{i:>2}. {d["iou"]:.3f}  {d["name"]:<22} ({d["area"]:,.0f} km^2)')

big = [d for d in rows if d['area'] >= 100_000]
print('\n=== substantial only (>= 100,000 km^2) ===')
for i, d in enumerate(big[:20], 1):
    print(f'{i:>2}. {d["iou"]:.3f}  {d["name"]:<22} ({d["area"]:,.0f} km^2)')

import csv
with open('circularity_iou.csv', 'w', newline='') as f:
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
        cen = np.array([d['cx'], d['cy']])
        pts = P - cen
        rad = max(np.sqrt((pts**2).sum(1)).max(), d['r'])
        sc = cell*0.33/rad
        path = 'M'+'L'.join(f'{cx0+sc*px:.1f},{cy0-sc*py:.1f}' for px, py in pts)+'Z'
        svg.append(f'<circle cx="{cx0:.1f}" cy="{cy0:.1f}" r="{d["r"]*sc:.1f}" fill="#2563eb" fill-opacity="0.12" stroke="#1d4ed8" stroke-width="1.8" stroke-dasharray="5 3"/>')
        svg.append(f'<path d="{path}" fill="#dc2626" fill-opacity="0.32" stroke="#7f1d1d" stroke-width="1.6"/>')
        svg.append(f'<text x="{cx0}" y="{r*cell+40+16}" font-size="13" text-anchor="middle" font-weight="bold" fill="#7f1d1d">{k+1}. {d["name"]}</text>')
        svg.append(f'<text x="{cx0}" y="{(r+1)*cell+40-8}" font-size="12" text-anchor="middle" fill="#555">IoU={d["iou"]:.3f}</text>')
    svg.append('</svg>')
    open(fname, 'w').write('\n'.join(svg))
    cairosvg.svg2png(url=fname, write_to=fname.replace('.svg', '.png'), scale=2.0)
    print('written', fname.replace('.svg', '.png'))

make_svg(rows[:10], 'Top 10 most circular countries (red = country, blue = best-fit circle; IoU 1.0 = perfect)', 'top_circular.svg')
make_svg(big[:10], 'Top 10 most circular substantial countries (>= 100,000 km^2)', 'top_circular_substantial.svg')
