#!/usr/bin/env python3
"""Rank countries by circularity (isoperimetric quotient 4*pi*A/P^2; 1=circle).
Computed on each country's local equal-area projection. Reports both the
mainland-only and whole-country (all islands) versions."""
import json, math
from shapely.geometry import shape
from shapely.ops import unary_union
from pyproj import Transformer

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

def ipq(poly):
    return 4 * math.pi * poly.area / (poly.length ** 2)

rows = []
for name, gs in geoms.items():
    g = unary_union(gs)
    lon, lat = g.centroid.x, g.centroid.y
    tf = Transformer.from_crs('EPSG:4326',
        f'+proj=laea +lat_0={lat} +lon_0={lon} +ellps=WGS84', always_xy=True)
    from shapely.ops import transform
    gp = transform(lambda x, y, z=None: tf.transform(x, y), g)
    main = gp if gp.geom_type == 'Polygon' else max(gp.geoms, key=lambda x: x.area)
    rows.append((name, ipq(main), ipq(gp), main.area / 1e6, gp.area / 1e6))

print('=== MOST CIRCULAR (mainland polygon only) ===')
for r, x in enumerate(sorted(rows, key=lambda r: -r[1])[:15], 1):
    print(f'{r:>2}. {x[1]:.3f}  {x[0]:<22} ({x[3]:,.0f} km^2)')

print('\n=== MOST CIRCULAR (whole country incl. all islands) ===')
for r, x in enumerate(sorted(rows, key=lambda r: -r[2])[:15], 1):
    print(f'{r:>2}. {x[2]:.3f}  {x[0]:<22} ({x[4]:,.0f} km^2)')

print('\n=== MOST CIRCULAR, substantial only (mainland, >= 100,000 km^2) ===')
big = [r for r in rows if r[3] >= 100_000]
for r, x in enumerate(sorted(big, key=lambda r: -r[1])[:15], 1):
    print(f'{r:>2}. {x[1]:.3f}  {x[0]:<22} ({x[3]:,.0f} km^2)')

# ---- visualization ---------------------------------------------------------
import numpy as np
from shapely.ops import transform

def projected_main(name, gs):
    g = unary_union(gs)
    lon, lat = g.centroid.x, g.centroid.y
    tf = Transformer.from_crs('EPSG:4326',
        f'+proj=laea +lat_0={lat} +lon_0={lon} +ellps=WGS84', always_xy=True)
    gp = transform(lambda x, y, z=None: tf.transform(x, y), g)
    return gp if gp.geom_type == 'Polygon' else max(gp.geoms, key=lambda x: x.area)

def draw_panel(svg, name, score, poly, cx, cy, box, label_y_top, label_y_bot):
    pts = np.array(poly.exterior.coords)
    c = pts.mean(0)
    pts = pts - c
    rad = np.sqrt((pts**2).sum(1)).max()
    sc = (box*0.34) / rad
    d = 'M' + 'L'.join(f'{cx+sc*x:.1f},{cy-sc*y:.1f}' for x, y in pts) + 'Z'
    # equal-area reference circle (dashed)
    r_eq = math.sqrt(poly.area / math.pi) * sc
    svg.append(f'<circle cx="{cx}" cy="{cy}" r="{r_eq:.1f}" fill="none" '
               f'stroke="#9ca3af" stroke-width="1.2" stroke-dasharray="4 3"/>')
    svg.append(f'<path d="{d}" fill="#3b82f6" fill-opacity="0.40" stroke="#1e3a8a" stroke-width="1.6"/>')
    svg.append(f'<text x="{cx}" y="{label_y_top}" font-size="13" text-anchor="middle" font-weight="bold" fill="#1e3a8a">{name}</text>')
    svg.append(f'<text x="{cx}" y="{label_y_bot}" font-size="12" text-anchor="middle" fill="#555">IPQ={score:.3f}</text>')

def make_svg(items, title, fname):
    cols, cell = 5, 220
    rows_n = (len(items)+cols-1)//cols
    W, H = cols*cell, rows_n*cell + 40
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="sans-serif" style="background:#fff">']
    svg.append(f'<text x="10" y="24" font-size="16" font-weight="bold">{title}</text>')
    for k, (name, score, poly) in enumerate(items):
        r, col = divmod(k, cols)
        cx, cy = col*cell+cell/2, r*cell+cell/2+40
        draw_panel(svg, name, score, poly, cx, cy, cell, r*cell+40+16, (r+1)*cell+40-8)
    svg.append('</svg>')
    open(fname, 'w').write('\n'.join(svg))
    import cairosvg
    cairosvg.svg2png(url=fname, write_to=fname.replace('.svg', '.png'), scale=2.0)
    print('written', fname)

main_polys = {name: projected_main(name, gs) for name, gs in geoms.items()}

top_overall = sorted(rows, key=lambda r: -r[1])[:10]
make_svg([(x[0], x[1], main_polys[x[0]]) for x in top_overall],
         'Most circular countries (dashed grey = equal-area reference circle; IPQ 1.0 = perfect circle)',
         'most_circular.svg')

top_big = sorted([r for r in rows if r[3] >= 100_000], key=lambda r: -r[1])[:10]
make_svg([(x[0], x[1], main_polys[x[0]]) for x in top_big],
         'Most circular substantial countries (>= 100,000 km^2)',
         'most_circular_substantial.svg')
