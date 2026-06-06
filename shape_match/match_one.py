#!/usr/bin/env python3
"""Rank all countries by outline similarity to one TARGET country and render
an aligned-overlay visualization of the top matches."""
import sys, json
import numpy as np
from shapely.geometry import shape
from shapely.ops import orient
from pyproj import Transformer

TARGET = sys.argv[1] if len(sys.argv) > 1 else 'Somalia'
TOPN = int(sys.argv[2]) if len(sys.argv) > 2 else 10
AREA_MIN = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0   # km^2 filter
N, K = 400, 12

data = json.load(open('countries_50m.geojson'))
SKIP = {'Antarctica'}

def largest_polygon(geom):
    g = shape(geom)
    return g if g.geom_type == 'Polygon' else max(g.geoms, key=lambda p: p.area)

def resample(pts, n):
    pts = np.asarray(pts, float)
    if np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    seg = np.r_[pts, pts[:1]]
    d = np.sqrt((np.diff(seg, axis=0)**2).sum(1))
    s = np.r_[0, np.cumsum(d)]
    targets = np.linspace(0, s[-1], n, endpoint=False)
    j = np.clip(np.searchsorted(s, targets) - 1, 0, len(seg) - 2)
    frac = (targets - s[j]) / (s[j+1] - s[j] + 1e-12)
    return seg[j] + frac[:, None] * (seg[j+1] - seg[j])

def descriptor(pts):
    Z = np.fft.fft(pts[:, 0] + 1j * pts[:, 1])
    Z[0] = 0
    f0 = abs(Z[1]) + 1e-12
    idx = list(range(1, K+1)) + list(range(-K, 0))
    return np.array([abs(Z[k % N]) for k in idx]) / f0

# Natural Earth splits some de-facto states out of their parent country.
# Merge them back so we compare full, internationally-recognized borders.
MERGE = {'Somaliland': 'Somalia'}

geoms = {}   # display name -> list of feature geometries to union
for feat in data['features']:
    p = feat['properties']
    name = p.get('NAME') or p.get('ADMIN')
    if name in SKIP or p.get('TYPE') not in ('Sovereign country', 'Country', 'Dependency'):
        continue
    name = MERGE.get(name, name)
    geoms.setdefault(name, []).append(shape(feat['geometry']))

from shapely.ops import unary_union
countries = []
for name, gs in geoms.items():
    merged = unary_union(gs)
    poly = merged if merged.geom_type == 'Polygon' else max(merged.geoms, key=lambda x: x.area)
    lon, lat = poly.centroid.x, poly.centroid.y
    tf = Transformer.from_crs('EPSG:4326',
        f'+proj=laea +lat_0={lat} +lon_0={lon} +ellps=WGS84', always_xy=True)
    px, py = tf.transform(*poly.exterior.coords.xy)
    proj = orient(shape({'type': 'Polygon', 'coordinates': [list(zip(px, py))]}), 1.0)
    pts = resample(np.array(proj.exterior.coords), N)
    countries.append({'name': name, 'area': proj.area/1e6, 'pts': pts,
                      'd': descriptor(pts), 'dm': descriptor(pts*np.array([-1,1]))})

tgt = next(c for c in countries if c['name'] == TARGET)
ranked = []
for c in countries:
    if c['name'] == TARGET or c['area'] < AREA_MIN:
        continue
    dist = min(np.linalg.norm(tgt['d'] - c['d']), np.linalg.norm(tgt['d'] - c['dm']))
    ranked.append((dist, c))
ranked.sort(key=lambda x: x[0])

print(f'Top {TOPN} countries most like {TARGET}:')
for r, (d, c) in enumerate(ranked[:TOPN], 1):
    print(f'{r:>2}. {d:.3f}  {c["name"]}  ({c["area"]:,.0f} km^2)')

# ---- aligned-overlay visualization ----------------------------------------
def normalize(pts):
    z = pts[:, 0] + 1j*pts[:, 1]; z -= z.mean(); z /= np.sqrt((np.abs(z)**2).mean()); return z

def best_align(za, zb):
    best = None
    for zr in (zb, np.conj(zb)):
        corr = np.fft.ifft(np.fft.fft(za) * np.conj(np.fft.fft(zr)))
        m = np.argmax(np.abs(corr)); ph = corr[m]/(abs(corr[m])+1e-12)
        cand = np.roll(zr, -m) * ph
        cost = np.abs(za - cand).sum()
        if best is None or cost < best[0]: best = (cost, cand)
    return best[1]

def path(z, cx, cy, sc):
    return 'M' + 'L'.join(f'{cx+sc*v.real:.1f},{cy-sc*v.imag:.1f}' for v in z) + 'Z'

za = normalize(tgt['pts'])
cols, cell = 5, 220
rows = (TOPN + cols - 1)//cols
W, H = cols*cell, rows*cell + 40
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="sans-serif" style="background:#fff">']
svg.append(f'<text x="10" y="24" font-size="16" font-weight="bold">Top {TOPN} countries most like {TARGET} '
           f'(grey = {TARGET}, red = match; aligned for scale/rotation/flip)</text>')
for k, (d, c) in enumerate(ranked[:TOPN]):
    r, col = divmod(k, cols)
    cx, cy = col*cell + cell/2, r*cell + cell/2 + 40
    zb = best_align(za, normalize(c['pts'])); sc = cell*0.30
    svg.append(f'<path d="{path(za,cx,cy,sc)}" fill="#9ca3af" fill-opacity="0.45" stroke="#374151" stroke-width="1.5"/>')
    svg.append(f'<path d="{path(zb,cx,cy,sc)}" fill="#ef4444" fill-opacity="0.30" stroke="#991b1b" stroke-width="1.8"/>')
    svg.append(f'<text x="{cx}" y="{r*cell+40+16}" font-size="13" text-anchor="middle" fill="#991b1b" font-weight="bold">{k+1}. {c["name"]}</text>')
    svg.append(f'<text x="{cx}" y="{(r+1)*cell+40-8}" font-size="12" text-anchor="middle" fill="#555">d={d:.3f}</text>')
svg.append('</svg>')
suffix = f'_min{int(AREA_MIN/1000)}k' if AREA_MIN else ''
out = f'like_{TARGET.replace(" ","_")}{suffix}.svg'
open(out, 'w').write('\n'.join(svg))
import cairosvg
cairosvg.svg2png(url=out, write_to=out.replace('.svg', '.png'), scale=2.0)
print('written:', out, 'and', out.replace('.svg', '.png'))
