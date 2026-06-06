#!/usr/bin/env python3
"""Alternative 'most like a w' ranking using the OUTLINE-shape metric:
scale/rotation/reflection-invariant Fourier descriptors of the boundary
(same metric used for country-vs-country matching). Compares each country's
outline to the lowercase 'w' glyph outline."""
import json, math
import numpy as np
from shapely.geometry import shape, Polygon
from shapely.ops import unary_union, transform
from pyproj import Transformer

N, K = 400, 24      # more harmonics: the 'w' has fine zig-zag detail

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

# target glyph
GP = np.array(json.load(open('w_glyph.json'))[0])
gpts = resample(GP, N)
gdesc = descriptor(gpts)

data = json.load(open('countries_50m.geojson'))
SKIP = {'Antarctica'}; MERGE = {'Somaliland': 'Somalia'}
geoms = {}
for feat in data['features']:
    p = feat['properties']
    name = p.get('NAME') or p.get('ADMIN')
    if name in SKIP or p.get('TYPE') not in ('Sovereign country', 'Country', 'Dependency'):
        continue
    geoms.setdefault(MERGE.get(name, name), []).append(shape(feat['geometry']))

rows = []
for name, gs in geoms.items():
    g = unary_union(gs)
    lon, lat = g.centroid.x, g.centroid.y
    tf = Transformer.from_crs('EPSG:4326',
        f'+proj=laea +lat_0={lat} +lon_0={lon} +ellps=WGS84', always_xy=True)
    gp = transform(lambda x, y, z=None: tf.transform(x, y), g)
    poly = gp if gp.geom_type == 'Polygon' else max(gp.geoms, key=lambda x: x.area)
    pts = resample(np.array(poly.exterior.coords), N)
    d = descriptor(pts)
    dm = descriptor(pts * np.array([-1, 1]))
    dist = min(np.linalg.norm(gdesc - d), np.linalg.norm(gdesc - dm))
    rows.append({'name': name, 'pts': pts, 'dist': dist, 'area': poly.area/1e6})

rows.sort(key=lambda d: d['dist'])
print("=== MOST LIKE 'w' by OUTLINE descriptor (lower = closer) ===")
for i, d in enumerate(rows[:20], 1):
    print(f'{i:>2}. {d["dist"]:.3f}  {d["name"]:<22} ({d["area"]:,.0f} km^2)')
big = [d for d in rows if d['area'] >= 100_000]
print('\n=== substantial only (>= 100,000 km^2) ===')
for i, d in enumerate(big[:15], 1):
    print(f'{i:>2}. {d["dist"]:.3f}  {d["name"]:<22} ({d["area"]:,.0f} km^2)')

# visualization: align each country outline to the w glyph
def norm(pts):
    z = pts[:, 0] + 1j*pts[:, 1]; z -= z.mean(); z /= np.sqrt((np.abs(z)**2).mean()); return z
def align(za, zb):
    best = None
    for zr in (zb, np.conj(zb)):
        c = np.fft.ifft(np.fft.fft(za) * np.conj(np.fft.fft(zr)))
        m = np.argmax(np.abs(c)); ph = c[m]/(abs(c[m])+1e-12)
        cand = np.roll(zr, -m)*ph; cost = np.abs(za-cand).sum()
        if best is None or cost < best[0]: best = (cost, cand)
    return best[1]

import cairosvg
zw = norm(gpts)
def make_svg(items, title, fname):
    cols, cell = 5, 230; rn = (len(items)+cols-1)//cols
    W, H = cols*cell, rn*cell+40
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="sans-serif" style="background:#fff">']
    svg.append(f'<text x="10" y="24" font-size="16" font-weight="bold">{title}</text>')
    for k, d in enumerate(items):
        r, col = divmod(k, cols); cx0, cy0 = col*cell+cell/2, r*cell+cell/2+40
        zb = align(zw, norm(d['pts'])); sc = cell*0.32
        def pth(z): return 'M'+'L'.join(f'{cx0+sc*v.real:.1f},{cy0-sc*v.imag:.1f}' for v in z)+'Z'
        svg.append(f'<path d="{pth(zw)}" fill="none" stroke="#1d4ed8" stroke-width="2"/>')
        svg.append(f'<path d="{pth(zb)}" fill="#dc2626" fill-opacity="0.34" stroke="#7f1d1d" stroke-width="1.5"/>')
        svg.append(f'<text x="{cx0}" y="{r*cell+40+16}" font-size="13" text-anchor="middle" font-weight="bold" fill="#7f1d1d">{k+1}. {d["name"]}</text>')
        svg.append(f'<text x="{cx0}" y="{(r+1)*cell+40-8}" font-size="12" text-anchor="middle" fill="#555">dist={d["dist"]:.3f}</text>')
    svg.append('</svg>'); open(fname, 'w').write('\n'.join(svg))
    cairosvg.svg2png(url=fname, write_to=fname.replace('.svg', '.png'), scale=2.0)
    print('written', fname.replace('.svg', '.png'))

make_svg(rows[:10], "Top 10 most 'w'-like outlines (blue = 'w', red = country)", 'top_w_outline.svg')
make_svg(big[:10], "Top 10 most 'w'-like substantial countries (outline)", 'top_w_outline_substantial.svg')
