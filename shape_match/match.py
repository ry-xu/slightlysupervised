#!/usr/bin/env python3
"""Find the countries whose outlines look most alike, invariant to
translation, scale, rotation and reflection.

Pipeline:
  1. Take the largest polygon (mainland) of each country.
  2. Project it into a local Lambert Azimuthal Equal-Area frame centered on
     the country's own centroid (equal-area => projected area is true km^2,
     and local shape is undistorted).
  3. Resample the outline to N points equally spaced by arc length (CCW).
  4. Fourier descriptor: FFT of z=x+iy; drop DC (translation), divide by the
     fundamental magnitude (scale), compare HARMONIC MAGNITUDES (invariant to
     rotation + start point). Reflection handled by also comparing the mirror.
  5. Rank every country pair by descriptor distance.
"""
import json
import numpy as np
from shapely.geometry import shape
from shapely.ops import orient
from pyproj import Transformer

N = 400          # resampled contour points
K = 12           # harmonics kept on each side of the spectrum
AREA_MIN = 100_000   # km^2 ; "substantial" countries for the headline ranking

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

countries = []
for feat in data['features']:
    p = feat['properties']
    name = p.get('NAME') or p.get('ADMIN')
    if name in SKIP or p.get('TYPE') not in ('Sovereign country', 'Country', 'Dependency'):
        continue
    poly = largest_polygon(feat['geometry'])
    lon, lat = poly.centroid.x, poly.centroid.y
    tf = Transformer.from_crs('EPSG:4326',
        f'+proj=laea +lat_0={lat} +lon_0={lon} +ellps=WGS84', always_xy=True)
    px, py = tf.transform(*poly.exterior.coords.xy)
    proj = orient(shape({'type': 'Polygon', 'coordinates': [list(zip(px, py))]}), 1.0)
    pts = resample(np.array(proj.exterior.coords), N)
    countries.append({
        'name': name, 'cont': p.get('CONTINENT'),
        'area': proj.area / 1e6,            # km^2 (equal-area projection)
        'pts': pts,
        'd': descriptor(pts),
        'dm': descriptor(pts * np.array([-1, 1])),
    })

print(f'{len(countries)} countries processed')

def rank(subset):
    out = []
    for i in range(len(subset)):
        a = subset[i]
        for j in range(i+1, len(subset)):
            b = subset[j]
            dist = min(np.linalg.norm(a['d'] - b['d']),
                       np.linalg.norm(a['d'] - b['dm']))
            out.append((dist, a, b))
    out.sort(key=lambda x: x[0])
    return out

big = [c for c in countries if c['area'] >= AREA_MIN]
all_pairs = rank(countries)
big_pairs = rank(big)

print(f'\n=== TOP 30: SUBSTANTIAL COUNTRIES (land area >= {AREA_MIN:,} km^2, n={len(big)}) ===')
print(f'{"#":>3}  {"dist":>6}  {"country A":<20} {"country B":<20}')
for r, (d, a, b) in enumerate(big_pairs[:30], 1):
    print(f'{r:>3}  {d:6.3f}  {a["name"]:<20} {b["name"]:<20} ({a["cont"]}/{b["cont"]})')

import csv
for fname, pairs in [('rankings_all.csv', all_pairs), ('rankings_substantial.csv', big_pairs)]:
    with open(fname, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rank', 'distance', 'country_a', 'area_a_km2', 'country_b', 'area_b_km2'])
        for r, (d, a, b) in enumerate(pairs, 1):
            w.writerow([r, f'{d:.5f}', a['name'], f'{a["area"]:.0f}', b['name'], f'{b["area"]:.0f}'])
print(f'\nCSVs written: rankings_all.csv ({len(all_pairs)}), rankings_substantial.csv ({len(big_pairs)})')

# ---- visualization: align each top pair and emit one combined SVG ----------
def normalize(pts):
    z = pts[:, 0] + 1j * pts[:, 1]
    z -= z.mean()
    z /= np.sqrt((np.abs(z)**2).mean())
    return z

def best_align(za, zb):
    """Return zb rotated/shifted/reflected to best overlay za (complex Procrustes)."""
    best = None
    for zr in (zb, np.conj(zb)):                 # try both reflections
        corr = np.fft.ifft(np.fft.fft(za) * np.conj(np.fft.fft(zr)))
        m = np.argmax(np.abs(corr))
        phase = corr[m] / (abs(corr[m]) + 1e-12)
        cand = np.roll(zr, -m) * phase           # apply shift + rotation
        cost = np.abs(za - cand).sum()
        if best is None or cost < best[0]:
            best = (cost, cand)
    return best[1]

def path(z, cx, cy, scale):
    p = [f'{cx + scale*v.real:.1f},{cy - scale*v.imag:.1f}' for v in z]
    return 'M' + 'L'.join(p) + 'Z'

TOP = 12
cols, cell, pad = 3, 230, 16
rows = (TOP + cols - 1) // cols
W, H = cols*cell, rows*cell + 30
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'font-family="sans-serif" style="background:#fff">']
svg.append(f'<text x="10" y="20" font-size="15" font-weight="bold">'
           f'Top {TOP} look-alike country pairs (aligned for scale/rotation/flip)</text>')
for k, (d, a, b) in enumerate(big_pairs[:TOP]):
    r, c = divmod(k, cols)
    cx, cy = c*cell + cell/2, r*cell + cell/2 + 30
    za, zb = normalize(a['pts']), normalize(b['pts'])
    zb = best_align(za, zb)
    sc = cell*0.30
    svg.append(f'<path d="{path(za,cx,cy,sc)}" fill="#3b82f6" fill-opacity="0.35" '
               f'stroke="#1e40af" stroke-width="1.5"/>')
    svg.append(f'<path d="{path(zb,cx,cy,sc)}" fill="#ef4444" fill-opacity="0.30" '
               f'stroke="#991b1b" stroke-width="1.5"/>')
    svg.append(f'<text x="{cx}" y="{r*cell+30+18}" font-size="12" text-anchor="middle" '
               f'fill="#1e40af">{a["name"]}</text>')
    svg.append(f'<text x="{cx}" y="{(r+1)*cell+30-8}" font-size="12" text-anchor="middle" '
               f'fill="#991b1b">{b["name"]} &#183; d={d:.3f}</text>')
svg.append('</svg>')
open('top_pairs.svg', 'w').write('\n'.join(svg))
print('visualization written: top_pairs.svg')
