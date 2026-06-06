# Country shape look-alikes

Finds the country outlines that most resemble each other, invariant to
translation, scale, rotation and reflection.

## Method
1. **Data:** Natural Earth 50m admin-0 borders (`countries_50m.geojson`).
2. **Shape:** largest polygon (mainland) per country; Antarctica excluded.
3. **Projection:** each country reprojected into its own local Lambert
   azimuthal *equal-area* frame (centered on its centroid) to remove lat/lon
   distortion and make projected area = true km².
4. **Descriptor:** outline resampled to 400 arc-length-even points → FFT of
   `z = x + iy` → drop DC (translation), divide by the fundamental (scale),
   compare harmonic magnitudes (rotation + start-point). Reflection handled
   by also comparing the mirrored descriptor.
5. **Distance:** Euclidean between descriptors; lower = more alike.

## Run
    pip install numpy scipy shapely pyproj cairosvg
    python3 match.py

## Outputs
- `top_pairs.svg` / `top_pairs.png` — aligned overlays of the top 12 pairs.
- `rankings_substantial.csv` — all pairs among countries ≥ 100,000 km².
- `rankings_all.csv` — all 26,796 pairs (includes micro-states).

## Headline result
Among substantial countries, the most alike pair is **Saudi Arabia ↔ Bulgaria**
(d = 0.055), followed by **Nigeria ↔ Ethiopia** and **Mongolia ↔ Honduras**.
Without a size filter the ranking is dominated by near-circular micro-states
(e.g. Wallis & Futuna ↔ Comoros), whose blob outlines trivially resemble each
other.
