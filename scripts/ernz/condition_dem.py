"""Hydraulic conditioning of the 10 m Ernz Noire DEM.
1. Base: block mean of the 0.5 m LiDAR (dem10_mean.tif).
2. Channel burn: along the AGE watercourse lines (primary, secondary, temporary; geoportail collections 749/23-25) the cell
   takes the block MINIMUM of the 0.5 m LiDAR (dem10_min.tif), so the channel bed is kept when 400 fine cells are averaged.
3. Crossing removal: along each line, a cell that stands >= 0.5 m above both its neighbours within 3 cells (30 m) is a bridge
   deck, culvert roof or road embankment left in the DTM; it is replaced by the linear interpolation of those neighbours.
   Natural steps (Mullerthal waterfalls) are monotonic and are not touched.
4. Stream burn: along each line, ordered from its high end to its low end, the profile is made non increasing (running
   minimum), so the channel has no pits; neighbouring cells are never raised.
5. Residual channel pits at line junctions are raised to their lowest neighbour + 1 mm (iterated).
usage: python condition_dem.py <domain dir> <lines geojson ...> ; writes inputs/dem10_cond.asc and a report"""
import sys, json, numpy as np, geopandas as gpd, rasterio, pandas as pd
from rasterio.features import rasterize
from shapely.geometry import LineString
D = sys.argv[1]; LINES = sys.argv[2:]
with rasterio.open(f"{D}/dem10_mean.tif") as s: dem = s.read(1).astype("float64"); tf = s.transform
with rasterio.open(f"{D}/dem10_min.tif") as s: dmin = s.read(1).astype("float64")
NR, NC = dem.shape; act = dem > -9000; C = tf.a
g = pd.concat([gpd.read_file(f).to_crs(2169).assign(src=f.split("/")[-1]) for f in LINES], ignore_index=True)
g = g.explode(index_parts=False, ignore_index=True); g = g[g.geometry.type == "LineString"]   # multipart stems were dropped before 25 Oct
out = dem.copy(); n_burn = 0; n_fix = 0; fixes = []
def cells(line):
    """ordered unique (row, col) along a line, sampled every 5 m"""
    n = max(int(line.length / 5), 2); pts = [line.interpolate(i / n, normalized=True) for i in range(n + 1)]
    rc = [(int((tf.f - p.y) / C), int((p.x - tf.c) / C)) for p in pts]
    seq = [rc[0]] + [b for a, b in zip(rc, rc[1:]) if b != a]
    return [(r, c) for r, c in seq if 0 <= r < NR and 0 <= c < NC and act[r, c]]
for line in g.geometry:
    rc = cells(line)
    if len(rc) < 2: continue
    r, c = zip(*rc); r, c = np.array(r), np.array(c)
    z = dmin[r, c]; n_burn += len(rc)
    for _ in range(3):                                           # iterate: adjacent crossing cells
        for i in range(1, len(z) - 1):
            lo, hi = max(0, i - 3), min(len(z), i + 4)
            zl, zr = z[lo:i].min(), z[i + 1:hi].min()
            if z[i] - max(zl, zr) >= 0.5:
                znew = 0.5 * (zl + zr); fixes.append((float(tf.c + (c[i] + .5) * C), float(tf.f - (r[i] + .5) * C), float(z[i]), float(znew)))
                z[i] = znew; n_fix += 1
    if z[0] < z[-1]: z = z[::-1]; r = r[::-1]; c = c[::-1]           # order the line from its high end to its low end
    z = np.minimum.accumulate(z)                                      # 4. stream burn: non increasing downstream (no pits)
    out[r, c] = np.minimum(out[r, c], z)
# 5. residual channel pits (line junctions, inconsistent line directions): a channel cell with no lower 8-neighbour is raised
#    to its lowest neighbour + 1 mm, iterated; nodata (domain edge) counts as lower, so the outlet is never touched
ch = (dem - out > 0.01) & act
for it in range(500):
    zi = np.where(act, out, -np.inf); low = np.full_like(out, np.inf)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr or dc: low = np.minimum(low, np.roll(np.roll(zi, dr, 0), dc, 1))
    pit = ch & (low >= out)
    if not pit.any(): break
    out[pit] = low[pit] + 0.001
print(f"channel pit raising: {it} iterations, remaining channel pits {int(pit.sum())}")
out[~act] = -9999
hdr = f"ncols {NC}\nnrows {NR}\nxllcorner {tf.c}\nyllcorner {tf.f - NR * C}\ncellsize {C}\nNODATA_value -9999"
np.savetxt(f"{D}/../inputs/dem10_cond.asc", out, fmt="%.3f", header=hdr, comments="")
with rasterio.open(f"{D}/dem10_cond.tif", "w", driver="GTiff", height=NR, width=NC, count=1, dtype="float32", crs="EPSG:2169", transform=tf, nodata=-9999, compress="deflate") as d: d.write(out.astype("float32"), 1)
pd.DataFrame(fixes, columns=["x", "y", "z_before", "z_after"]).to_csv(f"{D}/crossing_fixes.csv", index=False)
d = (dem - out)[act]
print(f"lines {len(g)} ({g.src.value_counts().to_dict()}); channel cells {n_burn}; crossings lowered {n_fix}; "
      f"mean lowering on channel cells {d[d > 0].mean():.2f} m, max {d.max():.1f} m; cells lowered {(d > 0).sum()}")
