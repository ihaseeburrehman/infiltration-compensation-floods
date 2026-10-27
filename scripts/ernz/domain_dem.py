"""Ernz Noire 2018 case: domain polygon from the AGE 'Bassins versants' hierarchy and 10 m DEM from the 0.5 m LiDAR.
Domain = all sub-basins under EZGID 1021 (Ernz Noire to the Sure at Grundhof, 101.8 km2, includes the Hallerbach) and
EZGID 951 (Ernz Blanche to the Sure at Reisdorf, 101.4 km2, Larochette gauge of CAMELS-LUX), the two Ernz valleys.
DEM: block mean and block min of the 0.5 m MNT (20 x 20 cells) on a 10 m EPSG:2169 grid aligned to 10 m, bbox + 300 m.
usage: python domain_dem.py <bassins_versants.geojson> <MNT_Lidar2024.tif> <outdir>"""
import sys, collections, numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.windows import Window
from rasterio.transform import from_origin
from rasterio.features import rasterize
GJ, MNT, O = sys.argv[1:4]
ROOTS = (1021, 951); BUF = 300.0; C = 10.0; K = 20
g = gpd.read_file(GJ).to_crs(2169)
kids = collections.defaultdict(list)
for e, u in zip(g.EZGID, g.UEZGID):
    if u == u: kids[int(u)].append(int(e))
def tree(root):
    sel, st = [root], [root]
    while st:
        c = st.pop(); sel += kids.get(c, []); st += kids.get(c, [])
    return sel
parts = []
for root, name in zip(ROOTS, ("ernz_noire", "ernz_blanche")):
    p = g[g.EZGID.isin(tree(root))].dissolve(); p["name"] = name; parts.append(p)
    p.to_file(f"{O}/{name}.geojson", driver="GeoJSON"); print(name, "km2 %.1f" % (p.area.iloc[0] / 1e6))
dom = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True)).dissolve(); dom["name"] = "Ernz Noire and Ernz Blanche to the Sure (AGE EZGID 1021 + 951)"
dom.to_file(f"{O}/domain.geojson", driver="GeoJSON")
hal, st = [1013], [1013]
while st:
    c = st.pop(); hal += kids.get(c, []); st += kids.get(c, [])
g[g.EZGID.isin(hal)].dissolve().to_file(f"{O}/hallerbach.geojson", driver="GeoJSON")
poly = dom.geometry.iloc[0]; print("domain area km2 %.1f" % (poly.area / 1e6), "bounds", poly.bounds, flush=True)
x0 = np.floor((poly.bounds[0] - BUF) / C) * C; x1 = np.ceil((poly.bounds[2] + BUF) / C) * C
y0 = np.floor((poly.bounds[1] - BUF) / C) * C; y1 = np.ceil((poly.bounds[3] + BUF) / C) * C
NC, NR = int((x1 - x0) / C), int((y1 - y0) / C); tf = from_origin(x0, y1, C, C)
print("grid", NR, NC, "cells", NR * NC, "origin", x0, y0, flush=True)
mask = rasterize([(poly, 1)], out_shape=(NR, NC), transform=tf, fill=0, dtype="uint8").astype(bool)
mean = np.full((NR, NC), np.nan, "float32"); mn = np.full((NR, NC), np.nan, "float32")
with rasterio.open(MNT) as s:
    assert abs(s.res[0] - 0.5) < 1e-6
    col0 = int(round((x0 - s.bounds.left) / 0.5)); row0 = int(round((s.bounds.top - y1) / 0.5))
    STEP = 100                                                   # 10 m rows per stripe = 2000 fine rows
    for r in range(0, NR, STEP):
        n = min(STEP, NR - r)
        a = s.read(1, window=Window(col0, row0 + r * K, NC * K, n * K)).astype("float64")
        a[a <= -9000] = np.nan
        b = a.reshape(n, K, NC, K).transpose(0, 2, 1, 3).reshape(n, NC, K * K)
        with np.errstate(all="ignore"):
            mean[r:r + n] = np.nanmean(b, axis=2); mn[r:r + n] = np.nanmin(b, axis=2)
        if r % 500 == 0: print("stripe", r, "/", NR, flush=True)
for name, arr in (("dem10_mean", mean), ("dem10_min", mn)):
    out = np.where(mask & np.isfinite(arr), arr, -9999).astype("float32")
    with rasterio.open(f"{O}/{name}.tif", "w", driver="GTiff", height=NR, width=NC, count=1, dtype="float32", crs="EPSG:2169",
                       transform=tf, nodata=-9999, compress="deflate") as d: d.write(out, 1)
np.save(f"{O}/domain_mask.npy", mask)
v = mean[mask & np.isfinite(mean)]
print("done; valid cells", int((mask & np.isfinite(mean)).sum()), "elev min/max %.1f %.1f" % (v.min(), v.max()),
      "mask cells without LiDAR", int((mask & ~np.isfinite(mean)).sum()))
