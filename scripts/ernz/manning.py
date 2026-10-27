"""Manning map for the Ernz grid from the Luxembourg Land Cover 2018 geodatabase (LC_2018.gdb, layer LC2018_Luxembourg_country,
field LABEL_en), following the user's creating_manning_lux_10m.py (same lookup, bbox read, rasterize, DEM mask) plus the
fill_manning_gaps.py step (unmapped cells inside the DEM get 0.035).
usage: python manning.py <dem10 tif or asc> <outdir> [gdb path]"""
import sys, numpy as np, geopandas as gpd, rasterio
from rasterio.features import rasterize
DEM, O = sys.argv[1:3]; GDB = sys.argv[3] if len(sys.argv) > 3 else "/Users/haseeb.rehman/Downloads/landcover2018-vector/LC_2018.gdb"
N = {"Bare soil": 0.02, "Building": 0.20, "Bush": 0.07, "Other constructed area": 0.016, "Permanent herbaceous vegetation": 0.035,
     "Seasonal herbaceous vegetation": 0.04, "Tree": 0.1, "Water": 0.03, "Vineyard": 0.045}; DEFAULT = 0.035
with rasterio.open(DEM) as s: dem = s.read(1); tf = s.transform; NR, NC = dem.shape; nd = s.nodata if s.nodata is not None else -9999
act = (dem != nd) & (dem > -9000)
bbox = (tf.c, tf.f + NR * tf.e, tf.c + NC * tf.a, tf.f)
g = gpd.read_file(GDB, layer="LC2018_Luxembourg_country", bbox=bbox, engine="pyogrio")
print("land cover polygons in bbox", len(g), "classes", sorted(g.LABEL_en.unique()))
g["n"] = g.LABEL_en.map(N).fillna(DEFAULT)
n = rasterize(((geom, float(v)) for geom, v in zip(g.geometry, g.n) if geom is not None), out_shape=(NR, NC), transform=tf, fill=-9999, dtype="float32")
gaps = act & (n == -9999); print("unmapped active cells filled with %.3f: %d (%.2f %%)" % (DEFAULT, gaps.sum(), 100 * gaps.sum() / act.sum())); n[gaps] = DEFAULT
cls = rasterize(((geom, i + 1) for i, (geom) in enumerate(g.geometry) if geom is not None), out_shape=(NR, NC), transform=tf, fill=0, dtype="int32")
for lab, v in N.items(): print("  %-32s n=%.3f  %5.1f %% of domain" % (lab, v, 100 * ((n == v) & act).sum() / act.sum()))
hdr = f"ncols {NC}\nnrows {NR}\nxllcorner {tf.c}\nyllcorner {tf.f + NR * tf.e}\ncellsize {tf.a}\nNODATA_value -9999"
np.savetxt(f"{O}/manning.n.ascii", np.where(act, n, -9999), fmt="%.4f", header=hdr, comments="")
print("written", f"{O}/manning.n.ascii", "mean n on domain %.4f" % n[act].mean())
