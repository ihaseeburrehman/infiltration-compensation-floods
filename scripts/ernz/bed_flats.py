"""Flat reaches in the conditioned Ernz bed, without the line files: burnt cells are those where the conditioned DEM lies below
the block mean DEM; a burnt cell is 'flat' when no burnt 8-neighbour lies more than 1 mm below it (no downstream fall).
Connected flat cells form a reach; reports count, longest and total length, and the largest reaches with their location and
elevation, plus the bed slope distribution along burnt cells. usage: python bed_flats.py"""
import numpy as np, rasterio, json
from scipy import ndimage
with rasterio.open("domain/dem10_cond.tif") as s: dem = s.read(1).astype(float); tf = s.transform
with rasterio.open("domain/dem10_mean.tif") as s: raw = s.read(1).astype(float)
act = dem > -9000; burnt = act & (dem < raw - 0.01)
print("burnt cells", burnt.sum(), f"({burnt.sum()*10/1000:.0f} km of channel at one cell per 10 m); mean lowering {np.mean(raw[burnt]-dem[burnt]):.2f} m")
# lowest burnt neighbour
big = 1e9; lowest = np.full(dem.shape, big)
for dr in (-1, 0, 1):
    for dc in (-1, 0, 1):
        if dr == dc == 0: continue
        sh = np.full(dem.shape, big); src = np.where(burnt, dem, big)
        r0, r1 = max(dr, 0), dem.shape[0] + min(dr, 0); c0, c1 = max(dc, 0), dem.shape[1] + min(dc, 0)
        sh[r0:r1, c0:c1] = src[r0 - dr:r1 - dr, c0 - dc:c1 - dc]; lowest = np.minimum(lowest, sh)
drop = np.where(burnt, dem - lowest, np.nan)   # fall to the lowest burnt neighbour (negative where the cell is itself the lowest = pit)
flat = burnt & (drop < 0.001) & (lowest < big)
pit = burnt & (drop < -0.001)
print(f"burnt cells with no fall to a burnt neighbour (flat): {flat.sum()} ({100*flat.sum()/burnt.sum():.1f} %); local pits (all burnt neighbours higher): {pit.sum()}")
lab, n = ndimage.label(flat, structure=np.ones((3, 3))); sizes = ndimage.sum(flat, lab, range(1, n + 1)) * 10
print(f"flat reaches: {n}; longest {sizes.max():.0f} m; total {sizes.sum()/1000:.1f} km; reaches >=100 m: {(sizes>=100).sum()}, >=300 m: {(sizes>=300).sum()}, >=1000 m: {(sizes>=1000).sum()}")
T = json.load(open("scripts/transects_ernz.json"))
for i in np.argsort(sizes)[::-1][:8]:
    rr, cc = np.where(lab == i + 1); x = tf.c + (cc.mean() + 0.5) * tf.a; y = tf.f - (rr.mean() + 0.5) * tf.a
    near = min(T, key=lambda g: (T[g]["x"] - x) ** 2 + (T[g]["y"] - y) ** 2); dist = np.hypot(T[near]["x"] - x, T[near]["y"] - y)
    print(f"  reach {sizes[i]:5.0f} m at ({x:.0f}, {y:.0f}), bed {dem[rr, cc].mean():.1f} m, nearest transect {near} {dist/1000:.1f} km away")
d = drop[burnt & (drop >= 0) & (lowest < big)]
print("fall per 10 m cell along burnt cells: median %.3f m, p10 %.3f, p90 %.3f; share with fall < 1 cm: %.1f %%, < 5 cm: %.1f %%" % (np.median(d), np.percentile(d, 10), np.percentile(d, 90), 100 * (d < 0.01).mean(), 100 * (d < 0.05).mean()))
