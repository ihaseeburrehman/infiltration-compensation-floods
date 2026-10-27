"""Sentinel-1 (EMSN139, 15 July 05:51 UTC) against the model near each gauge. For windows of radius R around the gauge:
SAR flooded area, model flooded area at hour 36 (depth>0.1 m, start-wet cells excluded, EMSN139 areas only), and the lateral
distance of the flooded cells from the channel (start-wet cells), as max and 95th percentile. Also checks the SAR mask
against the channel. Runs on the HPC (hipims/sar_window.py); writes sar_width.csv."""
import sys, json, csv, numpy as np
sys.path.insert(0, "/ec/res4/scratch/lux0804/hipims"); sys.path.insert(0, "/ec/res4/scratch/lux0804/radar_scale_2026/code")
from xs_fullvalley import readers, NR, X0, YTOP, CS
from scipy import ndimage
S = "/ec/res4/scratch/lux0804/"; M = np.load(S + "radar_scale_2026/code/sar_masks.npz"); T = json.load(open("transects_fullvalley.json"))
ROWS = list(range(NR))
def grid(rd, t):
    d = rd(t, "h", ROWS); return np.array([d[r] for r in ROWS], dtype=np.float32)
runs = {"LF_none": S + "radar_scale_2026/runs/L84_R_1km_5min_noinf_dp_1g", "LF_GA": S + "radar_scale_2026/runs/L84_R_1km_5min_ga3b_dp"}
G = {}
for k, run in runs.items():
    ts, rd = readers(run, "lf"); G[k] = grid(rd, 36)
    if "perm" not in G: G["perm"] = grid(rd, 0) > 0.1
perm = G["perm"]; sar = M["sar15"]; aoi = M["aoi01"] | M["aoi02"]
print("channel cells seen by SAR: %d of %d (%.0f %%); SAR cells %d, within 300 m of the river %d" % ((sar & perm).sum(), perm.sum(), 100 * (sar & perm).sum() / perm.sum(), sar.sum(), (sar & M["riv300"]).sum()))
dist = ndimage.distance_transform_edt(~perm) * CS
yy, xx = np.mgrid[0:NR, 0:perm.shape[1]]; X = X0 + (xx + 0.5) * CS; Y = YTOP - (yy + 0.5) * CS
def st(m):
    return (round(m.sum() / 100, 1), round(float(dist[m].max())) if m.any() else "", round(float(np.percentile(dist[m], 95))) if m.any() else "")
rows = []
for g, v in T.items():
    for R in (300, 500, 1000):
        W = ((X - v["x"]) ** 2 + (Y - v["y"]) ** 2) <= R * R; cov = (W & aoi).sum() / W.sum()
        r = dict(gauge=g, radius_m=R, aoi_cover=round(cov, 3)); r["sar_ha"], r["sar_dmax"], r["sar_d95"] = st(sar & W & ~perm)
        for k in runs: r[k + "_ha"], r[k + "_dmax"], r[k + "_d95"] = st((G[k] > 0.1) & W & ~perm & aoi)
        rows.append(r); print(r)
with open("sar_width.csv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
