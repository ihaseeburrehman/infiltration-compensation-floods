"""Extract 1.2 km windows around each gauge (Sentinel-1 EMSN139 mask of 15 July, EMSN139 areas, start-wet cells, LISFLOOD-FP
radar depth at hour 36 without and with infiltration) to one npz for the appendix map. Runs on the HPC (hipims/)."""
import sys, json, numpy as np
sys.path.insert(0, "/ec/res4/scratch/lux0804/hipims"); sys.path.insert(0, "/ec/res4/scratch/lux0804/radar_scale_2026/code")
from xs_fullvalley import readers, NR, X0, YTOP, CS
S = "/ec/res4/scratch/lux0804/"; M = np.load(S + "radar_scale_2026/code/sar_masks.npz"); T = json.load(open("transects_fullvalley.json"))
ROWS = list(range(NR))
def grid(rd, t):
    d = rd(t, "h", ROWS); return np.array([d[r] for r in ROWS], dtype=np.float32)
runs = {"none": S + "radar_scale_2026/runs/L84_R_1km_5min_noinf_dp_1g", "ga": S + "radar_scale_2026/runs/L84_R_1km_5min_ga3b_dp"}
G = {}
for k, run in runs.items():
    ts, rd = readers(run, "lf"); G[k] = grid(rd, 36)
    if "perm" not in G: G["perm"] = grid(rd, 0) > 0.1
out = {}; H = 60   # half window in cells (600 m)
for g, v in T.items():
    c = int((v["x"] - X0) // CS); r = int((YTOP - v["y"]) // CS); sl = (slice(r - H, r + H), slice(c - H, c + H))
    out[g + "_sar"] = M["sar15"][sl]; out[g + "_aoi"] = (M["aoi01"] | M["aoi02"])[sl]; out[g + "_perm"] = G["perm"][sl]
    out[g + "_none"] = G["none"][sl]; out[g + "_ga"] = G["ga"][sl]; out[g + "_origin"] = np.array([X0 + (c - H) * CS, YTOP - (r - H) * CS])
np.savez_compressed("sar_windows.npz", **out); print("saved", list(out)[:6])
