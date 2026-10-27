"""Simulated water level at the Walferdange gauge (77256.0, 81572.0 LUREF) for all 36 runs.
Level = highest water surface (bed + depth) among wet cells (depth > 1 cm) whose centre lies within 15 m of the gauge,
because the gauge sits at the channel edge and its own cell can be dry at low flow. Writes stage_walferdange.csv."""
import sys, os, csv, numpy as np, multiprocessing as mp
sys.path.insert(0, "/ec/res4/scratch/lux0804/hipims"); sys.path.insert(0, "/ec/res4/scratch/lux0804/radar_scale_2026/code")
import seg_extract as se
from xs_fullvalley import readers, X0, YTOP, CS
GX, GY, RAD = 77256.0, 81572.0, 15.0
S = "/ec/res4/scratch/lux0804/"
C = ["R_1km_5min", "R_1km_6h", "F_graphcast", "F_aifs", "F_wrf1p3", "F_wrf12", "F_wrf12noda", "VM_graphcast", "VM_aifs", "VM_wrf1p3"]
RUNS = []
for c in C:
    RUNS += [(c, "LISFLOOD-FP", "none", "lf", S + f"radar_scale_2026/runs/L84_{c}_noinf_dp_1g"), (c, "LISFLOOD-FP", "GA", "lf", S + f"radar_scale_2026/runs/L84_{c}_ga3b_dp"),
             (c, "HiPIMS", "none", "hp", S + f"hipims/runs/{c}_mg4"), (c, "HiPIMS", "GA", "hp", S + f"hipims/runs/{c}_inf_mg4")]
for c in ("R_1km_5min", "F_graphcast"):
    for k in ("ks05", "ks2"):
        RUNS += [(c, "LISFLOOD-FP", f"GA_{k}", "lf", S + f"radar_scale_2026/runs/L84_{c}_ga3b_{k}_dp"), (c, "HiPIMS", f"GA_{k}", "hp", S + f"hipims/runs/{c}_inf_{k}_mg4")]
r0, c0 = int((YTOP - GY) // CS), int((GX - X0) // CS)
CELLS = [(r, c) for r in range(r0 - 3, r0 + 4) for c in range(c0 - 3, c0 + 4)
         if np.hypot(X0 + (c + .5) * CS - GX, YTOP - (r + .5) * CS - GY) <= RAD]
ROWS = sorted({r for r, _ in CELLS})
Z = se.read_rows(S + "radar_scale_2026/inputs/Alzette_sub_basin_10m_bridge_burn.asc", ROWS)
def one(job):
    case, model, inf, kind, run = job
    ts, rd = readers(run, kind); out = []
    for t in ts:
        h = rd(t, "h", ROWS)
        wse = [Z[r][c] + h[r][c] for r, c in CELLS if h[r][c] > 0.01]
        out.append((case, model, inf, t if kind == "lf" else t // 3600, round(max(wse), 3) if wse else ""))
    return out
if __name__ == "__main__":
    print("cells", len(CELLS), "bed", [round(float(Z[r][c]), 2) for r, c in CELLS])
    with mp.Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", 4))) as p: res = p.map(one, RUNS)
    with open(S + "hipims/stage_walferdange.csv", "w") as f:
        w = csv.writer(f); w.writerow(["case", "model", "infiltration", "h", "wse"]); [w.writerows(r) for r in res]
    print("done", sum(len(r) for r in res))
