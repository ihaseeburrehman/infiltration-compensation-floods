"""Simulated flood extent against the Sentinel-1 flood masks of EMSN139 (15 July 05:51 UTC -> model hour 36;
16 July -> hour 60), for both models. Wet = depth > thr; cells wet at the start (river) are excluded; scores within
the EMSN139 areas, either near the river (riv500, 500 m) or the whole area. Writes sar_scores.csv."""
import sys, os, csv, numpy as np, multiprocessing as mp
sys.path.insert(0, "/ec/res4/scratch/lux0804/hipims"); sys.path.insert(0, "/ec/res4/scratch/lux0804/radar_scale_2026/code")
from xs_fullvalley import readers, NR
S = "/ec/res4/scratch/lux0804/"; M = np.load(S + "radar_scale_2026/code/sar_masks.npz")
C = ["R_1km_5min", "R_1km_6h", "F_graphcast", "F_aifs", "F_fuxi", "F_wrf12noda", "F_wrf12", "VM_graphcast", "VM_aifs", "VM_fuxi", "VM_wrf12noda", "VM_wrf12", "F_wrf1p3", "VM_wrf1p3"]
RUNS = []
for c in C:
    dp = "_dpn" if "fuxi" in c else "_dp"
    RUNS += [(c, "LISFLOOD-FP", "none", "lf", S + f"radar_scale_2026/runs/L84_{c}_noinf{dp}_1g"), (c, "LISFLOOD-FP", "GA", "lf", S + f"radar_scale_2026/runs/L84_{c}_ga3b{dp}"),
             (c, "HiPIMS", "none", "hp", S + f"hipims/runs/{c}_mg4"), (c, "HiPIMS", "GA", "hp", S + f"hipims/runs/{c}_inf_mg4")]
for c in ("R_1km_5min", "F_graphcast", "F_aifs", "F_fuxi"):
    RUNS += [(c, "LISFLOOD-FP", "GA_sm", "lf", S + f"radar_scale_2026/runs/L84_{c}_ga3b_sm_dpn"), (c, "HiPIMS", "GA_sm", "hp", S + f"hipims/runs/{c}_inf_sm_mg4")]
ROWS = list(range(NR))
def grid(rd, t):
    d = rd(t, "h", ROWS); return np.array([d[r] for r in ROWS], dtype=np.float32)
def one(job):
    case, model, inf, kind, run = job
    ts, rd = readers(run, kind); k = 3600 if kind == "hp" else 1
    t0 = 0 if kind == "lf" else min(t for t in [0] + ts)
    perm = grid(rd, 0) > 0.1; out = []
    for h, key in ((36, "sar15"), (60, "sar16")):
        wd = grid(rd, h * k); O = M[key]
        for zone in ("riv500", "aoi"):
            Z = (M["aoi01"] | M["aoi02"]) & ~perm
            if zone == "riv500": Z &= M["riv500"]
            for thr in (0.1, 0.3):
                Sx = wd > thr; H = (Sx & O & Z).sum(); Mi = (~Sx & O & Z).sum(); F = (Sx & ~O & Z).sum()
                out.append(dict(case=case, model=model, infiltration=inf, time=key, zone=zone, thr=thr, hit_rate=round(H / (H + Mi), 4),
                                FAR=round(F / (H + F), 4), CSI=round(H / (H + Mi + F), 4), bias=round((H + F) / (H + Mi), 4),
                                sim_ha=(Sx & Z).sum() / 100, obs_ha=(O & Z).sum() / 100))
    print(case, model, inf, flush=True); return out
if __name__ == "__main__":
    with mp.Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", 4)) // 4) as p: res = p.map(one, RUNS)
    rows = [r for x in res for r in x]
    with open(S + "hipims/sar_scores.csv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("done", len(rows))
