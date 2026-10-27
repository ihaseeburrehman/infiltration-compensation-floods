"""Water balance (million m3) of the LISFLOOD-FP runs from the .mass files (massint 600 s).
Qout is a rate (m3/s) integrated over each 600 s interval. The last column, Rain-(Inf+Evap), is cumulative and excludes
the dynamic rainfall, so it equals minus the infiltrated volume. The applied rainfall is taken from the closed balance of
the matching no-infiltration run (same rainfall file); the Green-Ampt residual then checks closure."""
import glob, numpy as np, csv
S = "/ec/res4/scratch/lux0804/radar_scale_2026/runs/"
C = ["R_1km_5min", "R_1km_6h", "F_graphcast", "F_aifs", "F_fuxi", "F_wrf12noda", "F_wrf12", "VM_graphcast", "VM_aifs", "VM_fuxi", "VM_wrf12noda", "VM_wrf12", "F_wrf1p3", "VM_wrf1p3"]
def wb(run):
    f = glob.glob(S + run + "/results/*.mass") or glob.glob(S + run.replace("_dp", "_dpn") + "/results/*.mass")
    a = np.loadtxt(f[0], skiprows=1)
    t = a[:, 0]; dt = np.diff(np.r_[0, t])
    return dict(v0=a[0, 5], v1=a[-1, 5], qout=(a[:, 8] * dt).sum(), inf=-a[-1, 11])
rows = []
for c in C:
    if not (glob.glob(S + f"L84_{c}_noinf_dp*_1g/results/*.mass") and glob.glob(S + f"L84_{c}_ga3b_dp*/results/*.mass")): print("missing", c); continue
    n, g = wb(f"L84_{c}_noinf_dp_1g"), wb(f"L84_{c}_ga3b_dp")
    rain = n["v1"] + n["qout"] - n["v0"]; inf = g["inf"]
    runs = [("none", n, 0.0), ("GA", g, inf)]
    if glob.glob(S + f"L84_{c}_ga3b_sm_dpn/results/*.mass"):
        x = wb(f"L84_{c}_ga3b_sm_dpn"); runs.append(("GA_sm", x, x["inf"]))
    for k in ("ks05", "ks2"):
        if glob.glob(S + f"L84_{c}_ga3b_{k}_dp/results/*.mass"):
            x = wb(f"L84_{c}_ga3b_{k}_dp"); runs.append((f"GA_{k}", x, x["inf"]))
    for lab, w, i in runs:
        rows.append(dict(case=c, infiltration=lab, rain=rain / 1e6, infiltrated=i / 1e6, outflow=w["qout"] / 1e6, storage0=w["v0"] / 1e6,
                         storage1=w["v1"] / 1e6, residual=(w["v0"] + rain - i - w["qout"] - w["v1"]) / 1e6))
w = csv.DictWriter(open("/ec/res4/scratch/lux0804/hipims/water_balance.csv", "w"), fieldnames=list(rows[0])); w.writeheader()
for r in rows: w.writerow({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}); print(*[f"{v:.2f}" if isinstance(v, float) else f"{v:12s}" for v in r.values()])
