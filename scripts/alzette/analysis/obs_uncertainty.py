"""Robustness of the ranking to the discharge observations, which at the 2021 peaks rest on rating curves extrapolated beyond
the gaugings (AGE rating curves not available). The observed series of every gauge is rescaled by a factor f and every run
is rescored (station section, scores as in score_all.py). Two error shapes: "uniform" multiplies the whole series by f,
"high" multiplies only the part above half the observed peak (Q' = Qb + f (Q - Qb) for Q > Qb, Qb = 0.5 peak), where the
rating curve is extrapolated. usage: python obs_uncertainty.py <paper_v4 dir>   ->  data/obs_uncertainty.csv"""
import sys, numpy as np, pandas as pd
sys.argv += []; P = sys.argv[1].rstrip("/")
hy = pd.read_csv(P + "/data/hydrographs.csv"); ob = pd.read_csv(P + "/data/obs_utc.csv")
hy = hy[(hy.width == "full") & (hy.infiltration.isin(["none", "GA"]))]
F = ["F_graphcast", "F_aifs", "F_fuxi", "F_wrf12noda", "F_wrf12"]; hy = hy[hy.case.isin(["R_1km_5min"] + F)]
def score(qs, qo):                                   # NSE and peak ratio as in score_all.score
    h = [k for k in range(1, 97) if k in qs and not np.isnan(qo.get(k, np.nan))]
    s = np.array([qs[k] for k in h]); o = np.array([qo[k] for k in h])
    return 1 - ((s - o) ** 2).sum() / ((o - o.mean()) ** 2).sum(), s.max() / o.max()
rows = []
for shape in ("uniform", "high"):
    for f in (0.7, 0.8, 0.9, 1.0, 1.1, 1.15, 1.2, 1.3):
        obs = {}
        for g, d in ob.groupby("gauge"):
            q = d.set_index("h").Q.astype(float); qb = 0.5 * q[(q.index >= 1) & (q.index <= 96)].max()
            obs[g] = (q * f if shape == "uniform" else q.where(q <= qb, qb + f * (q - qb))).to_dict()
        for (c, m, inf, g), d in hy.groupby(["case", "model", "infiltration", "gauge"]):
            nse, pr = score(dict(zip(d.h, d.Q)), obs[g]); rows.append(dict(shape=shape, f=f, case=c, model=m, infiltration=inf, gauge=g, nse=nse, peak_ratio=pr))
r = pd.DataFrame(rows); r.round(4).to_csv(P + "/data/obs_uncertainty.csv", index=False)
m = r.groupby(["shape", "f", "infiltration", "case", "model"])[["nse", "peak_ratio"]].mean().reset_index()
m["grp"] = np.where(m.case == "R_1km_5min", "radar", "forecasts")
t = m.groupby(["shape", "f", "infiltration", "grp"]).agg(nse_min=("nse", "min"), nse_max=("nse", "max"), pk_min=("peak_ratio", "min"), pk_max=("peak_ratio", "max")).round(2)
print(t.to_string())
print("\nradar NSE above every forecast (min radar > max forecasts), per model:")
for (sh, f, inf, mo), d in m.groupby(["shape", "f", "infiltration", "model"]):
    rad = d[d.grp == "radar"].nse.iloc[0]; fc = d[d.grp == "forecasts"].nse
    print(f"  {sh:7s} f={f:.1f} {inf:4s} {mo:12s} radar {rad:5.2f}  forecasts {fc.min():5.2f}..{fc.max():5.2f}  radar best: {rad > fc.max()}  radar worst: {rad < fc.min()}")
