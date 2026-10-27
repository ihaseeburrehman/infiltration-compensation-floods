"""Every score derived number the Alzette text quotes, from data/scores.csv and data/hydrographs.csv (station section = width
"full"; "valley" and "surveyed" are the sensitivities). usage: python numbers_dump.py <paper_v4 dir>"""
import sys, numpy as np, pandas as pd
P = sys.argv[1].rstrip("/"); pd.set_option("display.width", 250); pd.set_option("display.max_rows", 1000)
s = pd.read_csv(P + "/data/scores.csv"); hy = pd.read_csv(P + "/data/hydrographs.csv"); ob = pd.read_csv(P + "/data/obs_utc.csv")
F = ["F_graphcast", "F_aifs", "F_fuxi", "F_wrf12noda", "F_wrf12"]; VM = ["VM_graphcast", "VM_aifs", "VM_fuxi", "VM_wrf12noda", "VM_wrf12"]
def means(width, inf, cases, excl=None):
    d = s[(s.width == width) & (s.infiltration == inf) & (s.case.isin(cases))]
    if excl: d = d[d.gauge != excl]
    return d.groupby(["case", "model"])[["peak_ratio", "beta", "nse", "kge", "dt_peak"]].mean()
def rng(df, col, nd=2): return f"{df[col].min():.{nd}f} to {df[col].max():.{nd}f}"
for width in ("full", "valley", "surveyed"):
    print(f"\n######## width = {width}")
    for inf in ("none", "GA"):
        print(f"--- {inf}: gauge means")
        for name, cs in (("radar", ["R_1km_5min"]), ("forecasts", F), ("rescaled", VM)):
            m = means(width, inf, cs); print(f"  {name:10s} peak {rng(m,'peak_ratio')}  beta {rng(m,'beta')}  NSE {rng(m,'nse')}  KGE {rng(m,'kge')}  dt {rng(m,'dt_peak',1)}")
        m = means(width, inf, ["R_1km_5min"] + F + VM); print(m.round(2).to_string())
        print("  excluding Livange:"); 
        for name, cs in (("radar", ["R_1km_5min"]), ("forecasts", F), ("rescaled", VM)):
            m = means(width, inf, cs, "Livange"); print(f"  {name:10s} peak {rng(m,'peak_ratio')}  NSE {rng(m,'nse')}")
print("\n######## per gauge, station section")
g = s[(s.width == "full") & (s.infiltration.isin(["none", "GA"])) & (s.case.isin(["R_1km_5min"] + F + VM))]
print(g.pivot_table(index=["infiltration", "case", "model"], columns="gauge", values=["peak_ratio", "nse", "dt_peak"]).round(2).to_string())
print("\n######## per gauge ranges over the five forecasts as issued (none)")
for gg in ["Pfaffenthal", "Hesperange", "Steinsel", "Livange"]:
    d = g[(g.gauge == gg) & (g.case.isin(F)) & (g.infiltration == "none")]; r = g[(g.gauge == gg) & (g.case == "R_1km_5min")]
    print(f"  {gg:12s} forecasts peak {rng(d,'peak_ratio')} NSE {rng(d,'nse')} | radar none peak {rng(r[r.infiltration=='none'],'peak_ratio')} NSE {rng(r[r.infiltration=='none'],'nse')} | radar GA peak {rng(r[r.infiltration=='GA'],'peak_ratio')} NSE {rng(r[r.infiltration=='GA'],'nse')}")
print("\n######## model agreement (station section): |LF-HP| mean peak ratio per case/inf, beta diff, per gauge peak % diff")
m = s[(s.width == "full")].groupby(["case", "infiltration", "model"])[["peak_ratio", "beta", "nse"]].mean().unstack("model")
d = (m["peak_ratio"]["HiPIMS"] - m["peak_ratio"]["LISFLOOD-FP"]).dropna(); print("  max |d peak| over pairs:", round(d.abs().max(), 3), " n pairs", len(d))
db = (m["beta"]["HiPIMS"] - m["beta"]["LISFLOOD-FP"]).dropna(); print(f"  beta diff HP-LF: {db.min():.3f} to {db.max():.3f}")
pg = s[(s.width == "full") & (s.infiltration.isin(["none", "GA"])) & (s.case.isin(["R_1km_5min"] + F + VM))].pivot_table(index=["case", "infiltration", "gauge"], columns="model", values=["peak", "nse"])
pct = (pg["peak"]["HiPIMS"] / pg["peak"]["LISFLOOD-FP"] - 1) * 100; dn = pg["nse"]["HiPIMS"] - pg["nse"]["LISFLOOD-FP"]
for gg in ["Livange", "Pfaffenthal", "Hesperange", "Steinsel"]:
    x = pct.xs(gg, level="gauge"); y = dn.xs(gg, level="gauge"); yg = y.xs("GA", level="infiltration")
    print(f"  {gg:12s} HP peak vs LF: {x.min():+.0f} to {x.max():+.0f} %;  NSE diff GA {yg.min():+.2f} to {yg.max():+.2f}")
print("\n######## sensitivity sets (station section)")
for inf in ("GA_ks05", "GA_ks2", "GA_sm"):
    d = s[(s.width == "full") & (s.infiltration == inf)].groupby(["case", "model"])[["peak_ratio", "beta", "nse", "kge"]].mean(); print(inf); print(d.round(2).to_string())
    for cs in (["R_1km_5min"], ["F_graphcast"], ["F_aifs", "F_fuxi"]):
        dd = d.loc[[c for c in cs if c in d.index.get_level_values(0)]]
        if len(dd): print(f"   {cs} peak {rng(dd,'peak_ratio')} NSE {rng(dd,'nse')} KGE {rng(dd,'kge')} beta {rng(dd,'beta')}")
    sm = s[(s.width == "full") & (s.infiltration == inf) & (s.case == "R_1km_5min")]; print("   radar per gauge peak:", sm.groupby("gauge").peak_ratio.agg(["min", "max"]).round(2).to_dict("index"))
d = s[(s.width == "full") & (s.case == "R_1km_6h")].groupby(["infiltration", "model"])[["peak_ratio", "nse", "dt_peak", "beta"]].mean(); print("R_1km_6h"); print(d.round(2).to_string())
print("\n######## surveyed: fraction of the station-section peak excluded")
pk = s[(s.infiltration.isin(["none", "GA"])) & (s.case.isin(["R_1km_5min"] + F + VM))].pivot_table(index=["case", "infiltration", "model", "gauge"], columns="width", values="peak")
ex = (1 - pk["surveyed"] / pk["full"]) * 100
for gg in ["Pfaffenthal", "Hesperange", "Steinsel", "Livange"]: x = ex.xs(gg, level="gauge"); print(f"  {gg:12s} excluded {x.min():.0f} to {x.max():.0f} %")
print("\n######## event volume with the antecedent flow (hour 1 value) subtracted, radar GA, station section, Steinsel and all gauges")
for gg in ["Pfaffenthal", "Hesperange", "Steinsel", "Livange"]:
    o = ob[ob.gauge == gg].set_index("h").Q; o = o[(o.index >= 1) & (o.index <= 96)]
    for m in ("LISFLOOD-FP", "HiPIMS"):
        q = hy[(hy.case == "R_1km_5min") & (hy.model == m) & (hy.infiltration == "GA") & (hy.gauge == gg) & (hy.width == "full")].set_index("h").Q; q = q[(q.index >= 1) & (q.index <= 96)]
        print(f"  {gg:12s} {m:12s} beta {q.sum()/o.sum():.2f}  event volume ratio (minus h1 flow) {(q - q.iloc[0]).sum()/(o - o.iloc[0]).sum():.2f}  obs h1 {o.iloc[0]:.1f} m3/s, antecedent share of obs volume {100*o.iloc[0]*96/o.sum():.0f} %")
