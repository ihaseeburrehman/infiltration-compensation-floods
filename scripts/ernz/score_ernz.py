"""Scores at Larochette for the Ernz runs: hourly means of the CAMELS-LUX 15 min record (UTC = stamp - 1 h, interval end) against the
instantaneous hourly simulated discharge, as score_all.py does for the Alzette (NSE, KGE', peak ratio, beta, peak timing).
usage: python score_ernz.py <start YYYY-mm-ddTHH UTC> <hours> <out csv> <hydro csv ...>   (hydro csv named hydro_<case>_<inf>.csv)"""
import sys, os, re, numpy as np, pandas as pd
T0 = pd.Timestamp(sys.argv[1]); H = int(sys.argv[2]); OUT = sys.argv[3]; FILES = sys.argv[4:]
obs = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"../data/camels_larochette_15min_{T0.year}.csv"), parse_dates=["Date"])
obs["utc"] = obs.Date - pd.Timedelta(hours=1)
# hourly mean of the four 15 min values whose interval ends in (h, h+1]: the mean flow of hour [h, h+1); compared with the simulated value at h+1? No:
# as in the Alzette, the observation of model hour h is the mean of the native records within the hour centred on h: stamps in (h-30min, h+30min].
oh = []
for k in range(H + 1):
    t = T0 + pd.Timedelta(hours=k); w = obs[(obs.utc > t - pd.Timedelta(minutes=30)) & (obs.utc <= t + pd.Timedelta(minutes=30))]
    oh.append(w.Q.mean())
o = np.array(oh); ho = int(np.nanargmax(o))
rows = []
for f in FILES:
    m = re.match(r"hydro_(.+)_(noinf|ga)\.csv", os.path.basename(f)); case, inf = m.groups()
    h = pd.read_csv(f); s = np.interp(np.arange(H + 1), h.hour, h.Larochette)
    r = np.corrcoef(s, o)[0, 1]; beta = s.mean() / o.mean(); gam = (s.std() / s.mean()) / (o.std() / o.mean())
    nse = 1 - ((s - o) ** 2).sum() / ((o - o.mean()) ** 2).sum(); kge = 1 - np.sqrt((r - 1) ** 2 + (beta - 1) ** 2 + (gam - 1) ** 2)
    rows.append(dict(case=case, infiltration=inf, peak=s.max(), peak_obs=o.max(), peak_ratio=s.max() / o.max(), t_peak=int(s.argmax()), dt_peak=int(s.argmax()) - ho,
                     beta=beta, r=r, gamma=gam, nse=nse, kge=kge, vol_sim_Mm3=s.sum() * 3600 / 1e6, vol_obs_Mm3=o.sum() * 3600 / 1e6,
                     mullerthal_peak=h.Mullerthal.max(), grundhof_peak=h.Grundhof.max(), beaufort_peak=h.Beaufort.max()))
df = pd.DataFrame(rows).round(3); df.to_csv(OUT, index=False)
pd.set_option("display.width", 200); print(f"observed peak {o.max():.2f} m3/s at hour {ho} ({T0 + pd.Timedelta(hours=ho):%d %b %H} UTC), mean {o.mean():.2f}")
print(df[["case", "infiltration", "peak", "peak_ratio", "dt_peak", "beta", "nse", "kge", "mullerthal_peak", "grundhof_peak"]].to_string(index=False))
