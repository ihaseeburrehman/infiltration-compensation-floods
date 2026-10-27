"""Time to peak and peak at Larochette and Mullerthal: default roughness against the vegetated classes set to n = 0.05 (tree
0.10 and bush 0.07), radar forcing, without and with Green-Ampt. Rain peak hour from rain_totals (hour 6 = 00:00 to 01:00 UTC)."""
import pandas as pd
o = pd.read_csv("results/obs_larochette_utc.csv"); ho = int(o.loc[o.Q.idxmax(), "h"]); print(f"observed Larochette: peak {o.Q.max():.1f} m3/s at hour {ho} (rain peak hour 6, so {ho-6} h after the rain)")
for inf in ("noinf", "ga"):
    for tag, f in (("default n", f"results/hydro_R_1km_5min_{inf}.csv"), ("veg n=0.05", f"results/hydro_R_1km_5min_veg05_{inf}.csv")):
        h = pd.read_csv(f); L = h.Larochette; M = h.Mullerthal; tp = int(h.loc[L.idxmax(), "hour"])
        nse = 1 - ((L.values[:len(o)] - o.Q.values[:len(L)]) ** 2).sum() / ((o.Q - o.Q.mean()) ** 2).sum()
        print(f"{inf:5s} {tag:11s}: Larochette peak {L.max():6.1f} at hour {tp} ({tp-6} h after rain, {tp-ho:+d} h vs observed), Q at hour 9 {L[h.hour==9].iloc[0]:5.1f}, at hour 16 {L[h.hour==16].iloc[0]:5.1f}, NSE {nse:6.2f}; "
              f"Mullerthal peak {M.max():6.1f} at hour {int(h.loc[M.idxmax(),'hour'])}")
