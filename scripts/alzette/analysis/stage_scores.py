"""Stage at Walferdange in m above a local zero, for observations and simulations.
Observed stage = AGE reading above the gauge zero (datum 225.32 m). Simulated stage = simulated water surface minus the
lowest DEM bed among the cells within 15 m of the gauge (225.87 m), because the 10 m DEM does not resolve the gauge zero.
Scores: peak stage error (m), peak timing error (h), rise from h1 to the peak (m). usage: stage_scores.py <paper_v4 dir>"""
import sys, pandas as pd
P = sys.argv[1].rstrip("/") + "/data/"; BED = 225.87
st = pd.read_csv(P + "stage_walferdange.csv"); st["stage"] = st.wse - BED; st.to_csv(P + "stage_walferdange.csv", index=False)
so = pd.read_csv(P + "stage_obs_walferdange.csv").set_index("h").stage; po = so.loc[1:96]; ho = po.idxmax()
rows = []
for (c, m, i), d in st.groupby(["case", "model", "infiltration"]):
    s = d.set_index("h").stage.loc[1:96].astype(float)
    rows.append(dict(case=c, model=m, inf=i, peak=round(s.max(), 2), dpeak=round(s.max() - po.max(), 2), dt=int(s.idxmax() - ho),
                     h1=round(s.loc[1], 2), rise=round(s.max() - s.loc[1], 2)))
pd.DataFrame(rows).to_csv(P + "stage_scores.csv", index=False)
print("observed peak stage %.2f m at h%d, stage at h1 %.2f m" % (po.max(), ho, so.loc[1]))
