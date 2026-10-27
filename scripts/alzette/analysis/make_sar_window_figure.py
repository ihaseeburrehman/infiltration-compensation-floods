"""Appendix map: Sentinel-1 flood (EMSN139, 15 July 05:50 UTC) against the simulated flood at hour 36 (radar, LISFLOOD-FP,
no infiltration) in 1.2 km windows around Livange, Steinsel and Walferdange, with the gauge transect and the station section.
Also prints containment and lateral reach statistics within 300 m of each gauge. usage: python make_sar_window_figure.py <paper_v4>"""
import sys, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch; from matplotlib.lines import Line2D
P = sys.argv[1].rstrip("/"); z = np.load(P + "/data/sar_windows.npz"); T = json.load(open(P + "/scripts/extraction/transects_fullvalley.json"))
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.grid": False})
ST = {"Livange": (-200, 200), "Steinsel": (-400, 300)}; CS = 10.0
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9))
for ax, g in zip(axs, ["Livange", "Steinsel", "Walferdange"]):
    sar, aoi, perm, none = [z[g + "_" + k] for k in ("sar", "aoi", "perm", "none")]; x0, yt = z[g + "_origin"]; n = sar.shape[0]
    ext = (x0, x0 + n * CS, yt - n * CS, yt); wet = (none > 0.1) & ~perm
    img = np.zeros(sar.shape + (3,)); img[:] = 1.0
    img[wet] = (0.78, 0.86, 0.94); img[sar & ~perm] = (0.03, 0.32, 0.61); img[perm] = (0.35, 0.35, 0.35)
    ax.imshow(img, extent=ext, origin="upper", interpolation="none", rasterized=True)
    v = dict(T[g]); v["ux"], v["uy"] = -T[g]["uy"], T[g]["ux"]   # (ux, uy) in the file is the river tangent; the section runs along (-uy, ux)
    lo, hi = min(s[2] for s in v["segs"]), max(s[3] for s in v["segs"])
    ax.plot([v["x"] + v["ux"] * lo, v["x"] + v["ux"] * hi], [v["y"] + v["uy"] * lo, v["y"] + v["uy"] * hi], color="#d62728", lw=0.8, ls=":")
    if g in ST:
        a, b = ST[g]; ax.plot([v["x"] + v["ux"] * a, v["x"] + v["ux"] * b], [v["y"] + v["uy"] * a, v["y"] + v["uy"] * b], color="#d62728", lw=1.6)
    ax.plot(v["x"], v["y"], "o", color="#d62728", ms=4, mec="white", mew=0.6)
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_xticks([]); ax.set_yticks([])
    ax.plot([ext[0] + 60, ext[0] + 260], [ext[2] + 60, ext[2] + 60], color="black", lw=1.2); ax.text(ext[0] + 160, ext[2] + 85, "200 m", ha="center", fontsize=7)
    # statistics within 300 m of the gauge
    yy, xx = np.mgrid[0:n, 0:n]; X = x0 + (xx + 0.5) * CS; Y = yt - (yy + 0.5) * CS
    from scipy import ndimage
    dist = ndimage.distance_transform_edt(~perm) * CS
    for R in (300, 600):
        W = (X - v["x"]) ** 2 + (Y - v["y"]) ** 2 <= R * R; o = sar & ~perm & W & aoi; m = wet & W & aoi
        print(f"{g:12s} R={R}: SAR {o.sum()/100:5.1f} ha, of which in simulated flood {100*(o&m).sum()/max(o.sum(),1):4.0f} %; simulated {m.sum()/100:5.1f} ha; "
              f"reach from channel SAR max {dist[o].max() if o.any() else 0:4.0f} p95 {np.percentile(dist[o],95) if o.any() else 0:4.0f} m, model max {dist[m].max():4.0f} p95 {np.percentile(dist[m],95):4.0f} m")
    ax.set_title(f"({'abc'[list(axs).index(ax)]}) {g}", fontsize=8, loc="left")
hd = [Patch(color=(0.78, 0.86, 0.94), label="Simulated flood, hour 36 (depth > 0.1 m)"), Patch(color=(0.03, 0.32, 0.61), label="Sentinel-1 flood, 15 July 05:50 UTC"),
      Patch(color=(0.35, 0.35, 0.35), label="Channel (wet at the start)"), Line2D([], [], color="#d62728", lw=0.8, ls=":", label="Transect"),
      Line2D([], [], color="#d62728", lw=1.6, label="Station section"), Line2D([], [], color="#d62728", marker="o", lw=0, ms=4, label="Gauge")]
fig.legend(handles=hd, loc="lower center", ncol=3, frameon=False, fontsize=7); fig.tight_layout(rect=(0, 0.12, 1, 0.97))
fig.savefig(P + "/figures/figB2_sar_windows.png", dpi=300); print("saved")   # PNG only: the 2 bit indexed rasters of the PDF render striped in poppler
