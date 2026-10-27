"""SAR wet cells along each gauge transect: EMSN139 Sentinel-1 masks (15 July 05:51 UTC, 16 July) sampled at the transect
segment cells; prints per gauge the SAR-wet length, contiguous span around the gauge, extremes, and AOI coverage."""
import json, numpy as np
M = np.load("/ec/res4/scratch/lux0804/radar_scale_2026/code/sar_masks.npz"); T = json.load(open("transects_fullvalley.json"))
aoi = M["aoi01"] | M["aoi02"]
for key in ("sar15", "sar16"):
    S = M[key]
    for g, v in T.items():
        segs = v["segs"]; r = np.array([s[0] for s in segs]); c = np.array([s[1] for s in segs])
        s0 = np.array([s[2] for s in segs]); s1 = np.array([s[3] for s in segs]); L = s1 - s0; t = 0.5 * (s0 + s1)
        wet = S[r, c]; cov = aoi[r, c]
        # contiguous wet span containing or nearest to the gauge (t=0), gaps up to 30 m bridged
        idx = np.where(wet)[0]
        span = (np.nan, np.nan)
        if len(idx):
            groups = np.split(idx, np.where(np.diff(t[idx]) > 30)[0] + 1)
            best = min(groups, key=lambda gi: min(abs(t[gi].min()), abs(t[gi].max())) if not (t[gi].min() <= 0 <= t[gi].max()) else 0)
            span = (s0[best].min(), s1[best].max())
        print(f"{key} {g:12s} aoi cover {100*L[cov].sum()/L.sum():5.1f} %  SAR wet {L[wet].sum():6.1f} m of {L.sum():.0f} m  "
              f"extremes {s0[wet].min() if wet.any() else np.nan:7.1f}..{s1[wet].max() if wet.any() else np.nan:7.1f}  "
              f"span around gauge {span[0]:7.1f}..{span[1]:7.1f} ({(span[1]-span[0]) if wet.any() else 0:.0f} m)")
        print("   wet string:", "".join("#" if w else ("." if cv else " ") for w, cv in zip(wet, cov)))
