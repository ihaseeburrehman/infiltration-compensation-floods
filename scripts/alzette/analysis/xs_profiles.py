"""Cross section profiles along each transect: water depth and unit discharge q.n per segment at several hours before and at
the peak, from the per segment flux files. Also the share of the wetted width that carries almost no flow (standing water).
usage: python xs_profiles.py <npz file> <out pdf> [gauges...]"""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
z = np.load(sys.argv[1]); OUT = sys.argv[2]; gauges = sys.argv[3:] or sorted({k[:-2] for k in z.files if k.endswith("_t")})
steps = z["steps"]; OFF = [-24, -12, -6, 0, 6]
fig, ax = plt.subplots(2, len(gauges), figsize=(3.2 * len(gauges), 5.2), squeeze=False)
for j, g in enumerate(gauges):
    t, L, qn, wd = z[g + "_t"], z[g + "_L"], z[g + "_qn"], z[g + "_wd"]; Q = (qn * L).sum(1); ip = int(Q.argmax())
    for o in OFF:
        i = ip + o
        if i < 0 or i >= len(steps): continue
        lab = f"peak {o:+d} h" if o else f"peak (h {steps[ip]}, {Q[ip]:.0f} m3/s)"
        ax[0, j].plot(t, wd[i], lw=1.0 if o else 1.8, label=lab); ax[1, j].plot(t, qn[i], lw=1.0 if o else 1.8)
    wet = wd[ip] > 0.05; f = np.abs(qn[ip] * L); standing = wet & (f < 0.001 * f.sum())
    ax[0, j].set_title(f"{g}: wet {L[wet].sum():.0f} m, of which {L[standing].sum():.0f} m carry <0.1% each", fontsize=8)
    ax[0, j].set_ylabel("depth (m)"); ax[1, j].set_ylabel("q.n (m2/s)"); ax[1, j].set_xlabel("distance along section (m)")
    ax[0, j].legend(fontsize=6, frameon=False)
    print(f"{g}: peak h {steps[ip]} Q {Q[ip]:.1f}; wet width {L[wet].sum():.0f} m; standing (<0.1 % of Q each) {L[standing].sum():.0f} m, together {100*f[standing].sum()/f.sum():.2f} % of Q; "
          f"share of Q in cells with depth<0.3 m: {100*f[wet & (wd[ip]<0.3)].sum()/f.sum():.1f} %; negative (upstream) flux share {100*abs((qn[ip]*L)[(qn[ip]*L)<0].sum())/f.sum():.1f} %")
fig.tight_layout(); fig.savefig(OUT); print("saved", OUT)
