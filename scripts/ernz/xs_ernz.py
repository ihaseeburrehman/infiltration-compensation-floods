"""Discharge through the Ernz transects (transects_ernz.json) for a LISFLOOD-FP run: Q = sum(q.n * L) over the full transect,
as xs_fullvalley.py of paper_v4 (full wetted section). Grid constants from the run's DEM header.
usage: python xs_ernz.py <run dir> [transects json]   -> <run dir>/segflux_fv.npz and hydrograph csv"""
import sys, os, glob, re, json, numpy as np, multiprocessing as mp
sys.path.insert(0, "/ec/res4/scratch/lux0804/radar_scale_2026/code"); import seg_extract as se
RUN = sys.argv[1].rstrip("/"); TJ = sys.argv[2] if len(sys.argv) > 2 else "/ec/res4/scratch/lux0804/ernz2018/code/transects_ernz.json"
T = se.load(TJ); NEED = se.rows_needed(T); R = RUN + "/results/"
def one(t):
    f = lambda v: R + f"1hr-{t:04d}." + v
    QX, QY, WD = (se.read_rows(f(v), NEED) for v in ("Qx", "Qy", "wd"))
    r = se.seg_flux(T, QX, QY, WD); return {n: (float((v["qn"] * v["L"]).sum()), float(np.nanmax(v["wd"]))) for n, v in r.items()}
if __name__ == "__main__":
    TS = sorted(int(re.search(r"-(\d+)\.Qx", p).group(1)) for p in glob.glob(R + "1hr-*.Qx"))
    with mp.Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", 4))) as pool: RES = pool.map(one, TS)
    names = list(T); Q = np.array([[r[n][0] for n in names] for r in RES]); H = np.array([[r[n][1] for n in names] for r in RES])
    np.savez(RUN + "/segflux_fv.npz", steps=np.array(TS), names=names, Q=Q, hmax=H)
    with open(RUN + "/hydrograph.csv", "w") as f:
        f.write("hour," + ",".join(names) + "\n"); [f.write(f"{t}," + ",".join(f"{q:.3f}" for q in row) + "\n") for t, row in zip(TS, Q)]
    print(RUN.split("/")[-1], {n: "peak %.1f m3/s at %d h" % (Q[:, i].max(), TS[Q[:, i].argmax()]) for i, n in enumerate(names)})
