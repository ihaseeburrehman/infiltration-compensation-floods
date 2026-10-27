"""Per-segment q.n and depth on the Ernz transects for every hourly output of a LISFLOOD-FP run (for the width convergence check).
usage: python xs_ernz_seg.py <run dir>  -> <run dir>/segflux_seg.npz (same layout as xs_fullvalley.py of paper_v4)"""
import sys, os, glob, re, numpy as np, multiprocessing as mp
sys.path.insert(0, "/ec/res4/scratch/lux0804/radar_scale_2026/code"); import seg_extract as se
RUN = sys.argv[1].rstrip("/"); T = se.load("/ec/res4/scratch/lux0804/ernz2018/code/transects_ernz.json"); NEED = se.rows_needed(T); R = RUN + "/results/"
def one(t):
    f = lambda v: R + f"1hr-{t:04d}." + v
    QX, QY, WD = (se.read_rows(f(v), NEED) for v in ("Qx", "Qy", "wd")); return se.seg_flux(T, QX, QY, WD)
if __name__ == "__main__":
    TS = sorted(int(re.search(r"-(\d+)\.Qx", p).group(1)) for p in glob.glob(R + "1hr-*.Qx"))
    with mp.Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", 4))) as pool: RES = pool.map(one, TS)
    save = {"steps": np.array(TS)}
    for n in T:
        save[f"{n}_t"] = RES[-1][n]["t"]; save[f"{n}_L"] = RES[-1][n]["L"]
        save[f"{n}_qn"] = np.array([r[n]["qn"] for r in RES]); save[f"{n}_wd"] = np.array([r[n]["wd"] for r in RES])
    np.savez_compressed(RUN + "/segflux_seg.npz", **save); print(RUN, len(TS), "steps")
