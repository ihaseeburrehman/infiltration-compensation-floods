"""Full-valley transects: same centres/normals as transects400.json, but the along-line extent covers the whole
flooded floodplain at the peak (Livange, Walferdange were 400 m and missed overbank flow), stopping short of
neighbouring tributaries. Then extract q.n for LISFLOOD-FP 8.4 and HiPIMS (1-GPU or 4-GPU) runs in parallel.
usage: xs_fullvalley.py build | lf <run_dir> | hp <run_dir>"""
import sys, os, glob, re, json, numpy as np, multiprocessing as mp
sys.path.insert(0, "/ec/res4/scratch/lux0804/radar_scale_2026/code")
import seg_extract as se
CODE = "/ec/res4/scratch/lux0804/radar_scale_2026/code/"
X0, YTOP, CS, NR, NC = 59799.3, 49254.0 + 4362 * 10.0, 10.0, 4362, 2440
EXT = {"Livange": (-700, 400), "Walferdange": (-400, 800), "Steinsel": (-400, 300),
       "Pfaffenthal": (-200, 200), "Hesperange": (-200, 200)}

def build():
    T = se.load(CODE + "transects400.json"); N = {}
    for n, v in T.items():
        a, b = EXT[n]; ds = 0.05; s = np.arange(a, b, ds) + ds / 2
        px, py = v["x"] - v["uy"] * s, v["y"] + v["ux"] * s
        r = np.floor((YTOP - py) / CS).astype(int); c = np.floor((px - X0) / CS).astype(int)
        cut = np.flatnonzero((np.diff(r) != 0) | (np.diff(c) != 0)) + 1
        segs = [[int(r[i]), int(c[i]), float(s[i] - ds / 2), float(s[j - 1] + ds / 2)]
                for i, j in zip(np.r_[0, cut], np.r_[cut, len(s)])]
        N[n] = dict(v, segs=segs)
    json.dump(N, open(CODE + "transects_fullvalley.json", "w"))
    print({n: (len(v["segs"]), round(sum(q[3] - q[2] for q in v["segs"]))) for n, v in N.items()})

T = se.load(CODE + "transects_fullvalley.json") if os.path.exists(CODE + "transects_fullvalley.json") else None

def readers(run, kind):
    if kind == "lf":
        R = run + "/results/"
        ts = sorted(int(re.search(r"-(\d+)\.Qx", p).group(1)) for p in glob.glob(R + "1hr-*.Qx"))
        f = lambda t, v: R + f"1hr-{t:04d}." + {"x": "Qx", "y": "Qy", "h": "wd"}[v]
        return ts, (lambda t, v, need: se.read_rows(f(t, v), need))
    MG = sorted(glob.glob(run + "/[0-9]/output/")) or [run + "/output/"]
    name = {"x": "hUx", "y": "hUy", "h": "h"}
    DOM = []
    for d in MG:
        with open(glob.glob(d + "h_*.asc")[0]) as fh: h = {k: float(x) for k, x in (fh.readline().split() for _ in range(6))}
        nr = int(h["nrows"]); off = int(round((YTOP - (h["yllcorner"] + nr * CS)) / CS))
        DOM.append((d, off, 0 if off == 0 else 1, nr if off + nr == NR else nr - 1))   # skip halo rows
    def rd(t, v, need):
        out = {}
        for d, off, lo, hi in DOM:
            loc = [g - off for g in need if lo <= g - off < hi]
            if loc: out.update({l + off: x for l, x in se.read_rows(f"{d}{name[v]}_{t}.asc", loc).items()})
        return out
    ts = sorted(int(re.search(r"hUx_(\d+)\.asc", p).group(1)) for p in glob.glob(DOM[0][0] + "hUx_*.asc"))
    return [t for t in ts if t % 3600 == 0 and t > 0], rd

def one(args):
    t, = args
    QX, QY, WD = (RD(t, v, NEED) for v in "xyh")
    return se.seg_flux(T, QX, QY, WD)

if __name__ == "__main__":
    if sys.argv[1] == "build": build(); sys.exit()
    kind, run = sys.argv[1], sys.argv[2].rstrip("/")
    NEED = se.rows_needed(T); TS, RD = readers(run, kind)
    with mp.Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", 1))) as pool: RES = pool.map(one, [(t,) for t in TS])
    save = {"steps": np.array([t if kind == "lf" else t // 3600 for t in TS])}
    for n in T:
        save[f"{n}_t"] = RES[-1][n]["t"]; save[f"{n}_L"] = RES[-1][n]["L"]
        save[f"{n}_qn"] = np.array([r[n]["qn"] for r in RES]); save[f"{n}_wd"] = np.array([r[n]["wd"] for r in RES])
    np.savez_compressed(run + "/segflux_fv.npz", **save)
    print(run, len(TS), "steps")
