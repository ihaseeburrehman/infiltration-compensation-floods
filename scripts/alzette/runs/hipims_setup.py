"""Build a HiPIMS case for one rainfall case, using the LISFLOOD-FP inputs unchanged.
usage: hipims_setup.py CASE [NGPU]   (R_1km_5min supported: radar pixels as rain sources)"""
import sys, os, numpy as np, pandas as pd
if not hasattr(np, "trapz"): np.trapz = np.trapezoid   # hipims_io 0.6.3 predates NumPy 2
os.environ.setdefault("MPLCONFIGDIR", "/ec/res4/scratch/lux0804/hipims/.mpl")
from hipims_io import InputHipims, Raster
from pyproj import Transformer
CASE = sys.argv[1]; NG = int(sys.argv[2]) if len(sys.argv) > 2 else 1
import re
m_ = re.match(r"(.*?)_inf(?:_(\w+))?$", CASE); INF = bool(m_); BASE = m_.group(1) if INF else CASE; KS = m_.group(2) if INF else None
W = "/ec/res4/scratch/lux0804/radar_scale_2026"; H = "/ec/res4/scratch/lux0804/hipims"
dem = Raster(W + "/inputs/Alzette_sub_basin_10m_bridge_burn.asc")
NR, NC = dem.array.shape; print("DEM", dem.array.shape, dem.header, flush=True)
rd = lambda f: pd.read_csv(f, skiprows=6, sep=r"\s+", header=None, dtype=np.float32, engine="c").values
h0 = rd(W + "/inputs/13_july_18hr.wd"); man = rd(W + "/inputs/manning.n.ascii")
act = ~np.isnan(dem.array) if np.isnan(dem.array).any() else dem.array > -9000
case = InputHipims(dem, num_of_sections=NG, case_folder=f"{H}/runs/{CASE}" + (f"_mg{NG}" if NG > 1 else ""))
case.set_initial_condition("h0", np.where(act, h0, 0.0))
case.set_grid_parameter(manning=np.where(act, man, 0.035))
case.set_boundary_condition(outline_boundary="open")
case.set_runtime([0, 345600, 3600, 345600])
case.set_device_no(list(range(NG)) if NG > 1 else 0)
if BASE == "R_1km_5min":
    d = np.load(W + "/radclim_crop.npz")
    r5 = np.nan_to_num(d["r5"], nan=0.0); h1 = np.nan_to_num(d["r1h"], nan=0.0); rx, ry = d["x3812"], d["y3812"]
    ny, nx = h1.shape[1:]
    r5h = r5.reshape(96, 12, ny, nx); s = r5h.sum(axis=1)
    w = np.where(s[:, None] > 0, r5h / np.where(s > 0, s, 1)[:, None], 1.0 / 12)
    r5c = (w * h1[:, None]).reshape(1152, ny * nx)                 # mm per 5 min, identical to LISFLOOD R1
    x = 59799.282445142424 + (np.arange(NC) + .5) * 10; y = 49254.0 + NR * 10 - (np.arange(NR) + .5) * 10
    X, Y = np.meshgrid(x, y); px, py = Transformer.from_crs("EPSG:2169", "EPSG:3812", always_xy=True).transform(X.ravel(), Y.ravel())
    flat = (np.floor(((ry[0] + 500) - py) / 1000).astype(int) * nx + np.floor((px - (rx[0] - 500)) / 1000).astype(int)).reshape(NR, NC)
    used, mask = np.unique(flat[act], return_inverse=True)
    rain_mask = np.zeros((NR, NC), np.int32); rain_mask[act] = mask
    rate = r5c[:, used] / 1000.0 / 300.0                             # m/s, constant within each 5 min interval
    t0 = np.arange(1152) * 300.0
    t = np.ravel(np.column_stack([t0, t0 + 299.0])); v = np.repeat(rate, 2, axis=0)
    t = np.append(t, [345600.0, 345601.0]); v = np.vstack([v, np.zeros((2, len(used)))])
    case.set_rainfall(rain_mask=rain_mask, rain_source=np.column_stack([t, v]))
    tot = (rate.sum(0) * 300 * 1000)[mask].mean(); print("sources", len(used), "basin mean rain %.2f mm" % tot, flush=True)
else:
    import netCDF4
    src = {"F_graphcast": "inputs_18utc/rain_graphcast_v2_96h.nc", "F_aifs": "inputs_18utc/rain_aifs_96h.nc",
           "VM_graphcast": "forcing/rain_VM_graphcast.nc", "R_1km_6h": "forcing/rain_R_1km_6h.nc", "VM_aifs": "forcing/rain_VM_aifs.nc",
           "F_wrf1p3": "inputs_18utc/rain_wrf1p3_afterda_96h.nc", "F_wrf12": "inputs_18utc/rain_wrf12_afterda_96h.nc",
           "F_wrf12noda": "inputs_18utc/rain_wrf12_beforeda_96h.nc", "VM_wrf1p3": "forcing/rain_VM_wrf1p3.nc", "F_fuxi": "inputs_18utc/rain_fuxi_noflip_96h.nc", "VM_fuxi": "forcing/rain_VM_fuxi2.nc",
           "VM_wrf12": "forcing/rain_VM_wrf12.nc", "VM_wrf12noda": "forcing/rain_VM_wrf12noda.nc"}[BASE]
    nc = netCDF4.Dataset(W + "/" + src); tt = np.asarray(nc["time"][:]) * 3600.0; v = nc["rainfall_depth"]
    B = 5                                                             # 50 m blocks of 5 x 5 cells
    bi, bj = np.arange(NR) // B, np.arange(NC) // B
    blk = (bi[:, None] * (NC // B + 1) + bj[None, :])
    used, mask = np.unique(blk[act], return_inverse=True)
    rain_mask = np.zeros((NR, NC), np.int32); rain_mask[act] = mask
    cnt = np.bincount(mask)
    rows = []
    for i in range(len(tt) - 1):
        a_ = np.asarray(v[i], dtype=np.float64)
        if a_.shape != (NR, NC): a_ = a_.T
        depth = np.bincount(mask, a_[act]) / cnt                      # block-mean depth (mm), volume preserving
        rows.append(depth / 1000.0 / (tt[i + 1] - tt[i]))             # m/s
    rate = np.array(rows); dt = np.diff(tt)
    t = np.ravel(np.column_stack([tt[:-1], tt[1:] - 1.0])); vv = np.repeat(rate, 2, axis=0)
    t = np.append(t, [tt[-1], tt[-1] + 1.0]); vv = np.vstack([vv, np.zeros((2, len(used)))])
    case.set_rainfall(rain_mask=rain_mask, rain_source=np.column_stack([t, vv]))
    tot = ((rate * dt[:, None]).sum(0) * 1000)[mask].mean(); print("sources", len(used), "basin mean rain %.2f mm" % tot, flush=True)
if INF:   # Green-Ampt v3b: EU-SoilHydroGrids Ks and porosity (0-30 cm), v2 wet mask, satellite soil moisture 13 July 2021
    ga = np.load(H + "/soil/green_ampt_10m_v3b" + ("_" + KS if KS else "") + ".npz"); print("soil file", KS or "v3b", flush=True)
    case.set_grid_parameter(hydraulic_conductivity=np.where(act, ga["Ks"], 0.0),
                            capillary_head=np.where(act, ga["psi"], 0.0),
                            water_content_diff=np.where(act, ga["dtheta"], 0.0))
    print("infiltration on: basin mean Ks %.2e m/s" % ga["Ks"][act].mean(), flush=True)
case.write_input_files()
print("written", case.get_case_folder(), flush=True)
