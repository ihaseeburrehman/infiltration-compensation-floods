"""Fig. 1: location map and the two study areas side by side, in one style.
(a) Luxembourg and its neighbours (Natural Earth 10 m) with the two model domains shaded red.
(b) Upper Alzette, radar total 13 July 18:00 to 17 July 18:00 UTC 2021 (RADCLIM), gauges, Sentinel-1 flood of 15 July 05:50 UTC.
(c) Ernz valleys, corrected radar total 31 May 18:00 to 2 June 18:00 UTC 2018 (RADKLIM), Larochette gauge and catchment, rain gauges.
Both maps: hillshaded 10 m terrain, radar total in the same colour scale, rivers in the same blue, gauges as the same symbol.
usage: python make_fig01_study_areas.py <paper_v4 dir> <ernz2018 dir> <alzette dem asc> <alzette river shp> <emsn139 gdb> <camels catchments shp>"""
import sys, json, numpy as np, geopandas as gpd, rasterio, netCDF4 as nc, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MultipleLocator
from pyproj import Transformer
from rasterio import features
from shapely.geometry import shape
from shapely.ops import unary_union

P, E, DEM, RIV, GDB, CAM = [a.rstrip("/") for a in sys.argv[1:7]]
W = P + "/data_raw/fig01_inputs"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42, "font.family": "sans-serif"})
INK, RIVC, FLOOD, RED = "#1a1a1a", "#1f5fa8", "#c51b7d", "#d62728"
CMAP, VMIN, VMAX, ALPHA = "Blues", 30, 130, 0.6
lsrc = LightSource(azdeg=315, altdeg=45)
kmf = FuncFormatter(lambda v, p: f"{v / 1e3:.0f}")

def hill(dem, valid):
    hs = lsrc.hillshade(np.nan_to_num(dem, nan=np.nanmin(dem)), vert_exag=2, dx=10, dy=10)
    return np.where(valid, hs, np.nan)

# ---------- Alzette inputs (LUREF grid of the model, 10 m)
with open(DEM) as f: hdr = [next(f).split() for _ in range(6)]
NC_, NR_ = int(hdr[0][1]), int(hdr[1][1]); X0, Y0, dx = float(hdr[2][1]), float(hdr[3][1]), float(hdr[4][1])
demA = np.loadtxt(DEM, skiprows=6); demA[demA <= -9998] = np.nan
B = np.load(W + "/basin_mask.npy").astype(bool); assert B.shape == demA.shape
extA = [X0, X0 + NC_ * dx, Y0, Y0 + NR_ * dx]
z = np.load(W + "/radclim_crop.npz"); totR = z["r1h"].sum(0)
s = 5  # rainfall resampled to 50 m by nearest radar pixel
xx = X0 + dx * (np.arange(0, NC_, s) + .5); yy = Y0 + dx * (NR_ - np.arange(0, NR_, s) - .5)
GX, GY = np.meshgrid(xx, yy); LX, LY = Transformer.from_crs(2169, 3812, always_xy=True).transform(GX, GY)
ix = np.clip(np.rint((LX - z["x3812"][0]) / 1000).astype(int), 0, len(z["x3812"]) - 1)
iy = np.clip(np.rint((LY - z["y3812"][0]) / -1000).astype(int), 0, len(z["y3812"]) - 1)
rainA = np.where(B[::s, ::s], totR[iy, ix], np.nan)
M = np.load(W + "/sar_masks.npz"); T = json.load(open(W + "/transects400.json"))
riv = gpd.read_file(RIV); riv = riv[riv.river_name == "Alzette"].to_crs(2169)
aoi = gpd.read_file(GDB, layer="P00_aoi").to_crs(2169)
basinA = unary_union([shape(g) for g, v in features.shapes(B.astype(np.uint8), mask=B,
                      transform=rasterio.transform.from_origin(X0, Y0 + NR_ * dx, dx, dx)) if v == 1])

# ---------- Ernz inputs
with rasterio.open(E + "/domain/dem10_mean.tif") as r: demE = r.read(1).astype(float); tf = r.transform; NRe, NCe = demE.shape
act = np.load(E + "/domain/domain_mask.npy"); demE[~act] = np.nan
extE = [tf.c, tf.c + NCe * tf.a, tf.f + NRe * tf.e, tf.f]
d = nc.Dataset(E + "/forcing/rain_R_1km_5min.nc"); v = d["rainfall_depth"]; rainE = np.zeros((NRe, NCe))
for i in range(v.shape[0]): rainE += np.array(v[i])
rainE[~act] = np.nan
noire = gpd.read_file(E + "/domain/ernz_noire.geojson").to_crs(2169); blanche = gpd.read_file(E + "/domain/ernz_blanche.geojson").to_crs(2169)
cam = gpd.read_file(CAM); cam = cam[cam.gauge_id == 55].to_crs(2169)
rivE = gpd.read_file(E + "/domain/lines/col_2091.geojson").to_crs(2169); rivE = rivE[rivE.Name_2021.isin(["Ernz noire", "Ernz blanche", "Halerbaach"])]
gap = np.load(E + "/domain/gap_mask.npy")
domE = unary_union(list(noire.geometry) + list(blanche.geometry))

# ---------- countries
ne = gpd.read_file(W + "/ne_10m_admin_0_countries.shp")
ne = ne[ne.ADMIN.isin(["Luxembourg", "Belgium", "France", "Germany", "Netherlands"])].to_crs(2169)

# ---------- layout: location map | Alzette | Ernz, same map scale in (b) and (c)
fig = plt.figure(figsize=(7.2, 5.0))
H = 0.80  # map height (figure fraction) of the Alzette panel
hA = (extA[3] - extA[2]); scale = H * 5.0 / hA  # inches per metre
wA = (extA[1] - extA[0]) * scale / 7.2
wE = (extE[1] - extE[0]) * scale / 7.2; hE = (extE[3] - extE[2]) * scale / 5.0
x0a = 0.36; axA = fig.add_axes([x0a, 0.10, wA, H]); axE = fig.add_axes([x0a + wA + 0.075, 0.10 + H - hE, wE, hE])
axL = fig.add_axes([0.02, 0.40, 0.255, 0.50])

# (a) location
L = ne[ne.ADMIN == "Luxembourg"]
ne.plot(ax=axL, color="white", edgecolor="#8c8c8c", lw=0.6); L.plot(ax=axL, color="#ececec", edgecolor=INK, lw=0.8)
gpd.GeoSeries([basinA, domE], crs=2169).plot(ax=axL, color=RED, alpha=0.75, edgecolor=RED, lw=0.4)
axL.text(basinA.centroid.x - 9000, basinA.centroid.y - 2000, "b", fontsize=8, fontweight="bold", color=RED, ha="right")
axL.text(domE.centroid.x + 8000, domE.centroid.y, "c", fontsize=8, fontweight="bold", color=RED, ha="left")
for n, (lx, ly) in {"BELGIUM": (32000, 112000), "GERMANY": (124000, 82000), "FRANCE": (60000, 30000), "LUXEMBOURG": (72000, 142500)}.items():
    axL.text(lx, ly, n, fontsize=6 if n == "LUXEMBOURG" else 6.5, ha="center", color="#555555" if n != "LUXEMBOURG" else INK, rotation=90 if n == "GERMANY" else 0, va="center", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none") if n == "LUXEMBOURG" else None)
axL.set_xlim(20000, 135000); axL.set_ylim(22000, 148000); axL.set_aspect("equal")
axL.set_xticks([]); axL.set_yticks([])
for sp in axL.spines.values(): sp.set_linewidth(0.6); sp.set_color("#8c8c8c")
axL.set_title("(a) Location", loc="left")

def base(ax, dem, valid, rain, ext, title):
    ax.imshow(hill(dem, valid), extent=ext, cmap="gray", vmin=0, vmax=1.3, interpolation="nearest", rasterized=True)
    im = ax.imshow(rain, extent=ext, cmap=CMAP, alpha=ALPHA, vmin=VMIN, vmax=VMAX, interpolation="nearest", rasterized=True)
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
    ax.xaxis.set_major_formatter(kmf); ax.yaxis.set_major_formatter(kmf)
    ax.xaxis.set_major_locator(MultipleLocator(10000)); ax.yaxis.set_major_locator(MultipleLocator(10000))
    ax.set_title(title, loc="left"); ax.set_xlabel("Easting (km)")
    for sp in ax.spines.values(): sp.set_linewidth(0.6)
    return im

def gauge(ax, x, y, name, dxp, dyp, ha="left"):
    ax.plot(x, y, "o", ms=4.5, mfc="white", mec=INK, mew=0.9, zorder=6)
    ax.annotate(name, (x, y), xytext=(dxp, dyp), textcoords="offset points", fontsize=7, ha=ha, va="center", zorder=7,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.7))

# (b) Alzette
im = base(axA, demA, B, rainA, extA, "(b) Alzette, 13 to 17 July 2021")
gpd.GeoSeries([basinA], crs=2169).boundary.plot(ax=axA, color=INK, lw=0.9, zorder=4)
sar = np.where(B, M["sar15"], 0).astype(float)
axA.contourf(X0 + dx * (np.arange(0, NC_, 2) + .5), Y0 + dx * (NR_ - np.arange(0, NR_, 2) - .5), sar[::2, ::2], [.5, 1.5], colors=[FLOOD], zorder=5)
gpd.GeoSeries([basinA], crs=2169).intersection(aoi.unary_union.boundary).plot(ax=axA, color="#555555", lw=0.6, ls="--", zorder=4)
riv.plot(ax=axA, color=RIVC, lw=1.0, zorder=4)
OFF = dict(Walferdange=(-6, 0, "right"), Steinsel=(6, 3, "left"), Pfaffenthal=(6, 0, "left"), Livange=(-6, 0, "right"), Hesperange=(6, 0, "left"))
for n, vv in T.items(): gauge(axA, vv["x"], vv["y"], n, *OFF[n])
axA.set_ylabel("Northing (km, LUREF)")

# (c) Ernz
base(axE, demE, act, rainE, extE, "(c) Ernz, 31 May to 2 June 2018")
for g in (noire, blanche): g.boundary.plot(ax=axE, color=INK, lw=0.9, zorder=4)
cam.boundary.plot(ax=axE, color=INK, lw=0.8, ls="--", zorder=4)
rivE.plot(ax=axE, color=RIVC, lw=1.0, zorder=4)
axE.contour(np.linspace(extE[0] + 5, extE[1] - 5, NCe), np.linspace(extE[3] - 5, extE[2] + 5, NRe), gap.astype(float), levels=[0.5], colors="#555555", linewidths=0.8, linestyles=":", zorder=4)
tr = Transformer.from_crs(4326, 2169, always_xy=True)
for n, (lo, la) in {"Waldbillig": (6.2773, 49.79806), "Echternach": (6.44337, 49.8031), "Godbrange": (6.23569, 49.73637), "Bettendorf": (6.2095, 49.8741)}.items():
    x, y = tr.transform(lo, la)
    if extE[0] < x < extE[1] and extE[2] < y < extE[3]:
        axE.plot(x, y, "^", color=INK, ms=4.5, mec="white", mew=0.5, zorder=6)
        axE.annotate(n, (x, y), xytext=(4, 3), textcoords="offset points", fontsize=7, zorder=7, bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.7))
gauge(axE, 83689, 95240, "Larochette", 5, -6)
x, y = tr.transform(6.3055, 49.7896); axE.plot(x, y, "s", ms=4, mfc=INK, mec="white", mew=0.5, zorder=6)
axE.annotate("Mullerthal", (x, y), xytext=(5, -2), textcoords="offset points", fontsize=7, zorder=7, bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.7))

# shared colour bar and legend under the location map
cax = fig.add_axes([0.05, 0.27, 0.22, 0.022])
cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="both"); cb.set_label("Radar event rainfall (mm)"); cb.outline.set_linewidth(0.5)
handles = [Line2D([], [], color=INK, lw=0.9, label="Catchment"), Line2D([], [], color=INK, lw=0.8, ls="--", label="Larochette catchment"),
           Line2D([], [], color=RIVC, lw=1.0, label="River"), Line2D([], [], ls="", marker="o", ms=4.5, mfc="white", mec=INK, label="Discharge or stage gauge"),
           Line2D([], [], ls="", marker="^", ms=4.5, mfc=INK, mec="white", label="Rain gauge"),
           Patch(color=FLOOD, label="Sentinel-1 flood, 15 July 2021"), Line2D([], [], color="#555555", lw=0.6, ls="--", label="Sentinel-1 mapping area"),
           Line2D([], [], color="#555555", lw=0.8, ls=":", label="Limit of RADKLIM coverage")]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(x0a + wA + 0.06, 0.10 + H - hE - 0.06), frameon=False, handlelength=1.6, labelspacing=0.45)
out = P + "/figures/fig01_study_areas"
fig.savefig(out + ".pdf", dpi=300); fig.savefig(out + ".png", dpi=400); print("saved", out)
