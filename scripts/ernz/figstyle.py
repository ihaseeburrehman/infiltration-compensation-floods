"""One palette and text standard for every figure of the paper (imported by all plotting scripts).
Colour encodes the forcing by family where lines or points stand for forcings:
  radar black; WRF green family (before DA teal green, after DA dark green); data driven models in three distinct hues (GraphCast blue, AIFS amber, FuXi magenta);
  rescaled forecasts (VM) dashed or open in the colour of their forecast.
Where rows already name the forcing (scorecards, sensitivity), colour encodes the infiltration setting in greys:
  no infiltration light grey, Green-Ampt dark grey; conductivity variants on a grey ramp, darker = larger losses.
Observations are always red dots. Minimum pairwise OKLab distance (x100): 19.7 normal vision, 9.9 deuteranopia, 11.9 protanopia; observations also differ by marker (dots)."""
OBS = "#e8473c"
RADAR = "#000000"
FAM = {"F_graphcast": "#2166ac", "F_aifs": "#f0a400", "F_fuxi": "#ae017e", "F_wrf12noda": "#5ab4ac", "F_wrf12": "#00441b"}
INF = {"none": "#a6a6a6", "GA": "#3d3d3d"}
RC = {"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
      "font.family": "sans-serif", "axes.linewidth": 0.6, "pdf.fonttype": 42, "savefig.dpi": 600}

def col(case):
    """Colour of a forcing key such as R_1km_5min, F_aifs, VM_wrf12, F_fuxi_ksdeep."""
    c = case.replace("VM_", "F_").replace("_ksdeep", "")
    return FAM.get(c, RADAR)

def inf_col(inf):
    return INF["none"] if inf in ("none", "noinf") else INF["GA"]

def panel(a, i, x=0.0, y=1.02):
    """Panel label (a), (b), ... just above the upper left corner of axes a (HESS: lower case letters in brackets)."""
    a.text(x, y, f"({'abcdefghijklmnopqrstuvwxyz'[i]})", transform=a.transAxes, ha="left", va="bottom", fontsize=8)
