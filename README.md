# Infiltration and compensating errors in forecast driven flood simulations

Code, run files and derived data for the paper

> Rehman, H. ur, Teferle, F. N. R., and Schumann, G.: Right for the wrong reasons: data driven and physics based rainfall forecasts match two Luxembourg floods only when soil infiltration is omitted, manuscript for Hydrology and Earth System Sciences.

Two shallow water models (LISFLOOD-FP and HiPIMS) are forced with radar rainfall and with GraphCast, AIFS, FuXi and WRF forecasts, without and with Green–Ampt infiltration, for the Alzette flood of July 2021 and the Ernz flash flood of 1 June 2018 in Luxembourg.

## Layout

```
scripts/alzette/forcing     radar, forecast and rescaled rainfall forcings
scripts/alzette/soil        Green–Ampt parameter grids (Ks, suction head, porosity, initial water content)
scripts/alzette/runs        run submission for LISFLOOD-FP and HiPIMS
scripts/alzette/extraction  discharge across the gauge sections and stage at Walferdange
scripts/alzette/analysis    scores, water balance, Sentinel-1 extent, tables and figures
scripts/ernz                the Ernz 2018 case: domain, forcings, soil, runs, extraction, scores, figures
runs/alzette_lisflood       run.par of every Alzette LISFLOOD-FP run in the paper
runs/ernz_lisflood          run.par of every Ernz LISFLOOD-FP run
data/alzette                derived results: scores, simulated hydrographs and stage, water balance,
                            Sentinel-1 scores, section fluxes, width convergence, LiDAR section profiles
data/ernz                   derived results of the Ernz case
```

Run names: `R_1km_5min` radar, `F_<model>` forecast as issued, `VM_<model>` forecast rescaled to the radar event total; `noinf` without infiltration, `ga3b` (Alzette) or `ga` (Ernz) with Green–Ampt; `ks05`, `ks2`, `ksx2`, `ksx10`, `ksdeep` conductivity variants; `sm` wet start; `dp` double precision.

## Order of the work

1. Forcings: `scripts/alzette/forcing/radclim_crop.py`, `make_forcing.py`, `lf84_inputs*.py`; Ernz `radklim_forcing.py`, `eurad_forcing.py`, `wrf12_forcing.py`, `ai_to_lisflood_rain.py`.
2. Soil: `soil_v3a.py` then `soil_v3b.py` (+ `soil_sm.py`, `soil_wetmask.py`); Ernz `soil_ernz.py`, `soil_ernz_variants.py`.
3. Runs: `run_case84.sh` and the `*_submit.sh` wrappers (LISFLOOD-FP), `hipims_setup.py` and `run_hipims_mg.sh` (HiPIMS); Ernz `run_ernz.sh`, `submit_event.sh`.
4. Extraction: `xs_fullvalley.py` with `seg_extract.py` and the transect files, `stage_walferdange.py`; Ernz `xs_ernz.py`, `xs_ernz_seg.py`.
5. Scores and figures: `score_all.py` (or `rescore_local.py` from `data/alzette/hydrographs.csv`), `water_balance.py`, `stage_scores.py`, `sar_extent.py`, `obs_uncertainty.py`, then `make_figures.py`, `make_scorecards.py`, `make_perf_diagram.py`, `make_xs_figure.py`, `make_sar_window_figure.py`; Ernz `score_ernz.py`, `fig_ernz_hydro.py`, `fig_ernz_ks.py`. `numbers_dump.py` prints every score derived number quoted in the paper.

The scripts were run on the ECMWF high performance computer and keep its absolute paths; adapt the paths at the top of each script.

## Not included

| Item | Source |
|---|---|
| LISFLOOD-FP | 8.3 on Zenodo, https://doi.org/10.5281/zenodo.13121119. The extension used here (Green–Ampt infiltration, several GPUs) will be released with LISFLOOD-FP 8.4. |
| HiPIMS | from its developers (Xia et al., 2019) |
| Observed discharge and stage | Administration de la gestion de l'eau (AGE), Luxembourg; Larochette via CAMELS-LUX |
| Radar rainfall | RADCLIM, Royal Meteorological Institute of Belgium; RADKLIM, Deutscher Wetterdienst; EURADCLIM |
| Soil hydraulic data | EU-SoilHydroGrids, European Soil Data Centre (licence does not allow redistribution); SoilGrids 2.0 |
| Terrain | LiDAR terrain model of Luxembourg, Administration du cadastre et de la topographie |
| Sentinel-1 flood extent | Copernicus Emergency Management Service, EMSN139 (© European Union, 2021) |
| Model input grids and result grids | rebuilt from the sources above with the scripts; the run.par files give every setting |

## Requirements

Python 3.11 with numpy, pandas, scipy, matplotlib, xarray, netCDF4, h5py, rasterio, geopandas, shapely, pyproj, affine; `hipims_io` for HiPIMS. See `requirements.txt`.

## Licence

Code under the MIT licence (`LICENSE`). Derived data in `data/` and `runs/` under CC BY 4.0. Funded by the Luxembourg National Research Fund (FNR), Industrial Fellowship, project no. 17130773.
