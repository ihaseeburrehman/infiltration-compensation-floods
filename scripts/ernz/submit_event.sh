#!/bin/bash
# Submit the Ernz 2018 event runs (48 h from 31 May 2018 18 UTC) for all forcings, with and without infiltration.
# usage: bash submit_event.sh <startfile path or none>
START=${1:-none}; cd /ec/res4/scratch/lux0804/ernz2018/code
for CASE in R_1km_5min R_1km_1h_nocorr R_eurad_1h R_eurad_1h_corr F_graphcast F_aifs F_fuxi F_wrf12; do
  for INF in noinf ga; do sbatch --job-name=e18_${CASE}_${INF} run_ernz.sh $CASE $INF 48 $START; done
done
squeue -u lux0804 -o "%.10i %.28j %.8T %.10M" | grep e18
