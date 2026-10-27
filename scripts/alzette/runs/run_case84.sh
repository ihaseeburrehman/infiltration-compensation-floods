#!/bin/bash
# usage: sbatch --job-name=l84_<case> run_case84.sh <case> <noinf|ga> [sp|dp] [nGPU, default 4; sbatch --ntasks=n --gres=gpu:n]
#   LISFLOOD-FP 8.4 (single precision, CUDA FV1, MPI on 4 GPUs), inputs identical to the HiPIMS cases
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --qos=ng
#SBATCH --account=luuni
#SBATCH --output=/ec/res4/scratch/lux0804/radar_scale_2026/logs/run84_%x_%j.log
set -euo pipefail
CASE=$1; INF=$2; PREC=${3:-sp}; NG=${4:-4}
W=/ec/res4/scratch/lux0804/radar_scale_2026; I=/ec/res4/scratch/lux0804/lf84_cases/inputs
case $CASE in
  R_1km_5min) RAIN=$W/forcing/rain_R_1km_5min.nc ;;
  *)          RAIN=$I/rain_${CASE}_b50.nc ;;
esac
SUF=""; BIN=/ec/res4/scratch/lux0804/lf84_float/build/lisflood
if [ "$PREC" = dp ]; then SUF=_dp; BIN=/ec/res4/scratch/lux0804/LISFLOOD_versions/LISFLOOD-8.4-multigpu/build/lisflood; fi
[ "$NG" != 4 ] && SUF=${SUF}_${NG}g
R=$W/runs/L84_${CASE}_${INF}${SUF}; mkdir -p $R/results
cat > $R/run.par <<PAR
dirroot                  $R/results
resroot                  1hr
DEMfile                  $W/inputs/Alzette_sub_basin_10m_bridge_burn.asc
manningfile              $W/inputs/manning.n.ascii
dynamicrainfile          $RAIN
startfile                $W/inputs/13_july_18hr.wd
fv1
cuda
GPUs $NG
sim_time                 345600
initial_tstep            10
massint                  600
saveint                  3600
bcifile                  $W/inputs/Border.bci
nodata_elevation         -9999
drain_nodata
qoutput
PAR
GI=$I; [ "$INF" = ga2 ] && GI=/ec/res4/scratch/lux0804/lf84_cases/inputs_v2; [ "$INF" = ga3 ] && GI=/ec/res4/scratch/lux0804/lf84_cases/inputs_v3a; [ "$INF" = ga3b ] && GI=/ec/res4/scratch/lux0804/lf84_cases/inputs_v3b; case $INF in ga3b_ks*) GI=/ec/res4/scratch/lux0804/lf84_cases/inputs_v3b_${INF#ga3b_} ;; esac
if [ "$INF" = ga ] || [ "$INF" = ga2 ] || [ "$INF" = ga3 ] || [ "${INF%%_ks*}" = ga3b ]; then cat >> $R/run.par <<PAR
infiltration             greenampt
greenampt_ks             $GI/ga_ks_mmh.asc
greenampt_psi            $GI/ga_psi_mm.asc
greenampt_porosity       $GI/ga_porosity.asc
greenampt_theta_init     $GI/ga_theta_init.asc
PAR
fi
md5sum $W/inputs/Alzette_sub_basin_10m_bridge_burn.asc $W/inputs/manning.n.ascii $W/inputs/13_july_18hr.wd $W/inputs/Border.bci $RAIN > $R/inputs.md5
module load prgenv/gnu cuda/12.6 netcdf4/4.9.3 hpcx-openmpi/2.14.0-cuda
nvidia-smi --query-gpu=index,name --format=csv,noheader
/usr/bin/time -v -o $R/time.txt mpirun --bind-to none -np $NG $BIN $R/run.par
cp $R/results/*.mass $R/ 2>/dev/null || true
