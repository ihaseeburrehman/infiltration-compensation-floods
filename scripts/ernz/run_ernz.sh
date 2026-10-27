#!/bin/bash
# usage: sbatch --job-name=e18_<case>_<inf> run_ernz.sh <case> <noinf|ga> <hours> [startfile|none] [nGPU=1] [soil folder=soil_v3b] [run dir suffix]
#   LISFLOOD-FP 8.4 FV1, CUDA, double precision, Ernz Noire 2018 case. Rain file forcing/rain_<case>.nc.
#   e.g. ... 48 <startfile> 1 soil_v3b_ksx10 _ksx10 (sbatch at ECMWF does not pass environment variables, so they are arguments).
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --qos=ng
#SBATCH --account=luuni
#SBATCH --output=/ec/res4/scratch/lux0804/ernz2018/logs/run_%x_%j.log
set -euo pipefail
CASE=$1; INF=$2; HOURS=$3; START=${4:-none}; NG=${5:-1}; SOIL=${6:-soil_v3b}; SUFFIX=${7:-}
E=/ec/res4/scratch/lux0804/ernz2018; BIN=/perm/lux0804/LISFLOOD_versions/LISFLOOD-8.4/build_double/lisflood
R=$E/runs/E18_${CASE}_${INF}${SUFFIX}; mkdir -p $R/results
cat > $R/run.par <<PAR
dirroot                  $R/results
resroot                  1hr
DEMfile                  $E/inputs/dem10_cond.asc
manningfile              $E/inputs/manning.n.ascii
dynamicrainfile          $E/forcing/rain_${CASE}.nc
fv1
cuda
GPUs $NG
sim_time                 $((HOURS*3600))
initial_tstep            10
massint                  600
saveint                  3600
bcifile                  $E/inputs/Border.bci
nodata_elevation         -9999
drain_nodata
qoutput
PAR
[ "$START" != none ] && echo "startfile                $START" >> $R/run.par
if [ "$INF" = ga ]; then cat >> $R/run.par <<PAR
infiltration             greenampt
greenampt_ks             $E/inputs/$SOIL/ga_ks_mmh.asc
greenampt_psi            $E/inputs/$SOIL/ga_psi_mm.asc
greenampt_porosity       $E/inputs/$SOIL/ga_porosity.asc
greenampt_theta_init     $E/inputs/$SOIL/ga_theta_init.asc
PAR
fi
md5sum $E/inputs/dem10_cond.asc $E/inputs/manning.n.ascii $E/inputs/Border.bci $E/forcing/rain_${CASE}.nc > $R/inputs.md5
module load prgenv/gnu cuda/12.6 netcdf4/4.9.3 hpcx-openmpi/2.14.0-cuda
nvidia-smi --query-gpu=index,name --format=csv,noheader
/usr/bin/time -v -o $R/time.txt mpirun --map-by ppr:4:node --bind-to core -np $NG $BIN $R/run.par
cp $R/results/*.mass $R/ 2>/dev/null || true
