#!/bin/bash
# FuXi (and its rescaled version) in both models, and the ERA5-Land wet start (ga3b_sm / _inf_sm) for radar, GraphCast,
# AIFS and FuXi. LISFLOOD-FP uses the final 8.4 build (dpn); each run is chained to its extraction.
set -e
H=/ec/res4/scratch/lux0804/hipims; W=/ec/res4/scratch/lux0804/radar_scale_2026
fv() { sbatch --dependency=afterok:$1 --job-name=fv --qos=nf --account=luuni --cpus-per-task=32 --mem=32G --time=01:00:00 \
  --output=$H/logs/fv_$2.out --wrap="module load python3; cd $H && python3 xs_fullvalley.py $3 $4" >/dev/null; }
lf() {   # case inf ngpu dir
  j=$(sbatch --parsable --job-name=l84_$1_$2_dpn --ntasks=$3 --gres=gpu:$3 $W/code/run_case84.sh $1 $2 dpn $3)
  fv $j $4 lf $W/runs/$4; echo "lf $1 $2: $j"; }
hp() {   # hipims case
  s=$(sbatch --parsable --job-name=setup_$1 --qos=nf --account=luuni --cpus-per-task=4 --mem=64G --time=02:00:00 \
    --output=$H/logs/setup_$1_%j.log --wrap="cd $H && venv/bin/python hipims_setup.py $1 4")
  r=$(sbatch --parsable --dependency=afterok:$s --job-name=hpmg_$1 $H/run_hipims_mg.sh $1_mg4)
  fv $r $1_mg4 hp $H/runs/$1_mg4; echo "hp $1: $r"; }
for c in F_fuxi VM_fuxi; do
  lf $c noinf 1 L84_${c}_noinf_dpn_1g; lf $c ga3b 4 L84_${c}_ga3b_dpn; hp $c; hp ${c}_inf; done
for c in R_1km_5min F_graphcast F_aifs F_fuxi; do lf $c ga3b_sm 4 L84_${c}_ga3b_sm_dpn; hp ${c}_inf_sm; done
