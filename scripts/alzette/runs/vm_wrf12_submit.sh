#!/bin/bash
# Rescaled (VM) WRF 12 km before and after DA in both models: forcing -> LISFLOOD-FP inputs / HiPIMS setup -> runs -> extraction.
set -e
H=/ec/res4/scratch/lux0804/hipims; W=/ec/res4/scratch/lux0804/radar_scale_2026
f=$(sbatch --parsable --qos=nf --account=luuni --mem=64G --cpus-per-task=4 --time=02:00:00 --job-name=vmwrf_forcing --output=$W/logs/vmwrf_forcing_%j.log \
  --wrap="module load python3; python3 -u $W/code/make_forcing.py $W VM_wrf12 && python3 -u $W/code/make_forcing.py $W VM_wrf12noda")
i=$(sbatch --parsable --dependency=afterok:$f --qos=nf --account=luuni --mem=64G --cpus-per-task=4 --time=02:00:00 --job-name=vmwrf_inputs --output=$W/logs/vmwrf_inputs_%j.log \
  --wrap="module load python3; cd $H && python3 -u lf84_inputs_wrf.py VM_wrf12 VM_wrf12noda")
fv() { sbatch --dependency=afterok:$1 --job-name=fv --qos=nf --account=luuni --cpus-per-task=32 --mem=32G --time=01:00:00 \
  --output=$H/logs/fv_$2.out --wrap="module load python3; cd $H && python3 xs_fullvalley.py $3 $4" >/dev/null; }
for c in VM_wrf12 VM_wrf12noda; do
  j=$(sbatch --parsable --dependency=afterok:$i --job-name=l84_${c}_noinf --ntasks=1 --gres=gpu:1 $W/code/run_case84.sh $c noinf dp 1); fv $j L84_${c}_noinf_dp_1g lf $W/runs/L84_${c}_noinf_dp_1g
  k=$(sbatch --parsable --dependency=afterok:$i --job-name=l84_${c}_ga3b $W/code/run_case84.sh $c ga3b dp 4); fv $k L84_${c}_ga3b_dp lf $W/runs/L84_${c}_ga3b_dp
  echo "lf $c: $j $k"
  for cc in $c ${c}_inf; do
    s=$(sbatch --parsable --dependency=afterok:$f --job-name=setup_$cc --qos=nf --account=luuni --cpus-per-task=4 --mem=64G --time=02:00:00 --output=$H/logs/setup_${cc}_%j.log --wrap="cd $H && venv/bin/python hipims_setup.py $cc 4")
    r=$(sbatch --parsable --dependency=afterok:$s --job-name=hpmg_$cc $H/run_hipims_mg.sh ${cc}_mg4); fv $r ${cc}_mg4 hp $H/runs/${cc}_mg4
    echo "hp $cc: setup $s run $r"
  done
done
