#!/bin/bash
# Ks sensitivity: radar and GraphCast, Ks x0.5 and x2, LISFLOOD-FP (4 GPUs, dp) and HiPIMS (4 GPUs); each run chained to its extraction.
set -e
H=/ec/res4/scratch/lux0804/hipims; W=/ec/res4/scratch/lux0804/radar_scale_2026
for c in R_1km_5min F_graphcast; do for k in ks05 ks2; do
  j=$(sbatch --parsable --job-name=l84_${c}_ga3b_${k} $W/code/run_case84.sh $c ga3b_$k dp 4)
  sbatch --dependency=afterok:$j --job-name=fv --qos=nf --account=luuni --cpus-per-task=32 --mem=32G --time=01:00:00 \
    --output=$H/logs/fv_L84_${c}_ga3b_${k}_dp.out --wrap="module load python3; cd $H && python3 xs_fullvalley.py lf $W/runs/L84_${c}_ga3b_${k}_dp" >/dev/null
  s=$(sbatch --parsable --job-name=setup_${c}_inf_${k} --qos=nf --account=luuni --cpus-per-task=4 --mem=64G --time=02:00:00 \
    --output=$H/logs/setup_${c}_inf_${k}_%j.log --wrap="cd $H && venv/bin/python hipims_setup.py ${c}_inf_${k} 4")
  r=$(sbatch --parsable --dependency=afterok:$s --job-name=hpmg_${c}_inf_${k} $H/run_hipims_mg.sh ${c}_inf_${k}_mg4)
  sbatch --dependency=afterok:$r --job-name=fv --qos=nf --account=luuni --cpus-per-task=32 --mem=32G --time=01:00:00 \
    --output=$H/logs/fv_${c}_inf_${k}_mg4.out --wrap="module load python3; cd $H && python3 xs_fullvalley.py hp $H/runs/${c}_inf_${k}_mg4" >/dev/null
  echo "$c $k: lf=$j setup=$s hp=$r"
done; done
