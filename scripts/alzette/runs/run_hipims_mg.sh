#!/bin/bash
# usage: sbatch --job-name=hpmg_<case> run_hipims_mg.sh <case>_mg4   (HiPIMS multi-GPU solver, 4 GPUs)
#SBATCH --qos=ng
#SBATCH --account=luuni
#SBATCH --nodes=1 --gres=gpu:4 --ntasks=1 --cpus-per-task=16 --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/ec/res4/scratch/lux0804/hipims/logs/run_%x_%j.log
set -e
module load prgenv/gnu cuda/12.6 >/dev/null 2>&1
cd /ec/res4/scratch/lux0804/hipims/runs/$1
for d in */; do mkdir -p ${d}output; done
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
/usr/bin/time -v -o time.txt /ec/res4/scratch/lux0804/hipims/src/release/bin/hipims-flood-mgpus
