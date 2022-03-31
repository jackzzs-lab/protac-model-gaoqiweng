#!/usr/bin/bash -l
#SBATCH -J protac-gao-test
#SBATCH -o test.log
#SBATCH -e test.log
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=48

echo "Started at $(date)" >&2
__conda_init
source "env.super.sh"
python main.py -irec ./data/crystal/6BOY/receptor.pdb -ilig ./data/crystal/6BOY/target.pdb -site=65.05,47.07,48.65 -ismi ./data/crystal/6BOY/protac.smi -o ./data/crystal_results/6BOY -cpu 48 --skip-frodock
echo "Finished at $(date)" >&2
 