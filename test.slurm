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
python main.py -irec ./data/crystal/6BN7/receptor.pdb -ilig ./data/crystal/6BN7/target.pdb -site=65.02,47.86,49.05 -ismi ./data/crystal/6BN7/protac.smi -o ./data/crystal_results/6BN7 -cpu 48 --skip-frodock
echo "Finished at $(date)" >&2
 