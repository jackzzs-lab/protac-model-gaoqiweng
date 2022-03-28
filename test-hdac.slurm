#!/usr/bin/bash -l
#SBATCH -J protac-gao-hdac-HDAC7-TMP269-C6-Triazole-C2-CRBN-oO-refine
#SBATCH -o ./data/hdac_results/HDAC7-TMP269-C6-Triazole-C2-CRBN-oO/%j-refine.log
#SBATCH -e ./data/hdac_results/HDAC7-TMP269-C6-Triazole-C2-CRBN-oO/%j-refine.log
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=30

echo "Started at $(date)" >&2
__conda_init
source "env.super.sh"
python main.py -irec ./data/hdac/HDAC7-TMP269-C6-Triazole-C2-CRBN-oO/receptor.pdb -ilig ./data/hdac/HDAC7-TMP269-C6-Triazole-C2-CRBN-oO/target.pdb \
               -site=67.92,45.93,48.65 -ismi ./data/hdac/HDAC7-TMP269-C6-Triazole-C2-CRBN-oO/protac.smi \
               -ie3lig1 ./data/hdac/HDAC7-TMP269-C6-Triazole-C2-CRBN-oO/rec_lig_1.sdf -ie3lig2 ./data/hdac/HDAC7-TMP269-C6-Triazole-C2-CRBN-oO/rec_lig_2.sdf \
               -o ./data/hdac_results/HDAC7-TMP269-C6-Triazole-C2-CRBN-oO -cpu 30 -refine -refineonly
echo "Finished at $(date)" >&2
 