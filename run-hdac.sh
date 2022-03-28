#!/bin/bash
project=hdac
cpu_per_item=30

find "./data/${project}" -mindepth 1 -maxdepth 1 -type d -print0 | while IFS= read -r -d '' item; do
item_result_dir="./data/${project}_results/$(basename "$item")"
mkdir -p "$item_result_dir"
cat <<EOT
#!/usr/bin/bash -l
#SBATCH -J protac-gao-${project}-$(basename "$item")-refine
#SBATCH -o ${item_result_dir}/%j-refine.log
#SBATCH -e ${item_result_dir}/%j-refine.log
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=${cpu_per_item}

echo "Started at \$(date)" >&2
__conda_init
source "env.super.sh"
python main.py -irec ${item}/receptor.pdb -ilig ${item}/target.pdb -site=$(head -n1 "${item}/site_info.txt") \
               -ismi ${item}/protac.smi -ie3lig1 ${item}/rec_lig_1.sdf -ie3lig2 ${item}/rec_lig_2.sdf \
               -o ${item_result_dir} -cpu ${cpu_per_item} -refine
echo "Finished at \$(date)" >&2
EOT
done
