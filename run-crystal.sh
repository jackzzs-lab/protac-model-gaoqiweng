#!/bin/bash
project=crystal

find "./data/${project}" -mindepth 1 -maxdepth 1 -type d -not -name 'raw' -print0 | while IFS= read -r -d '' item; do
    if [ ! -f "${item}/protac.smi" ]; then
        continue
    fi
    item_result_dir="./data/${project}_results/$(basename "$item")"
    mkdir -p "$item_result_dir"
    if [ ! -d "${item_result_dir}/frodock_results" ]; then
            sbatch <<EOT
#!/usr/bin/bash -l
#SBATCH -J protac-gao-${project}-$(basename "$item")
#SBATCH -o ${item_result_dir}/%j.log
#SBATCH -e ${item_result_dir}/%j.log
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=8

echo "Started at \$(date)" >&2
__conda_init
source "env.super.sh"
python main.py -irec ${item}/receptor.pdb -ilig ${item}/target.pdb -site=$(head -n1 "${item}/site_info.txt") \
               -ismi ${item}/protac.smi -o ${item_result_dir} -cpu 9 -refine
echo "Finished at \$(date)" >&2
EOT
    fi
done
