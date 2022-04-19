#!/bin/bash
cpu_per_item=40

for item in "./data/example_2"; do
item_name="$(basename "$item")"
item_result_dir="./data/$(basename "$item")_results"
mkdir -p "$item_result_dir"
sbatch <<EOT
#!/usr/bin/bash -l
#SBATCH -J protac-gao-${item_name}
#SBATCH -o ${item_result_dir}/%j.log
#SBATCH -e ${item_result_dir}/%j.log
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=${cpu_per_item}

echo "Started at \$(date)" >&2
__conda_init
source "env.super.sh"

if [[ "${item_name}" == "example_1" ]]; then
    python main.py -irec ${item}/receptor.pdb -ilig ${item}/target.pdb -site=$(head -n1 "${item}/site_info.txt") \
                -ismi ${item}/protac.smi -o ${item_result_dir} -cpu ${cpu_per_item} -refine
elif [[ "${item_name}" == "example_2" ]]; then
    python main.py -irec ${item}/receptor.pdb -ilig ${item}/target.pdb -site=$(head -n1 "${item}/site_info.txt") \
                -ismi ${item}/protac.smi -ie3lig1 ${item}/rec_lig_1.sdf -ie3lig2 ${item}/rec_lig_2.sdf \
                -o ${item_result_dir} -cpu ${cpu_per_item} -refine
fi
echo "Finished at \$(date)" >&2
EOT
done
