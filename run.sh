#!/bin/bash

help_msg() {
    echo "Usage: $(basename "$0") [-c <cpus>] [-cr <cpus>] [-p <project>] [-r/-ro/-nf] [items]" >&2
    echo "       -c/--cpu: cpu to allocate when submitting, default=6." >&2
    echo "       -cr/--cpu-run: cpu to use when running, default=cpu+2." >&2
    echo "       -d/--dry: dry run only, run cat instead of sbatch." >&2
    echo "       -p/--project: project path containing multiple items." >&2
    echo "       -r/--refine: run rosetta phase." >&2
    echo "       -ro/--refine-only: skip frodock and frodock-clustering phase." >&2
    echo "       -nf/--skip-frodock: skip frodock phase." >&2
    echo "Note:  the protocol can run frodock->frodock-clustering->rosetta," >&2
    echo "       by default only frodock->frodock-clustering will be run." >&2
}

ITEMS=()
PARAMS=()
while [[ $# -gt 0 ]]; do
    case $1 in
    -c | --cpu)
        CPUS="$2"
        shift # past argument
        shift # past value
        ;;
    -cs | --cpu-submit)
        CPUS_RUN="$2"
        shift # past argument
        shift # past value
        ;;
    -d | --dry)
        DRY="YES"
        shift # past argument
        ;;
    -p | --project)
        PROJECT="$2"
        shift # past argument
        shift # past value
        ;;
    -nf | --skip-frodock)
        FRODOCK_REQ="YES"
        FRODOCK_SKIP="NO"
        PARAMS+=("--skip-frodock")
        shift # past argument
        ;;
    -ro | --refine-only)
        FRODOCK_REQ="YES"
        PARAMS+=("--rosettadock-refinement-only")
        shift # past argument
        ;;
    -r | --refine)
        FRODOCK_SKIP="NO"
        ROSETTA_SKIP="YES"
        PARAMS+=("--rosettadock-refinement")
        shift # past argument
        ;;
    -h | --help)
        help_msg
        exit 0
        ;;
    -*)
        echo "Unknown option $1"
        help_msg
        exit 1
        ;;
    *)
        ITEMS+=("$1") # save positional arg
        shift         # past argument
        ;;
    esac
done

set -- "${ITEMS[@]}" # restore positional parameters

[ -n "$DRY" ] && cmd="cat" || cmd="sbatch"
[ -z "$CPUS" ] && CPUS=6
[ -z "$CPUS_RUN" ] && CPUS_RUN="$((CPUS + 2))"
if [ -d "$PROJECT" ]; then
    while IFS= read -r -d '' item; do
        ITEMS+=("$item")
    done < <(find "$PROJECT" -mindepth 1 -maxdepth 1 -type d -print0)
fi

for item in "${ITEMS[@]}"; do
    if [ ! -f "${item}/protac.smi" ] || [ ! -f "${item}/receptor.pdb" ] || [ ! -f "${item}/target.pdb" ]; then
        echo "Warning: '${item}' is not ready to be run, skipping." >&2
        continue
    fi
    project="$(basename "$(dirname "$item")")"
    item_result_dir="$(dirname "$item")_results/$(basename "$item")"
    mkdir -p "$item_result_dir"
    if [ -z "$FRODOCK_SKIP" ] && [ -d "${item_result_dir}/frodock_results" ]; then
        continue
    elif [ -n "$ROSETTA_SKIP" ] && [ -d "${item_result_dir}/rosetta_results" ]; then
        continue
    fi
    if [ -n "$FRODOCK_REQ" ] && [ ! -d "${item_result_dir}/frodock_results" ]; then
        echo "Warning: '${item}' requires frodock_results, which does not exist, skipping." >&2
        continue
    fi
    echo "Info: submitting '${item}'." >&2
    [ -n "$DRY" ] && echo "============================================================="
    "$cmd" <<EOT
#!/usr/bin/bash -l
#SBATCH -J protac-gao-${project}-$(basename "$item")
#SBATCH -o ${item_result_dir}/%j.log
#SBATCH -e ${item_result_dir}/%j.log
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --cpus-per-task=${CPUS}

echo "Started at \$(date)" >&2
__conda_init
source "env.super.sh"
python main.py -irec ${item}/receptor.pdb -ilig ${item}/target.pdb -site=$(head -n1 "${item}/site_info.txt") \
-ismi ${item}/protac.smi -o ${item_result_dir} -cpu ${CPUS_RUN} $(IFS=' ' ;echo "${PARAMS[*]}") \
&>>"${item_result_dir}/\${SLURM_JOB_ID}.log"
echo "Finished at \$(date)" >&2
EOT
[ -n "$DRY" ] && echo "============================================================="
done
