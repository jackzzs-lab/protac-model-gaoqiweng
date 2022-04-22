module load intel rosetta/3.13 autodock-vina/1.1.2 adfrsuite voronota fcc frodock/2.1 schrodinger
while [ "$CONDA_SHLVL" -ne 0 ]; do
    conda deactivate
done
conda activate main2
export ROSETTA="/opt/rosetta/3.13"
export ADFRSUITE="${ADFRSUITE_HOME}"
export FRODOCK="${FRODOCK_HOME}"
export VINA="${AUTODOCK_VINA_HOME}"
export VOROMQA="${VORONOTA_HOME}"
export FCC="${FCC_HOME}"
export MERGER="$HOME/proj/protac-ternary-modelling/gaoqiweng/utils/merger.py"
