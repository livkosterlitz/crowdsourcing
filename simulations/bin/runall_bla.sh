source activate crowdSims
## evolutionary simulations for the bla
mkdir -p results/
mkdir -p results/bla
Rscript src/GradSelSim_bla.R
conda deactivate