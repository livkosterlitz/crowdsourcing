source activate crowdSims
## evolutionary simulations for the bla
mkdir -p results/
mkdir -p results/bla
Rscript src/GradSelSim_bla.R
## evolutionary simulations for the DHFR
mkdir -p results/DHFR
Rscript src/GradSelSim_DHFR.R
conda deactivate