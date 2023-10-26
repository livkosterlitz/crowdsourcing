source activate crowdSims
### download the new .yml file and put it in the bin folder

## step 1: iterate through random seeds & calculate alignment metrics
## checked 2023-08-15 by OK
# mkdir -p results/
# mkdir -p results/random
# mkdir -p results/random/1_random_seeds
# Rscript src/Seeds_random.R > /dev/null 2>&1
# echo "step 1 finished"

## step 2: filter and select seeds
## checked 2023-08-15 by OK
# mkdir -p results/random/2_selected_seeds
# Rscript src/Seeds_select.R > /dev/null 2>&1
# echo "step 2 finished"

## step 3: run simulations in parallel
mkdir -p results/random/3_simulations
chmod a+x results/random/2_selected_seeds/step3_runall.sh
time results/random/2_selected_seeds/step3_runall.sh
echo "step 3 finished"

## step 4: collect simulations and determine evolutionary outcomes

## step 5: Permutation test to determine signficance between the landscape alignment and evolutionary outcomes
#Rscript src/GradSelSim_bla.R -i step3_output.csv
# BK job
conda deactivate