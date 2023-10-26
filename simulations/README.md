This workflow runs the evolutionary simulations, which corresponds to the description in the supplementary materials and methods section, entitled 'evolutionary simulations,' found within the accompanying manuscript. 

# File Architecture
The project has a hierarchical structure that separates data, results, and scripts into distinct folders for efficient management and reproducibility. The primary components of this file architecture are:

1. [data](https://github.com/livkosterlitz/crowdsourcing/tree/main/simulations/data): Raw input files required for running the simulations such as allele-specific resistance levels and treatment master files. These are placed into two sub-folders for each gene used in the study (i.e., bla and random). 
2. [results](https://github.com/livkosterlitz/crowdsourcing/tree/main/simulations/results): The ```.csv``` output files corresponding to each simulation are numbered according to the numbers in the treatment master files. These are used to generate the figures in the manuscript. These outputs are placed into two sub-folders, one for each gene used in the study (i.e., bla and DHFR). 
3. [src](https://github.com/livkosterlitz/crowdsourcing/tree/main/simulations/src): Custom R scripts for running the evolutionary simulations.
4. [bin](https://github.com/livkosterlitz/crowdsourcing/tree/main/simulations/bin): Utility files for creating the Anaconda environment and the master ```runall.sh``` shell script.

# Program Requirements
We provided a YAML [```.yml``` ](https://github.com/livkosterlitz/crowdsourcing/blob/main/simulations/bin/crowdsourcing_simulations_environment.yml) file to be used with Anaconda to provide an environment for this specific project. It defines the packages and their specific versions needed to run this code smoothly. This makes it easy to replicate this workflow across different machines and platforms. 

```bash
conda env create -f bin/crowdsourcing_simulations_environment.yml
```

# Analysis Workflow for evolutionary simulations on empirical landscapes
For the evolutionary simulations on the empirical landscapes, we have included a [```runall_bla.sh```](https://github.com/livkosterlitz/crowdsourcing/blob/main/simulations/bin/runall_bla.sh) script in our project to execute each simulation used in the paper according to the treatment files which specify the modeling parameters. Each simulation generates ```.csv``` files and is placed in the ```results/bla``` sub-folder. 


# Analysis Workflow for evolutionary simulations on randomly generated landscapes
For the evolutionary simulations on randomly generated landscapes, we have included a [```runall_random.sh```](https://github.com/livkosterlitz/crowdsourcing/blob/main/competition_analysis/bin/runall_random.sh) script in our project to streamline the process of executing these simulations. Each step (outlined in the table below) has a corresponding sub-folder placed in the ```results/random``` sub-folder. 

|Steps| Step description |
| :--- | :--- | 
| [Step 1](#Step-1) | Iterate through seeds to create randomly generated landscapes and calculate alignment metrics| 
| [Step 2](#Step-2) | Filter and select the seeds| 
| [Step 3](#Step-3) | Run the sims for each selected seed in parallel |
| [Step 4](#Step-4) | Collect simulations and determine evolutionary outcomes |
| [Step 5](#Step-5) | Permutation test to determine the significance between the landscape alignment and evolutionary outcome|
