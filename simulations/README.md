This workflow runs the evolutionary simulations, which corresponds to the description in the supplementary materials and methods section, entitled 'evolutionary simulations,' found within the accompanying manuscript. 

# File Architecture
The project has a hierarchical structure that separates data, results, and scripts into distinct folders for efficient management and reproducibility. The primary components of this file architecture are:

1. [data](https://github.com/livkosterlitz/crowdsourcing/tree/main/simulations/data): Raw input files required for running the simulations such as allele-specific resistance levels and treatment master files. These are placed into two sub-folders for each gene used in the study (i.e., bla and DHFR). 
2. [results](https://github.com/livkosterlitz/crowdsourcing/tree/main/simulations/results): The ```.csv``` output files corresponding to each simulation numbered according to the numbers in the treatment master files. These are used to generate the figures in the manuscript. These outputs are placed into two sub-folders, one for each gene used in the study (i.e., bla and DHFR). 
3. [src](https://github.com/livkosterlitz/crowdsourcing/tree/main/simulations/src): Custom R scripts for running the evolutionary simulations.
4. [bin](https://github.com/livkosterlitz/crowdsourcing/tree/main/simulations/bin): Utility files for creating the Anaconda environment and the master ```runall.sh``` shell script.

# Program Requirements
We provided a YAML [```.yml``` ](https://github.com/livkosterlitz/crowdsourcing/blob/main/simulations/bin/crowdsourcing_simulations_environment.yml) file to be used with Anaconda to provide an environment for this specific project. It defines the packages and their specific versions needed to run this code smoothly. This makes it easy to replicate this workflow across different machines and platforms. 

```bash
conda env create -f bin/crowdsourcing_simulations_environment.yml
```

# Analysis Workflow
We have included a [```runall.sh```](https://github.com/livkosterlitz/crowdsourcing/blob/main/simulations/bin/runall.sh) script in our project to execute each simulation used in the paper according to the treatment files which specify the modeling parameters. Each simulation generates ```.csv``` files and are placed in the ```results``` sub-folder. 
