This collection of folders are named according to the figure numbers in the manuscript and contain the data files and R scripts used to create the figures found within the accompanying manuscript. 

# Program Requirements
We provided a YAML [```.yml``` ](https://github.com/livkosterlitz/crowdsourcing/blob/main/figures/bin/crowdsourcing_figures_environment.yml) file to be used with Anaconda to define the packages and their specific versions needed to run the figure codes smoothly. This allows for easy replication of this workflow across different machines and platforms.

```bash
conda env create -f bin/crowdsourcing_figures_environment.yml
```

# Analysis Workflow
To regenerate the figures, first activate the environment. 

```bash
conda activate crowdFigs
```

Lastly, execute the appropriate Rscript. 

```bash
Rscript Fig.R
```
