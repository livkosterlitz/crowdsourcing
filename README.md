# Introduction

This repository accompanies the publication of **_Evolutionary "crowdsourcing" allows for cross-species adaption of a horizontally transferred gene_** by Kosterlitz et. al. This research project explores the cross-species adaptation of a horizontally transferred gene through evolutionary "crowdsourcing." The repository provides all relevant data, code, and figures associated with the publication, enabling users to replicate the results and explore the findings in-depth.

# Dependencies 
Dependencies for each part of the project are described in their respective subfolders. To install the required dependencies, locate the provided YAML file in the appropriate subfolder and follow the instructions for installation.

For example, if a subfolder contains a dependencies.yml file, you can install the dependencies using the following command:

```bash
conda env create -f dependencies.yml
```

This will create a new Conda environment with the necessary dependencies installed. Make sure to activate the environment before running any code:

```bash
conda activate <env_name>
```

Replace ```<env_name>``` with the name of the environment specified in the YAML file.

By following the instructions in each subfolder, you can ensure that all required dependencies are installed correctly before running the code.

# Usage
To replicate the results of the study, follow the steps provided in the respective folders for pooled competition analysis and simulations. Detailed instructions and workflows can be found within each folder.

# Navigating the repository

## **Data and Figures** 
The [```figures```](https://github.com/livkosterlitz/Crowdsourcing) folder contains all the raw data, calculations, codes, and final figures includued in the publication. To access the associated data and materials, refer to the figure numbers as ordered in the paper.

## **Pooled competition analysis to construct host-specific landscapes**
The [```competition_analysis```](https://github.com/livkosterlitz/crowdsourcing/tree/main/competition_analysis) folder contains the code, workflow, and data generated from analyzing the data for the pooled competition assays used to construct host-specific landscapes.

## **Simulations**
The [```simulations```](https://github.com/livkosterlitz/Crowdsourcing) folder contains the code, workflow, and data generated from the evolutionary simulations. 
