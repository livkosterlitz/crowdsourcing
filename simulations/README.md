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
For the evolutionary simulations on randomly generated landscapes, we have included a [```runall_random.sh```](https://github.com/livkosterlitz/crowdsourcing/blob/main/competition_analysis/bin/runall_random.sh) script in our project to streamline the process of executing ......

coming soon

|Steps| Step description |
| :--- | :--- | 
| [Step 1](#Step-1) | Iterate through seeds to create randomly generated landscapes and calculate alignment metrics| 
| [Step 2](#Step-2) | Filter and select the seeds| 
| [Step 3](#Step-3) | Run the sims for each selected seed in parallel |
| [Step 4](#Step-4) | |
| [Step 5](#Step-5) | |

## Step 1: Determine random seeds
<a id="Step-1"></a>

Utilize the decompressed [ ```.fastq``` ](https://github.com/livkosterlitz/crowdsourcing/tree/main/competition_analysis/data/sequencing) files to add a 4 character prefix to all of the reads with the following bash command. This step is essential for the subsequent bartender workflow. 

```bash
mkdir -p results/1_sequencing_prefix
find "./data/sequencing" -type f -name "*.fastq" -exec cp {} "./results/1_sequencing_prefix" \;
find "./results/1_sequencing_prefix" -type f -name "*.fastq" | xargs sed -i.bak 'n;s/^/AAAA/'
```

## Step 2: Run simulations
<a id="Step-2"></a>

For each seed, if the simulation does not have an output file indicating it already ran in the ```results/random/``` sub-folder (a check in the bash command below), then it runs the evolutionary simulation with that seed to generate the two host landscapes. 

uses the  extract the barcodes from each read using the barcode extractor from [bartender](https://github.com/LaoZZZZZ/bartender-1.1) with the bash command provided below. 

```bash
mkdir -p results/2_bartender_barcode_extract
for FILE in results/1_sequencing_prefix/*.fastq; 
do 
    FILENAME="${FILE##*/}"
    TREATMENT="${FILENAME/_*.fastq/}"
    bartender_extractor_com -f $FILE -o results/2_bartender_barcode_extract/${TREATMENT} -q ? -p AAAA[18]CGTA -m 2 -d both
done
```
## Step 3: ...
<a id="Step-3"></a>


## Step 4: ...
<a id="Step-4"></a>


## Step 5: ...
<a id="Step-5"></a>
