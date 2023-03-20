This workflow examines the pooled competition assays to create species-specific landscapes. This analysis corresponds to the supplementary materials and methods section, entitled 'library sequence analysis,' found within the accompanying manuscript.

# File Architecture
The project has a hierarchical structure that separates data, results, and scripts into distinct folders for efficient management and reproducibility. The primary components of this file architecture are:

1. [data](https://github.com/livkosterlitz/crowdsourcing/tree/main/competition_analysis/data): Raw input files required for the analysis such as decompressed fastq files, genotype-barcode maps, and treatment master files.
2. [results](https://github.com/livkosterlitz/crowdsourcing/tree/main/competition_analysis/results): Sub-folders corresponding to each step of the analysis workflow. Each sub-folder stores the output files (e.g., .csv and .pdf) generated at that specific stage.
3. [src](https://github.com/livkosterlitz/crowdsourcing/tree/main/competition_analysis/src): Custom scripts (e.g., R scripts) used throughout the analysis workflow.
4. [bin](https://github.com/livkosterlitz/crowdsourcing/tree/main/competition_analysis/bin): Utility files for creating the Anaconda environment and the master shell script.

# Program Requirements
We provided a YAML [```.yml``` ](https://github.com/livkosterlitz/crowdsourcing/blob/main/competition_analysis/bin/crowdsourcing_environment.yml) file to be used with Anaconda to provide an environment for this specific project. It defines the packages and their specific versions needed to run this code smoothly. This makes it easy to replicate this workflow across different machines and platforms. 

```bash
conda env create -f bin/crowdsourcing_environment.yml
```

# Analysis Workflow
We have included a [```runall.sh```](https://github.com/livkosterlitz/crowdsourcing/blob/main/competition_analysis/bin/runall.sh) script in our project to streamline the process of executing multiple commands and scripts in a specified sequence. The workflow consists of six main steps, each with various sub-steps denoted by alphabetical letters. Throughout the execution of these steps, the script generates a range of outputs, including ```.csv``` files and graphical plots saves as ```.pdf``` files. This automation helps ensure a consistent and efficient workflow, making it easier to manage and reproduce the project's outcomes. 

|Steps| Step description |
| :--- | :--- | 
| [Step 1](#Step-1)  | Add a 4 character prefix to all the reads in each uncompressed fastq file. A necessary step for the downstream bartender workflow.| 
| [Step 2](#Step-2) | Extract the 18bp barcode from each sequenced read from the FASTQ file using Bartender.| 
| [Step 3](#Step-3) | Cluster the 18bp barcodes into groups using Bartender.|
| [Step 4](#Step-4) | Calculate barcode frequency by dividing the barcode counts by using the total counts in each sample.|
| [Step 5](#Step-5) | Calculate the approximate growth rate for each barcode using the equation derived in supplemental section 1 of the accompanying manuscript.|
| [Step 6](#Step-6) | Identify the level of resistance of each allele by fitting a four-parameter log logistic dose-response curve.|  

## Step 1: Add read prefix
<a id="Step-1"></a>

Utilize the decompressed [ ```.fastq``` ](https://github.com/livkosterlitz/crowdsourcing/tree/main/competition_analysis/data/sequencing) files to add a 4 character prefix to all of the reads with the following bash command. This step is essential for the subsequent bartender workflow. 

```bash
mkdir -p results/1_sequencing_prefix
find "./data/sequencing" -type f -name "*.fastq" -exec cp {} "./results/1_sequencing_prefix" \;
find "./results/1_sequencing_prefix" -type f -name "*.fastq" | xargs sed -i.bak 'n;s/^/AAAA/'
```

## Step 2: Barcode extraction
<a id="Step-2"></a>

For each output ```.fastq``` file, extract the barcodes from each read using the barcode extractor from [bartender](https://github.com/LaoZZZZZ/bartender-1.1) with the bash command provided below. 

```bash
mkdir -p results/2_bartender_barcode_extract
for FILE in results/1_sequencing_prefix/*.fastq; 
do 
    FILENAME="${FILE##*/}"
    TREATMENT="${FILENAME/_*.fastq/}"
    bartender_extractor_com -f $FILE -o results/2_bartender_barcode_extract/${TREATMENT} -q ? -p AAAA[18]CGTA -m 2 -d both
done
```
## Step 3: Barcode clustering
<a id="Step-3"></a>

For each output ```.txt``` file from [Step 2](#Step-2), cluster the barcodes into groups using the barcode clustering from [bartender](https://github.com/LaoZZZZZ/bartender-1.1) using the bash command provided below. 

```bash
mkdir -p results/3_bartender_cluster
for FILE in results/2_bartender_barcode_extract/*.txt; 
do 
    FILENAME="${FILE##*/}"
    TREATMENT="${FILENAME/_*.txt/}"
    bartender_single_com -f $FILE -o results/3_bartender_cluster/${TREATMENT} -d 3
done
```
## Step 4: Calculate barcode frequency
<a id="Step-4"></a>

For each output ```cluster.csv``` file from [Step 3](#Step-3), analyze the clusters using the custom [```cluster_analysis.R```](https://github.com/livkosterlitz/crowdsourcing/blob/main/competition_analysis/src/cluster_analysis.R) script provided. The script performs the following steps: 
  * step a: Link the treatment's barcode clusters to the corresponding allele (e.g., single, double, and triple mutants, etc.) using the [genotype-barcode map](https://github.com/livkosterlitz/crowdsourcing/blob/main/competition_analysis/data/genotype_barcode_map.csv) verified through Sanger sequencing after plasmid cloning. 
  * step b: Calculate the barcode frequency by summing the total barcode counts in the treatment, then divides each barcode count by the total number of counts. 

```bash
mkdir -p results/4_barcode_frequency
mkdir -p results/4_barcode_frequency/a_linked_clusters
mkdir -p results/4_barcode_frequency/b_frequencies/csv_files
mkdir -p results/4_barcode_frequency/b_frequencies/plots
for FILE in results/3_bartender_cluster/*_cluster.csv; 
do 
    FILENAME="${FILE##*/}"
    TREATMENT="${FILENAME/_*.csv/}"
    Rscript src/cluster_analysis.R -o $TREATMENT -f $FILE -m data/genotype_barcode_map.csv
done
```

## Step 5: Calculate each barcode's approximate growth rate
<a id="Step-5"></a>

For each frequency output ```.csv``` file from [Step 4](#Step-4), calculate the approximate growth rate for each barcode using the custom [```treatment_analysis.R```](https://github.com/livkosterlitz/crowdsourcing/blob/main/competition_analysis/src/treatment_analysis.R) script provided. The script performs the following steps: 
 * step a: Determine the species-specific Cmax and the median approximate growth rate for Gprime at Cmax.
 * step b: Calculate the median Tc for Gprime at each concentration.
 * step c: Calculate the approximate growth rate for each barcode using the concentration-specific Tc.
 * step d: Identify and eliminate the most deviant barcode for each allele. 

```bash
mkdir -p results/5_approximate_growth_rate/
mkdir -p results/5_approximate_growth_rate/a_Cmax
mkdir -p results/5_approximate_growth_rate/b_Tc
mkdir -p results/5_approximate_growth_rate/c_AppGrowthRate
mkdir -p results/5_approximate_growth_rate/d_Outlier
Rscript src/treatment_analysis.R -t data/treatment_master.xlsx -f results/4_barcode_frequency/b_frequencies/csv_files/ -m data/genotype_barcode_map.csv
```
## Step 6: Identify the level of resistance for each allele
<a id="Step-6"></a>

Using the approximate growth rates for each barcode across the drug gradient, fit a four-parameter log logistic dose-response curve using the custom [```curve_analysis.R```](https://github.com/livkosterlitz/crowdsourcing/blob/main/competition_analysis/src/cluster_analysis.R) script provided. The script performs the following steps:

 * step a: Determine the average lower asymptote for each species by averaging the fitted lower asymptotes from the alleles with the lowest levels of resistance.
 * step b: Fit dose-response curves using the species-specific lower asymptotes, and create visual representations of the fitted curves.
 * step c: Calculate the level of resistance for each allele-species combination by averaging the inflection point parameters from each barcode, and report the resistance levels in two CSV output files. Both of these files are used in various figures in the accompanying manuscript. 
 * step d: Determine if a mutational step was beneficial, deleterious, or neutral by comparing the level of resistance for each neighboring allele, and report the differences between neighboring alleles in two CSV output files. Both of these files are used in various figures in the accompanying manuscript. 

```bash
mkdir -p results/6_dose_curves/
mkdir -p results/6_dose_curves/a_lowerAsymptote
mkdir -p results/6_dose_curves/b_modelFit
mkdir -p results/6_dose_curves/c_resistance
mkdir -p results/6_dose_curves/d_mutations
Rscript src/curve_analysis.R -g results/5_approximate_growth_rate/d_Outlier/Outliers.csv -b data/genotype_coding.csv -s data/mutational_steps.csv
```
