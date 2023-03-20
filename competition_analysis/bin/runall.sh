cd /Users/oliviakosterlitz/Dropbox/preTUMS/Competition_analysis
source activate crowd
## step 1: add read prefix
mkdir -p results/
mkdir -p results/1_sequencing_prefix
find "./data/sequencing" -type f -name "*.fastq" -exec cp {} "./results/1_sequencing_prefix" \;
find "./results/1_sequencing_prefix" -type f -name "*.fastq" | xargs sed -i.bak 'n;s/^/AAAA/'
## step 2: barcode extraction
mkdir -p results/2_bartender_barcode_extract
for FILE in results/1_sequencing_prefix/*.fastq; 
do 
    FILENAME="${FILE##*/}"
    TREATMENT="${FILENAME/_*.fastq/}"
    bartender_extractor_com -f $FILE -o results/2_bartender_barcode_extract/${TREATMENT} -q ? -p AAAA[18]CGTA -m 2 -d both
done
## step 3: barcode clustering
mkdir -p results/3_bartender_cluster
for FILE in results/2_bartender_barcode_extract/*.txt; 
do 
    FILENAME="${FILE##*/}"
    TREATMENT="${FILENAME/_*.txt/}"
    bartender_single_com -f $FILE -o results/3_bartender_cluster/${TREATMENT} -d 3
done
## step 4: calculate barcode frequency
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
## step 5: barcode approximate growth rate 
mkdir -p results/5_approximate_growth_rate/
mkdir -p results/5_approximate_growth_rate/a_Cmax
mkdir -p results/5_approximate_growth_rate/b_Tc
mkdir -p results/5_approximate_growth_rate/c_AppGrowthRate
mkdir -p results/5_approximate_growth_rate/d_Outlier
Rscript src/treatment_analysis.R -t data/treatment_master.xlsx -f results/4_barcode_frequency/b_frequencies/csv_files/ -m data/genotype_barcode_map.csv
## step 6: determine each genotypes level of resistance
mkdir -p results/6_dose_curves/
mkdir -p results/6_dose_curves/a_lowerAsymptote
mkdir -p results/6_dose_curves/b_modelFit
mkdir -p results/6_dose_curves/c_resistance
mkdir -p results/6_dose_curves/d_mutations
Rscript src/curve_analysis.R -g results/5_approximate_growth_rate/d_Outlier/Outliers.csv -b data/genotype_coding.csv -s data/mutational_steps.csv
conda deactivate