## Packages ----------------------------------------------------------------
library(argparse)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

## Command line ------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("-o", "--outputname", help="beginning of the file name")
parser$add_argument("-m", "--map", help="genotype to barcode map")
parser$add_argument("-f", "--file", help="sample cluster CSV output file")
args <- parser$parse_args()
output_filename <- args$outputname
filename_genotype_map <- args$map
filename_clusters <- args$file

## Troubleshooting variables ---------------------------------------------------------------
# filename_genotype_map <- "data/genotype_barcode_map.csv"
# filename_clusters <- "results/3_bartender_cluster/139_cluster.csv"
# output_filename <- "0"

## Load data ---------------------------------------------------------------
genotype_barcode_xl <- read_csv(filename_genotype_map)
clusters_csv <- read_csv(filename_clusters)

## step a: link genotype to barcode ------------------------------------------------
clusters <- clusters_csv %>% # make the barcode column. This is the reverse complement of the barcode in the map given that we are using read2 
  rename(ReverseBarcode = Center)
all_info <- full_join(genotype_barcode_xl, clusters)
write_csv(all_info, paste('results/4_barcode_frequency/a_linked_clusters/', output_filename, '_linked.csv', sep = ""))

## step b: frequencies --------------------------------------------------------
GenotypeCounts <- all_info %>%
  filter(!is.na(Genotype)) %>% 
  filter(!is.na(time_point_1)) %>%
  group_by(Shorthand) %>%
  mutate(BarcodeCounts = time_point_1,
         GenotypeCounts = sum(BarcodeCounts, na.rm = T))

BarcodeCountSummary <- GenotypeCounts %>%
  group_by(Shorthand) %>%
  summarise(N = n())

TotalCounts <- sum(GenotypeCounts$BarcodeCounts, na.rm = T)

JoinGenotypeCounts <- left_join(GenotypeCounts, BarcodeCountSummary)
JoinGenotypeCounts$TotalCounts <- TotalCounts

Frequencies <- JoinGenotypeCounts %>%
  mutate(GenotypeFrequency = GenotypeCounts/TotalCounts,
         BarcodeFrequency = BarcodeCounts/TotalCounts) %>%
  select(Genotype,
         Shorthand,
         Barcode,
         BarcodeCounts,
         GenotypeCounts,
         GenotypeFrequency,
         BarcodeFrequency,
         N)
genotype_barcode_xl_out <- genotype_barcode_xl %>%
  select(Genotype, Shorthand, Barcode)
Frequencies_out <- left_join(genotype_barcode_xl_out, Frequencies)
write_csv(Frequencies_out, paste('results/4_barcode_frequency/b_frequencies/csv_files/', output_filename, '_frequencies.csv', sep = ""))

## step b: diagnostic plots ---------------------------------------------------------------
pdf(paste('results/4_barcode_frequency/b_frequencies/plots/', output_filename, '_counts_vs_counts.pdf'), width = 7, height = 7)
ggplot(data = Frequencies, aes(x = GenotypeCounts, y = BarcodeCounts, color = Shorthand))+
       geom_point()
dev.off()

pdf(paste('results/4_barcode_frequency/b_frequencies/plots/', output_filename, '_freq_vs_freq.pdf'), width = 7, height = 7)
ggplot(data = Frequencies, aes(x = GenotypeFrequency, y = BarcodeFrequency, color = Shorthand))+
  geom_point()
dev.off()

pdf(paste('results/4_barcode_frequency/b_frequencies/plots/', output_filename, '_genotypes_by_barcode_proportion.pdf'), width = 15, height = 6)
ggplot(data = Frequencies, aes(x = Shorthand, y = BarcodeCounts, fill = Barcode)) +
  geom_histogram(stat = 'identity', position = 'fill', color = 'black') +
  theme(legend.text = element_text(size = 4), legend.title = element_text(size = 4))+
  theme(legend.position = "none")
dev.off()


