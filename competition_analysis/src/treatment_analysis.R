## Packages ----------------------------------------------------------------
library(argparse)
library(tidyverse)
library(broom)
library(readxl)
library(cowplot)
theme_set(theme_cowplot())
## Command line ------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("-t", "--treatmentsfile", help="the file with treatment information")
parser$add_argument("-m", "--map", help="genotype to barcode map")
parser$add_argument("-f", "--file", help="frequency CSV output file")
args <- parser$parse_args()
treatmentsfile <- args$treatmentsfile
filename_genotype_map <- args$map
frequency_csvfiles <- args$file

## Troubleshooting variables ---------------------------------------------------------------
# treatmentsfile <- "data/treatment_master.xlsx"
# filename_genotype_map <- "data/genotype_barcode_map.csv"
# frequency_csvfiles <- "results/4_barcode_frequency/b_frequencies/csv_files/"

## Treatment information ----------------------------------------------------------
treatment_description_file <- read_xlsx(treatmentsfile)

## Function Load data ---------------------------------------------------------------
Load_data <- function(datafolder, treatment_description) {
  datafiles <- list.files(datafolder, full.names = T)
  datafileslist <- lapply(datafiles, read.csv, header = T, stringsAsFactors = F)
  treatments <- strsplit(datafiles, "/")
  treatment_names <- rep(NA, length(treatments))
  for (i in 1:length(treatments)) {
    filename <- treatments[[i]][length(treatments[[i]])]
    treatment_names[i] <- strsplit(filename, '_')[[1]][1]
  }
  for (l in 1:length(datafileslist)) {
    datafileslist[[l]]['Treatment'] <- rep(treatment_names[l], nrow(datafileslist[[l]]))
    datafileslist[[l]]['Treatment_num'] <- rep(treatment_description['Treatment_num'][treatment_description['Treatment_num'] == treatment_names[l]], nrow(datafileslist[[l]]))
    datafileslist[[l]]['Species'] <- rep(treatment_description['Species'][treatment_description['Treatment_num'] == treatment_names[l]], nrow(datafileslist[[l]]))
    datafileslist[[l]]['Time'] <- rep(treatment_description['Time'][treatment_description['Treatment_num'] == treatment_names[l]], nrow(datafileslist[[l]]))
    #datafileslist[[l]]['Volume'] <- rep(treatment_description['Volume'][treatment_description['Treatment_num'] == treatment_names[l]], nrow(datafileslist[[l]]))
    datafileslist[[l]]['Concentration'] <- rep(treatment_description['Concentration'][treatment_description['Treatment_num'] == treatment_names[l]], nrow(datafileslist[[l]]))
    datafileslist[[l]]['Density'] <- rep(treatment_description['24h_density'][treatment_description['Treatment_num'] == treatment_names[l]], nrow(datafileslist[[l]]))
  }
  all_data <- bind_rows(datafileslist)
  all_data <- all_data %>%
    mutate(DensityPerBarcode = as.numeric(Density) * as.numeric(BarcodeFrequency))
  return(all_data)
}

## Load data ---------------------------------------------------------------
data <- Load_data(datafolder = frequency_csvfiles, treatment_description_file)
alldata <- data %>%
  select(Species, Genotype, Shorthand, Barcode, Time, Concentration, BarcodeCounts, BarcodeFrequency, Density, DensityPerBarcode)

## Restructure data  -------------------------------------------------------
FinalFrequency <- alldata %>% filter(Time != 0) %>% 
  select(-Time) %>%
  group_by(Species, Concentration) %>%
  mutate(Reads = sum(BarcodeCounts, na.rm = T)) %>%
  ungroup()
InitialFrequency <- alldata %>%
  filter(Time == 0) %>%
  select(Species, Shorthand, Barcode, BarcodeFrequency, Density, DensityPerBarcode) %>%
  rename(InitialBarcodeFrequency = BarcodeFrequency, 
         InitialDensity = Density,
         InitialDensityPerBarcode = DensityPerBarcode)
Frequency_Initial_Final <- left_join(FinalFrequency, InitialFrequency)
Frequency_Initial_Final$Density <- as.numeric(Frequency_Initial_Final$Density)
Frequency_Initial_Final$InitialDensity <- as.numeric(Frequency_Initial_Final$InitialDensity)

## step a: determine species-specific Cmax  -------------------------------------------------------
# calculate the approximate growth rate for Gprime at each tested concentration
Approx_growth_rate_Gprime <- Frequency_Initial_Final %>%
  filter(Shorthand == 'pAEGM') %>%
  mutate(approx_growth_rate = (1/24)*(log(Density/InitialDensity)+log(BarcodeFrequency/InitialBarcodeFrequency))) 
# find the concentration for each three barcodes for the best allele of each species where the approximate growth rate is maximal
Approx_growth_rate_Gprime_Cmax_byBarcode <- Approx_growth_rate_Gprime %>%
  group_by(Species, Barcode)  %>%
  filter(approx_growth_rate == max(approx_growth_rate)) %>%
  select(Species, Barcode, Concentration, approx_growth_rate)
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]}
# pick one concentration (i.e., the Cmax) from the three barcodes maximum concentrations. The one picked is either the most common concentration among the three or if they are all unique it is the minimum concentration. 
Approx_growth_rate_Gprime_Cmax <- Approx_growth_rate_Gprime_Cmax_byBarcode %>%
  group_by(Species)  %>%
  filter(if (n_distinct(Concentration) == 1) TRUE else Concentration == mode(Concentration)) %>%
  summarise(Cmax = if (n_distinct(Concentration) == 1) min(Concentration) else mode(Concentration))
Approx_growth_rate_Gprime <- left_join(Approx_growth_rate_Gprime, Approx_growth_rate_Gprime_Cmax)
Approx_growth_rate_Gprime_Cmax_median <- Approx_growth_rate_Gprime %>%
  filter(Concentration == Cmax) %>%
  group_by(Species, Cmax) %>%
  summarise(Gprime_maxC_growthrate = median(approx_growth_rate))
Approx_growth_rate_Gprime <- left_join(Approx_growth_rate_Gprime, Approx_growth_rate_Gprime_Cmax_median)
write_csv(Approx_growth_rate_Gprime, paste('results/5_approximate_growth_rate/a_Cmax/', 'Gprime_appGrowthRate.csv', sep = ""))
# plot the approximate growth rate for Gprime and the species-specific Cmax is shown with a vertical black line
pdf(paste('results/5_approximate_growth_rate/a_Cmax/', 'Gprime_appGrowthRate.pdf', sep = ""), width = 16, height = 4)
ggplot(Approx_growth_rate_Gprime, aes (x = Concentration, y = approx_growth_rate, color = Barcode)) +
  geom_point()+
  geom_line()+
  scale_x_continuous(trans = 'log2', minor_breaks = unique(Approx_growth_rate_Gprime$Concentration))+
  facet_wrap(~Species)+
  theme(panel.grid.minor.x =element_line(colour="grey", linetype = 'dashed'),
  panel.grid.major.x = element_line(colour="grey", linetype = 'dashed'))+
  geom_vline(aes(xintercept = Cmax), linetype = 1)+
  geom_point(aes(x = Cmax, y = Gprime_maxC_growthrate), color = 'black', size = 2)
dev.off()

## step b: calculate the approximate growth period (Tc) for each concentration  -------------------------------------------------------
Frequency_Initial_Final_maxC = left_join(Frequency_Initial_Final, Approx_growth_rate_Gprime_Cmax_median)
Tc_calculation <- Frequency_Initial_Final_maxC %>%
  filter(Shorthand == 'pAEGM') %>%
  mutate(Tc_per_barcode = (1/Gprime_maxC_growthrate)*(log(Density/InitialDensity)+log(BarcodeFrequency/InitialBarcodeFrequency)))
Tc_calculation_median <- Tc_calculation %>%
  group_by(Species, Concentration, Cmax) %>%
  summarise(Tc = median(Tc_per_barcode))
Tc_out <- left_join(Tc_calculation, Tc_calculation_median)
Tc_out <- Tc_out %>%
  select(Species, Genotype, Shorthand, Barcode, Concentration, Tc_per_barcode, Tc)
write_csv(Tc_out, paste('results/5_approximate_growth_rate/b_Tc/', 'Tc.csv', sep = ""))
pdf(paste('results/5_approximate_growth_rate/b_Tc/', 'Tc.pdf', sep = ""), width = 16, height = 4)
ggplot() +
  geom_line(data = Tc_calculation %>% filter(Cmax >= Concentration), aes (x = Concentration, y = Tc_per_barcode, group = Barcode, color = Barcode))+
  geom_point(data = Tc_calculation %>% filter(Cmax >= Concentration), aes (x = Concentration, y = Tc_per_barcode, group = Barcode, color = Barcode))+
  geom_line(data = Tc_calculation_median %>% filter(Cmax >= Concentration), aes(x = Concentration, y = Tc), color = 'black')+
  scale_x_continuous(trans = 'log2', minor_breaks = unique(Tc_calculation$Concentration))+
  theme(panel.grid.minor.x =element_line(colour="grey", linetype = 'dashed'),
        panel.grid.major.x = element_line(colour="grey", linetype = 'dashed'))+
  facet_grid(~Species)+
  geom_hline(yintercept = 24)
dev.off()

## step c: calculate the approximate growth rate for each barcode using concentration specific Tc  -------------------------------------------------------
Frequency_Initial_Final_Tc = left_join(Frequency_Initial_Final, Tc_calculation_median)

AppGrowthRate_genotype <- Frequency_Initial_Final_Tc %>%
  mutate(BarcodeFrequencyNA = ifelse(is.na(BarcodeFrequency), T, F),
         BarcodeCountsRecalculated = ifelse(is.na(BarcodeCounts), 1, BarcodeCounts),
         BarcodeFrequencyRecalculated = BarcodeCountsRecalculated/Reads) %>%
  mutate(AppGrowthRate = (1/Tc)*(log(Density/InitialDensity)+log(BarcodeFrequencyRecalculated/InitialBarcodeFrequency))) 
write_csv(Tc_out, paste('results/5_approximate_growth_rate/c_AppGrowthRate/', 'AppGrowthRate.csv', sep = ""))

pdf(paste('results/5_approximate_growth_rate/c_AppGrowthRate/', 'AppGrowthRate.pdf', sep = ""), width = 16, height = 6)
ggplot(data = AppGrowthRate_genotype %>% filter(Cmax >= Concentration), aes(x = Concentration, y = AppGrowthRate, group = Barcode, color = Shorthand))+
  geom_line()+
  geom_point(aes(fill = BarcodeFrequencyNA), shape = 22)+
  scale_x_continuous(trans = 'log2', minor_breaks = unique(AppGrowthRate_genotype$Concentration))+
  facet_grid(~Species)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.minor.x =element_line(colour="grey", linetype = 'dashed'),
        panel.grid.major.x = element_line(colour="grey", linetype = 'dashed'))+
  scale_fill_manual(values = c('white', 'black'))
dev.off()

pdf(paste('results/5_approximate_growth_rate/c_AppGrowthRate/', 'AppGrowthRate_byGenotype.pdf', sep = ""), width = 75, height = 6)
ggplot(data = AppGrowthRate_genotype %>% filter(Cmax >= Concentration), aes(x = Concentration, y = AppGrowthRate, group = Barcode, color = Shorthand))+
  geom_line()+
  geom_point(aes(fill = BarcodeFrequencyNA), shape = 22)+
  scale_x_continuous(trans = 'log2', minor_breaks = unique(AppGrowthRate_genotype$Concentration))+
  facet_grid(Species~Shorthand)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.major.x = element_line(colour="grey", linetype = 'dashed'))+
  scale_fill_manual(values = c('white', 'black'))
dev.off()

## step d: determine each allele's outlier using pair-wise sum squares--------
BarcodeReplicateLabel <- read.csv(filename_genotype_map)
BarcodeReplicateLabel <- BarcodeReplicateLabel %>%
  select(Shorthand, Barcode, ReplicateLabel)
AppGrowthRate_RepLabel_all <- left_join(AppGrowthRate_genotype, BarcodeReplicateLabel)
AppGrowthRate_RepLabel <- AppGrowthRate_RepLabel_all %>%
  filter(Cmax > Concentration)
AppGrowthRate_RepLabel_filter <- AppGrowthRate_RepLabel %>%
  select(Species, Genotype, Shorthand, Concentration, AppGrowthRate, ReplicateLabel) %>%
  pivot_wider(names_from = ReplicateLabel, values_from = AppGrowthRate) %>%
  mutate(A_B = (A - B)^2,
         A_C = (A - C)^2,
         B_C = (B - C)^2) %>%
  group_by(Species, Shorthand) %>%
  summarise(A_B = sum(A_B),
            A_C = sum(A_C),
            B_C = sum(B_C)) %>%
  mutate(A = A_B + A_C,
         B = A_B + B_C,
         C = A_C + B_C) %>%
  select(Species, Shorthand, A, B, C) %>%
  pivot_longer(cols = c('A', 'B', 'C'), names_to = 'ReplicateLabel', values_to = 'Pairwise_sum_squares') %>%
  group_by(Species, Shorthand) %>%
  mutate(max_pairwise = max(Pairwise_sum_squares), 
         Outlier = ifelse(Pairwise_sum_squares == max_pairwise, T, F)) %>%
  select(-max_pairwise)
AppGrowthRate_Outliers_plot <- left_join(AppGrowthRate_RepLabel, AppGrowthRate_RepLabel_filter)
AppGrowthRate_Outliers <- left_join(AppGrowthRate_RepLabel_all, AppGrowthRate_RepLabel_filter)
write_csv(AppGrowthRate_Outliers, paste('results/5_approximate_growth_rate/d_Outlier/', 'Outliers.csv', sep = ""))

pdf(paste('results/5_approximate_growth_rate/d_Outlier/', 'AppGrowthRate_Outliers.pdf', sep = ""), width = 16, height = 6)
ggplot(data = AppGrowthRate_Outliers_plot, aes(x = Concentration, y = AppGrowthRate, group = Barcode, color = Outlier))+
  geom_line()+
  geom_point(aes(fill = BarcodeFrequencyNA), shape = 22)+
  scale_x_continuous(trans = 'log2', minor_breaks = unique(AppGrowthRate_Outliers_plot$Concentration))+
  facet_grid(~Species)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.minor.x =element_line(colour="grey", linetype = 'dashed'),
        panel.grid.major.x = element_line(colour="grey", linetype = 'dashed'))+
  scale_fill_manual(values = c('white', 'black'))
dev.off()

pdf(paste('results/5_approximate_growth_rate/d_Outlier/', 'AppGrowthRate_Outliers_perGenotype.pdf', sep = ""), width = 75, height = 6)
ggplot(data = AppGrowthRate_Outliers_plot, aes(x = Concentration, y = AppGrowthRate, group = Barcode, color = Outlier))+
  geom_line()+
  geom_point(aes(fill = BarcodeFrequencyNA), shape = 22)+
  scale_x_continuous(trans = 'log2', minor_breaks = unique(AppGrowthRate_Outliers_plot$Concentration))+
  facet_grid(Species~Shorthand)+
  geom_hline(yintercept = 0)+
  theme(panel.grid.minor.x =element_line(colour="grey", linetype = 'dashed'),
        panel.grid.major.x = element_line(colour="grey", linetype = 'dashed'))+
  scale_fill_manual(values = c('white', 'black'))
dev.off()