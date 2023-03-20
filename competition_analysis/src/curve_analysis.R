## Packages ----------------------------------------------------------------
library(argparse)
library(drc)
library(tidyverse)
library(ggnewscale)
library(broom)
library(modelr)

## Command line ------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("-g", "--appGrowthRate")
parser$add_argument("-b", "--genotypeBinary")
parser$add_argument("-s", "--mutationalSteps")
args <- parser$parse_args()
AppGrowthRate_filename <- args$appGrowthRate
genotype_binary_filename <- args$genotypeBinary
steps_filename <- args$mutationalSteps

## Load data ---------------------------------------------------------------
AppGrowthRate_data <- read.csv(AppGrowthRate_filename)
genotype_binary <- read_csv(genotype_binary_filename)
steps_file <- read_csv(steps_filename)
# AppGrowthRate_data <- read.csv("results/5_approximate_growth_rate/d_Outlier/Outliers.csv")
# genotype_binary <- read_csv("data/genotype_coding.csv")
# steps_file <- read_csv(file = "data/mutational_steps.csv")

## Functions for 4 parameter model fitting and plotting  -------------------------------------------------------
drm.func <- function(x) {
  drm(AppGrowthRate ~ Concentration,
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
      data = x)}
predict.fun <- function(x) {
  add_predictions(data.frame(Concentration = c(seq(0.0156,1, 0.001),seq(1.001,515, 0.01))), x)}
predict.fun.untruncated <- function(x) {
  add_predictions(data.frame(Concentration = c(seq(0.0156,1, 0.001),seq(1.001, 100, 0.01),seq(100.001, 12000, 0.1))), x)}
coefs.fun <- function(x) {coef(x) %>% tidy}

## step a: determine the species-specific lower asymptote  -------------------------------------------------------
AppGrowthRate_AllData <- AppGrowthRate_data %>%
  select(Species, Genotype, Shorthand, Barcode, Concentration, Cmax, AppGrowthRate, Outlier)
Outlier_info <- AppGrowthRate_AllData %>%
  select(Species, Shorthand, Barcode, Outlier) %>%
  group_by(Species, Shorthand, Barcode, Outlier) %>%
  summarise()
AppGrowthRate_truncated <- AppGrowthRate_AllData %>%
  filter(Concentration <= Cmax) %>%
  filter(!(Shorthand %in% c('pAEGM'))) %>%
  select(Species, Shorthand, Barcode, Outlier, Concentration, AppGrowthRate) %>%
  mutate(Concentration = ifelse(Concentration == 0, 0.01575, Concentration))
df_truncated_4param <- AppGrowthRate_truncated %>%
  group_by(Species, Shorthand, Barcode) %>%
  nest() %>%
  mutate(drmod = map(data, drm.func),
         coefs = map(drmod, coefs.fun))
df_truncated_coefs_nested <- df_truncated_4param %>%
  dplyr::select(-drmod, -data)
df_truncated_coefs <- df_truncated_coefs_nested %>% unnest(coefs) %>% spread(names, x)
LowerConcentrationCutoff <- tibble(Species = c('Ec','Kp','Se'), LowerCutoff = c(2,1,2))
coeffsum_4param_filter <- left_join(df_truncated_coefs, LowerConcentrationCutoff)
coeffsum_4param_lowerasymptote <- coeffsum_4param_filter %>%
  filter(`ED50:(Intercept)` < LowerCutoff) %>%
  group_by(Species) %>%
  summarise(LowerAsymptoteMean = mean(`Lower Limit:(Intercept)`))
AppGrowthRate_AllData_wLowerAsymptote <- left_join(AppGrowthRate_AllData, LowerConcentrationCutoff)
AppGrowthRate_AllData_wLowerAsymptote <- left_join(AppGrowthRate_AllData_wLowerAsymptote, coeffsum_4param_lowerasymptote)
write_csv(AppGrowthRate_AllData_wLowerAsymptote, paste('results/6_dose_curves/a_lowerAsymptote/', 'lowerAsymptote.csv', sep = ""))

pdf(paste('results/6_dose_curves/a_lowerAsymptote/', 'lowerasymptote_modelfitLL4.pdf'), width = 20, height = 6)
ggplot() +
  geom_line(data = AppGrowthRate_AllData_wLowerAsymptote %>% filter(Concentration <= Cmax), aes(Concentration, AppGrowthRate, color = Shorthand, group = Barcode), size = 0.5) +
  geom_point(data = AppGrowthRate_AllData_wLowerAsymptote %>% filter(Concentration <= Cmax), aes(Concentration, AppGrowthRate, color = Shorthand, fill = Outlier), size = 1, shape = 22) +
  geom_vline(data = AppGrowthRate_AllData_wLowerAsymptote, aes(xintercept = LowerCutoff), color = 'black') +
  geom_hline(data = AppGrowthRate_AllData_wLowerAsymptote, aes(yintercept = LowerAsymptoteMean), color = 'black') +
  scale_x_continuous(trans = 'log2', limits = c(0.0156, 515)) +
  scale_y_continuous(limits = c(-0.1, 0.75)) +
  theme_bw() +
  facet_grid(~Species) +
  scale_fill_manual(values = c('white', 'black'))
dev.off()

## step b: curve fitting with species-specific fixed lower asymptote ---------------
untrunctated_data_aboveCmax <- AppGrowthRate_data %>%
  filter(Shorthand %in% c('pAEGM', 'pEGM')) %>%
  filter(Concentration > Cmax) %>%
  mutate(AppGrowthRate = (1/24)*(log(Density/InitialDensity)+log(BarcodeFrequency/InitialBarcodeFrequency))) %>%
  select(Species, Genotype, Shorthand, Barcode, Concentration, Cmax, AppGrowthRate, ReplicateLabel)
untruncated_data_belowCmax <- AppGrowthRate_data %>%
  filter(Concentration <= Cmax) %>%
  filter(Shorthand %in% c('pAEGM', 'pEGM')) %>%
  select(Species, Genotype, Shorthand, Barcode, Concentration, Cmax, AppGrowthRate, ReplicateLabel)
untruncated_all_data <- rbind(untruncated_data_belowCmax, untrunctated_data_aboveCmax)

ApproxGrowthRate_untruncated_outliers <- untruncated_all_data %>%
  dplyr::select(Species, Genotype, Shorthand, Concentration, AppGrowthRate, ReplicateLabel) %>%
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
  dplyr::select(Species, Shorthand, A, B, C) %>%
  pivot_longer(cols = c('A', 'B', 'C'), names_to = 'ReplicateLabel', values_to = 'Pairwise_sum_squares') %>%
  group_by(Species, Shorthand) %>%
  mutate(max_pairwise = max(Pairwise_sum_squares),
         Outlier = ifelse(Pairwise_sum_squares == max_pairwise, T, F)) %>%
  dplyr::select(-max_pairwise)
label_info <- untruncated_all_data %>%
  group_by(Species, Shorthand, Barcode, ReplicateLabel) %>%
  summarise()
untruncated_outliers <- left_join(ApproxGrowthRate_untruncated_outliers, label_info)
untruncated_outliers <- untruncated_outliers %>%
  select(Species, Shorthand, Barcode, Outlier)
untruncated_data <- left_join(untruncated_all_data, untruncated_outliers)
untruncated_data <- untruncated_data  %>%
  select(Species, Genotype, Shorthand, Barcode, Concentration, Cmax, AppGrowthRate, Outlier)

# Ec 
drm.func.Ec <- function(x) {
  drm(AppGrowthRate ~ Concentration,
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                 fixed = c(NA, 
                           coeffsum_4param_lowerasymptote$LowerAsymptoteMean[coeffsum_4param_lowerasymptote$Species == 'Ec'], #0.1135523
                           NA, NA)), data = x)}
Ec_truncated <- AppGrowthRate_truncated %>%
  filter(!(Shorthand %in% c('null', 'pEGM')) & Species == 'Ec') %>%
  mutate(Concentration = ifelse(Concentration == 0, 0.01575, Concentration))
Ec_trunc_4param <- Ec_truncated %>%
  group_by(Species, Shorthand, Barcode) %>%
  nest() %>%
  mutate(drmod = map(data, drm.func.Ec),
         pred = map(drmod, predict.fun),
         coefs = map(drmod, coefs.fun))
Ec_trunc_4param_outliers <- left_join(Ec_trunc_4param, Outlier_info)
Ec_truncated$truncated <- T
Ec_trunc_4param_outliers$truncated <- T

Ec_untruncated <- untruncated_data %>%
  filter(Species == 'Ec') %>%
  select(Species, Shorthand, Barcode, Outlier, Concentration, AppGrowthRate)
Ec_untrunc_4param <- Ec_untruncated %>%
  group_by(Species, Shorthand, Barcode) %>%
  nest() %>%
  mutate(drmod = map(data, drm.func.Ec),
         pred = map(drmod, predict.fun.untruncated),
         coefs = map(drmod, coefs.fun))
Ec_untrunc_4param_outliers <- left_join(Ec_untrunc_4param, untruncated_outliers)
Ec_untruncated$truncated <- F
Ec_untrunc_4param_outliers$truncated <- F

# Kp 
drm.func.Kp <- function(x) {
  drm(AppGrowthRate ~ Concentration,
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                 fixed = c(NA, 
                           coeffsum_4param_lowerasymptote$LowerAsymptoteMean[coeffsum_4param_lowerasymptote$Species == 'Kp'], #0.05838588 
                           NA, NA)), data = x)}
Kp_truncated <- AppGrowthRate_truncated %>%
  filter(!(Shorthand %in% c('null', 'pEGM')) & Species == 'Kp') %>%
  mutate(Concentration = ifelse(Concentration == 0, 0.01575, Concentration))
Kp_trunc_4param <- Kp_truncated %>%
  group_by(Species, Shorthand, Barcode) %>%
  nest() %>%
  mutate(drmod = map(data, drm.func.Kp),
         pred = map(drmod, predict.fun),
         coefs = map(drmod, coefs.fun))
Kp_trunc_4param_outliers <- left_join(Kp_trunc_4param, Outlier_info)
Kp_truncated$truncated <- T
Kp_trunc_4param_outliers$truncated <- T

Kp_untruncated <- untruncated_data %>%
  filter(Species == 'Kp') %>%
  select(Species, Shorthand, Barcode, Outlier, Concentration, AppGrowthRate)
Kp_untrunc_4param <- Kp_untruncated %>%
  group_by(Species, Shorthand, Barcode) %>%
  nest() %>%
  mutate(drmod = map(data, drm.func.Kp),
         pred = map(drmod, predict.fun.untruncated),
         coefs = map(drmod, coefs.fun))
Kp_untrunc_4param_outliers <- left_join(Kp_untrunc_4param, untruncated_outliers)
Kp_untruncated$truncated <- F
Kp_untrunc_4param_outliers$truncated <- F

# Se 
drm.func.Se <- function(x) {
  drm(AppGrowthRate ~ Concentration,
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                 fixed = c(NA, 
                           coeffsum_4param_lowerasymptote$LowerAsymptoteMean[coeffsum_4param_lowerasymptote$Species == 'Se'], 
                           NA, NA)), data = x)}
Se_truncated <- AppGrowthRate_truncated %>%
  filter(!(Shorthand %in% c('null')) & Species == 'Se') %>%
  mutate(Concentration = ifelse(Concentration == 0, 0.01575, Concentration)) %>%
  mutate(Concentration = ifelse(Shorthand == 'pEGM' & Concentration == 0.01575, 0, Concentration))
Se_trunc_4param <- Se_truncated %>%
  group_by(Species, Shorthand, Barcode) %>%
  nest() %>%
  mutate(drmod = map(data, drm.func.Se),
         pred = map(drmod, predict.fun),
         coefs = map(drmod, coefs.fun))
Se_trunc_4param_outliers <- left_join(Se_trunc_4param, Outlier_info)
Se_truncated$truncated <- T
Se_trunc_4param_outliers$truncated <- T

Se_untruncated <- untruncated_data %>%
  filter(Species == 'Se', Shorthand == 'pAEGM') %>%
  select(Species, Shorthand, Barcode, Outlier, Concentration, AppGrowthRate)
Se_untrunc_4param <- Se_untruncated %>%
  group_by(Species, Shorthand, Barcode) %>%
  nest() %>%
  mutate(drmod = map(data, drm.func.Se),
         pred = map(drmod, predict.fun.untruncated),
         coefs = map(drmod, coefs.fun))
Se_untrunc_4param_outliers <- left_join(Se_untrunc_4param, untruncated_outliers)
Se_untruncated$truncated <- F
Se_untrunc_4param_outliers$truncated <- F

# merge model fitting data
all_data <- rbind(Ec_truncated, Ec_untruncated, Kp_truncated, Kp_untruncated, Se_truncated, Se_untruncated)
all_4param <- rbind(Ec_trunc_4param_outliers, Ec_untrunc_4param_outliers, Kp_trunc_4param_outliers, Kp_untrunc_4param_outliers, Se_trunc_4param_outliers, Se_untrunc_4param_outliers)
all_4param_coeffs <- all_4param %>%
  dplyr::select(-drmod, -data, -pred)
all_4param_coefs <- all_4param_coeffs %>% unnest(coefs) %>% spread(names, x)

write_csv(all_4param_coefs, paste('results/6_dose_curves/b_modelFit/', 'modelFit_coefs.csv', sep = ""))

pdf(paste('results/6_dose_curves/b_modelFit/', 'modelfit_bySpecies.pdf'), width = 25, height = 6)
ggplot() +
  geom_line(data = all_data, aes(Concentration, AppGrowthRate, color = Shorthand, group = Barcode, linetype = Outlier, alpha = Outlier), size = 0.5) +
  new_scale_color() + 
  geom_point(data = all_data, aes(Concentration, AppGrowthRate, color = Outlier, fill = Shorthand, alpha = Outlier), size = 2, shape = 22) +
  scale_color_manual(values = c('white', 'black')) +
  new_scale_color() + 
  geom_line(aes(Concentration, pred, color = Shorthand, group = Barcode, alpha = Outlier), data = all_4param %>% unnest(pred), size = 1.5) +
  geom_vline(aes(xintercept = x, color = Shorthand, group = Barcode, alpha = Outlier), linetype = 5, data = all_4param %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)")) +
  scale_x_continuous(trans = 'log2', limits = c(0.0156, 515))+
  theme_bw()+
  scale_alpha_discrete(range=c(1, 0.3))+
  facet_grid(~Species)
dev.off()

pdf(paste('results/6_dose_curves/b_modelFit/', 'modelfit_byGenotype.pdf'), width = 70, height = 6)
ggplot() +
  geom_line(data = all_data, aes(Concentration, AppGrowthRate, color = Shorthand, group = Barcode, linetype = Outlier, alpha = Outlier), size = 0.5) +
  new_scale_color() + 
  geom_point(data = all_data, aes(Concentration, AppGrowthRate, color = Outlier, fill = Shorthand, alpha = Outlier), size = 2, shape = 22) +
  scale_color_manual(values = c('white', 'black')) +
  new_scale_color() + 
  geom_line(aes(Concentration, pred, color = Shorthand, group = Barcode, alpha = Outlier), data = all_4param %>% unnest(pred), size = 1) +
  geom_vline(aes(xintercept = x, color = Shorthand, group = Barcode, alpha = Outlier), linetype = 5, data = all_4param %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)")) +
  scale_x_continuous(trans = 'log2', limits = c(0.0156, 515))+
  theme_bw()+
  scale_alpha_discrete(range=c(1, 0.3))+
  facet_grid(Species~Shorthand)
dev.off()

## step c: calculate genotype-level MIC with sqrt 2 errorbar -------------------------------
genotype_binary_filter <- genotype_binary %>%
  dplyr::select(-Genotype)
coeff_sum_LogRoot2 <- all_4param_coefs %>%
  mutate(MIC_LogRoot2 = log(`ED50:(Intercept)`, base = sqrt(2))) %>%
  filter(Outlier == FALSE)
MIC_summary <- coeff_sum_LogRoot2 %>%
  group_by(Species, Shorthand) %>%
  summarise(N = n(),
            MIC_mean = mean(MIC_LogRoot2),
            MIC_sd = sd(MIC_LogRoot2),
            MIC_se = MIC_sd/sqrt(N)) %>%
  select(-MIC_sd, -N) 
MIC_output <- left_join(MIC_summary, genotype_binary_filter)
MIC_out <- MIC_output %>%
  rename(hosts = Species, fit.1 = MIC_mean) %>%
  mutate(Num_of_mutations = g + A + E + G + M) %>%
  select(hosts, g, A, E, M, G, fit.1, MIC_se, Num_of_mutations)
write_csv(MIC_out, 'results/6_dose_curves/c_resistance/resistance_levels.csv')

three_barcodes <- all_4param_coefs %>%
  ungroup()%>%
  select(Species, Shorthand, Outlier, `ED50:(Intercept)`) %>%
  group_by(Species, Shorthand, Outlier) %>%
  mutate(id = row_number(),
         id = ifelse(Outlier == T, 3, id), 
         id = recode(id,
                     '1' = 'fit.1',
                     '2' = 'fit.2',
                     '3' = 'outlier')) %>%
  ungroup() %>%
  select(-Outlier)
three_barcodes_wider <- three_barcodes %>%
  pivot_wider(names_from = id, values_from = `ED50:(Intercept)`) %>%
  select(Species, Shorthand, fit.1, fit.2, outlier)
MIC_merge <- MIC_output %>%
  rename(fit.log = MIC_mean) %>%
  select(-MIC_se) %>%
  mutate(fit_average = sqrt(2)^fit.log) %>%
  select(Species, Shorthand, fit_average, fit.log, g, A, E, G, M)
three_barcodes_out <- left_join(three_barcodes_wider, MIC_merge)
write_csv(three_barcodes_out, "results/6_dose_curves/c_resistance/MIC_simulation_file_3barcodes.csv")          

## step d: analyze the mutational steps in the resistance landscapes ---------------------
steps <- steps_file %>%
  mutate(step = paste(focal, mutant, sep = "->"))
focal_node <- steps %>%
  select(-mutant) %>%
  rename(Genotype = focal)
mutant_step <- steps  %>%
  select(-focal) %>%
  rename(Genotype = mutant)
data <- MIC_output %>%
  rename(Genotype = Shorthand)%>%
  mutate(num_of_mutations = g + A + E + G + M) %>%
  group_by(Species) %>%
  mutate(MIC_WT = MIC_mean[Genotype == 'WT'],
         MIC_rel = MIC_mean - MIC_WT)
focal <- left_join(data, focal_node)
focal$identity <- 'focal'
mutant <- left_join(data, mutant_step)
mutant$identity <- 'mutant'
full_data <- full_join(focal, mutant) 
full <- full_data %>%
  separate(col = step, sep = '->', into = c('focal', 'mutant'), remove = F) %>%
  filter(!(is.na(focal))) %>%
  group_by(Species, step) %>%
  mutate(focal_MIC = MIC_mean[Genotype == focal], 
         mutant_MIC = MIC_mean[Genotype == mutant],
         focal_se = MIC_se[Genotype == focal],
         mutant_se = MIC_se[Genotype == mutant],
         step_slope = (mutant_MIC) - (focal_MIC), 
         step_effect = case_when((mutant_MIC-mutant_se) > (focal_MIC+focal_se) ~ 'beneficial',
                                 (mutant_MIC+mutant_se) < (focal_MIC-focal_se) ~ 'deleterious',
                                 TRUE ~ 'neutral')) %>%
  group_by(step) %>%
  mutate(host_epistasis = ifelse(length(unique(step_effect)) == 1, F, T))
host_epistasis <- full %>%
  group_by(step) %>%
  filter(host_epistasis == T)
unique(host_epistasis$step)
full_rel <- full_data %>%
  separate(col = step, sep = '->', into = c('focal', 'mutant'), remove = F) %>%
  filter(!(is.na(focal))) %>%
  group_by(Species, step) %>%
  mutate(focal_MIC = MIC_rel[Genotype == focal], 
         mutant_MIC = MIC_rel[Genotype == mutant],
         focal_se = MIC_se[Genotype == focal],
         mutant_se = MIC_se[Genotype == mutant],
         step_slope = (mutant_MIC) - (focal_MIC), 
         step_effect = case_when((mutant_MIC-mutant_se) > (focal_MIC+focal_se) ~ 'beneficial',
                                 (mutant_MIC+mutant_se) < (focal_MIC-focal_se) ~ 'deleterious',
                                 TRUE ~ 'neutral')) %>%
  group_by(step) %>%
  mutate(host_epistasis = ifelse(length(unique(step_effect)) == 1, F, T))
step_data_rel <- full_rel %>%
  group_by(Species)%>%
  filter(!(is.na(step)),
         identity == 'focal')
step_output_data_rel <- step_data_rel %>%
  mutate(focal_num = nchar(focal),
         mutant_num = nchar(mutant),
         focal_num = ifelse(focal == 'WT', 0, focal_num)) %>%
  select(Species, step, step_slope, step_effect, focal_MIC, mutant_MIC, focal_num, mutant_num)
write.csv(step_output_data_rel, file = 'results/6_dose_curves/d_mutations/mutational_step_slopes_rel.csv')
step_data <- full %>%
  group_by(Species)%>%
  filter(!(is.na(step)),
         identity == 'focal')
step_output_data <- step_data %>%
  mutate(focal_num = nchar(focal),
         mutant_num = nchar(mutant),
         focal_num = ifelse(focal == 'WT', 0, focal_num)) %>%
  select(Species, step, step_slope, step_effect)
write.csv(step_output_data, file = 'results/6_dose_curves/d_mutations/mutational_step_slopes_only.csv')




