# Packages ---------------------------------------------------------------
library('tidyverse')
library("ggpubr")
library(car)
library(cowplot)
theme_set(theme_cowplot())

# Load data ---------------------------------------------------------------
data_file <- read.csv("DHFR.csv")
data <- data_file %>%
  rename(Genotype = Shorthand)%>%
  mutate(num_of_mutations = P + A + L) %>%
  group_by(Species) %>%
  mutate(MIC_WT = MIC_mean[Genotype == 'WT'],
         MIC_rel = MIC_mean - MIC_WT)

subset_genotypes <- data %>%
  select(Species, Genotype, P, A, L, num_of_mutations)


# landscape graph ------
steps_file <- read.csv('mutational_steps_DHFR.csv', header = T)
steps <- steps_file %>%
  mutate(step = paste(focal, mutant, sep = "->"))

focal_node <- steps %>%
  select(-mutant) %>%
  rename(Genotype = focal)
mutant_step <- steps  %>%
  select(-focal) %>%
  rename(Genotype = mutant)

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

# mutation effect ------
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
write.csv(step_output_data_rel, file = 'mutational_step_slopes_rel.csv')

step_data <- full %>%
  group_by(Species)%>%
  filter(!(is.na(step)),
         identity == 'focal')

step_output_data <- step_data %>%
  mutate(focal_num = nchar(focal),
         mutant_num = nchar(mutant),
         focal_num = ifelse(focal == 'WT', 0, focal_num)) %>%
  select(Species, step, step_slope, step_effect, focal_MIC, mutant_MIC, focal_num, mutant_num)
write.csv(step_output_data, file = 'mutational_step_slopes.csv')

