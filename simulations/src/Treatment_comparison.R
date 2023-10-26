# Packages ------
library(tidyverse)
library(ggpubr)
library(cowplot)
source("src/Random_Landscape_Functions.R")
source("src/Crowdsourcing_plotting_settings.R")
theme_set(theme_cowplot())
# Load ----
df <- read.csv("results/random/4_treatment_comparison/combined_data.csv")

# Overview plots----
bin_count_plot <- ggplot(data = df %>%
                           distinct(Partition_ID, bin, .keep_all = TRUE), 
                         aes(x = bin))+
                  geom_bar()
ggsave("results/random/4_treatment_comparison/bin_count.pdf", bin_count_plot)

# Assign significance--------
# Perform Wilcoxon test for each Partition_ID and calculate adjusted p-value, also calculate mean resistance for each treatment
mapping <- c('Ec_N' = 'Focal only', 'Ec_HGT' = 'Focal with HGT to transient')
result <- df %>%
  mutate(Treatment_ID = recode(Treatment_ID, !!!mapping)) %>%
  mutate(Treatment_ID = factor(Treatment_ID, levels = c('Focal only', 'Focal with HGT to transient')))%>%
  group_by(Partition_ID, bin) %>%
  nest() %>%
  mutate(p_value = map_dbl(data, ~ perform_wilcoxon(.x$Resistance[.x$Treatment_ID == 'Focal only'],
                                                    .x$Resistance[.x$Treatment_ID == 'Focal with HGT to transient'])),
         Ec_N_mean = map_dbl(data, ~ calculate_mean_resistance(.x$Resistance[.x$Treatment_ID == 'Focal only'])),
         Ec_HGT_mean = map_dbl(data, ~ calculate_mean_resistance(.x$Resistance[.x$Treatment_ID == 'Focal with HGT to transient']))) %>%
  ungroup() %>%
  group_by(bin) %>%
  mutate(num_bin_comparisons = n()) %>%
  ungroup() %>%
  mutate(corrected_alpha = 0.05 / num_bin_comparisons) %>% #Bonferonni correction based on the number of comparisons within a bin
  mutate(significant = p_value < corrected_alpha) %>%
  mutate(difference = Ec_N_mean - Ec_HGT_mean) %>%
  mutate(outcome = case_when(significant & difference < 0 ~ "outsourcing",
                             significant & difference > 0 ~ "insourcing",
                             TRUE ~ "crowdsourcing")) %>%
  select(-corrected_alpha)


# Save the result to a new CSV file
output_file <- "results/random/4_treatment_comparison/evo_p_values_bincorrected.csv"
write_csv(result, output_file)

# Load alignment data ##########
alignment_df <- read.csv(file = "results/random/1_random_seeds/seed_alignment_scores.csv")
alignment_df <- alignment_df %>%
  select(-X) %>%
  rename(Partition_ID = seeds)
result <- left_join(result, alignment_df)

bin_verse_outcome <- result %>%
  group_by(bin, outcome) %>%
  reframe(n = n(),
          freq_outcome = n/50)

bin_verse_outcome <- result %>%
  group_by(bin, outcome) %>%
  summarise(n = n(), .groups = "drop") %>%  # Calculate count, drop grouping
  complete(bin, outcome, fill = list(n = 0)) %>%  # Fill missing combinations with n = 0
  group_by(bin) %>%  # Regroup by 'bin' to calculate frequencies
  mutate(freq_outcome = n / sum(n))  # Calculate frequencies

plot_output_file <- "results/random/4_treatment_comparison/misalignment_versus_outcome_frequency.csv"
write_csv(bin_verse_outcome, plot_output_file)

crowdsourcing_color <- rgb(255, 128, 128, maxColorValue = 255)
insourcing_color <- rgb(180, 130, 180, maxColorValue = 255)
outsourcing_color <- rgb(255, 200, 140, maxColorValue = 255)

frequency <- ggplot(data = bin_verse_outcome, aes(x = bin, y = freq_outcome, color = outcome)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x)+
  geom_point(data = bin_verse_outcome %>% filter(outcome == 'crowdsourcing'), shape = 23, size = 4, fill = crowdsourcing_color)+
  geom_point(data = bin_verse_outcome %>% filter(outcome == 'outsourcing'), shape = 24, size = 4, fill = outsourcing_color)+
  geom_point(data = bin_verse_outcome %>% filter(outcome == 'outsourcing'), shape = 8, size = 1.5, color = 'black')+
  geom_point(data = bin_verse_outcome %>% filter(outcome == 'insourcing'), shape = 25, size = 4, fill = insourcing_color)+
  geom_point(data = bin_verse_outcome %>% filter(outcome == 'insourcing'), shape = 8, size = 1.5, color = 'black')+
  stat_cor(method = "pearson", label.x = 8) +
  ylab('frequency of evolutionary outcome') +
  xlab('misalignment score')+
  scale_color_manual(values = c("crowdsourcing" = crowdsourcing_color,
                                "outsourcing" = outsourcing_color,
                                "insourcing" = insourcing_color))+
  theme(legend.position = "none")

ggsave(filename = "results/random/4_treatment_comparison/frequency_outcome_versus_landscape_misalignment.pdf", 
       width = 6, height = 6)

# frequency <- ggplot(data = bin_verse_outcome, aes(x = bin, y = freq_outcome, color = outcome)) +
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'crowdsourcing'), shape = 5, size = 4)+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'outsourcing'), shape = 2, size = 4)+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'outsourcing'), shape = 8, size = 1.5)+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'insourcing'), shape = 6, size = 4)+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'insourcing'), shape = 8, size = 1.5)+
#   geom_smooth(method = "lm", se = FALSE, formula = y ~ x)+
#   stat_cor(method = "pearson") +
#   ylab('frequency of evolutionary outcome') +
#   xlab('misalignment score')+
#   scale_color_manual(values = c("crowdsourcing" = crowdsourcing_color,
#                                 "outsourcing" = 'light grey',
#                                 "insourcing" = 'dark grey')) 
# frequency
# frequency1 <- ggplot(data = bin_verse_outcome, aes(x = bin, y = freq_outcome, color = outcome)) +
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'crowdsourcing'), shape = 23, size = 4, fill = crowdsourcing_color)+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'outsourcing'), shape = 24, size = 4, fill = outsourcing_color)+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'outsourcing'), shape = 8, size = 1.5, color = 'black')+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'insourcing'), shape = 25, size = 4, fill = insourcing_color)+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'insourcing'), shape = 8, size = 1.5, color = 'black')+
#   geom_smooth(method = "lm", se = FALSE, formula = y ~ x)+
#   stat_cor(method = "pearson") +
#   ylab('frequency of evolutionary outcome') +
#   xlab('misalignment score')+
#   scale_color_manual(values = c("crowdsourcing" = crowdsourcing_color,
#                                 "outsourcing" = outsourcing_color,
#                                 "insourcing" = insourcing_color)) 
# frequency1
# outsourcing_color <- rgb(160, 100, 160, maxColorValue = 255)
# frequency2 <- ggplot(data = bin_verse_outcome, aes(x = bin, y = freq_outcome, color = outcome)) +
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'crowdsourcing'), shape = 23, size = 4, fill = insourcing_color)+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'outsourcing'), shape = 24, size = 4, fill = outsourcing_color)+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'outsourcing'), shape = 8, size = 1.5, color = 'black')+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'insourcing'), shape = 25, size = 4, fill = crowdsourcing_color)+
#   geom_point(data = bin_verse_outcome %>% filter(outcome == 'insourcing'), shape = 8, size = 1.5, color = 'black')+
#   geom_smooth(method = "lm", se = FALSE, formula = y ~ x)+
#   stat_cor(method = "pearson") +
#   ylab('frequency of evolutionary outcome') +
#   xlab('misalignment score')+
#   scale_color_manual(values = c("crowdsourcing" = insourcing_color,
#                                 "outsourcing" = outsourcing_color,
#                                 "insourcing" = crowdsourcing_color))+
#   theme(legend.position = "none")
# frequency2


# Partition sub-plots #######
n <- 32  # Number of colors you want to generate
default_scale <- scale_colour_hue()
palette_fun <- default_scale$palette
colors <- palette_fun(n)
fill_palette <- c(colors, '#C0C0C0')
fill_mapping <- setNames(fill_palette, 1:33)

plot_list <- result %>%
  group_split(Partition_ID) %>%
  map(create_partition_plots)

#create_partition_plots <- function(df) {
  # df <- result %>%
  #   filter(Partition_ID == 6622)
  # Figure_width = 8 # in width for one column figure in PNAS
  # Figure_height = 8
# 
# landscape_plot <- function(Pid, focal_host) {
#     landscapes <- random_landscape(seed = Pid, sd_value = 0.125)
#     master.landscape.df <- landscapes %>%
#       mutate(num_of_mutations = g + A + E + G + M) 
#     steps_file <- read.csv('data/random/mutational_steps.csv', header = T)
#     steps <- steps_file %>%
#       mutate(step = paste(focal, mutant, sep = "->"))
#     focal_node <- steps %>%
#       select(-mutant) %>%
#       rename(Genotype = focal)
#     mutant_step <- steps  %>%
#       select(-focal) %>%
#       rename(Genotype = mutant)
#     focal <- left_join(master.landscape.df, focal_node, multiple = "all")
#     focal$identity <- 'focal'
#     mutant <- left_join(master.landscape.df, mutant_step, multiple = "all")
#     mutant$identity <- 'mutant'
#     full_data <- full_join(focal, mutant) 
#     full <- full_data %>%
#       separate(col = step, sep = '->', into = c('focal', 'mutant'), remove = F) %>%
#       filter(!(is.na(focal))) %>%
#       group_by(Species, step) %>%
#       mutate(focal_MIC = mic[Genotype == focal], 
#              mutant_MIC = mic[Genotype == mutant],
#              step_slope = (mutant_MIC) - (focal_MIC), 
#              step_effect = case_when((mutant_MIC) > (focal_MIC) ~ 'beneficial',
#                                      (mutant_MIC) < (focal_MIC) ~ 'deleterious',
#                                      TRUE ~ 'neutral')) %>%
#       group_by(step) %>%
#       mutate(host_epistasis = ifelse(length(unique(step_effect)) == 1, 'noGxH', 'GxH')) %>%
#       group_by(Genotype, Species) %>%
#       mutate(peak = all((identity == 'focal' & step_effect == 'deleterious')|
#                           (identity == 'mutant' & step_effect == 'beneficial')),
#              peak_color = ifelse(peak == T, id, 33)) %>%
#       ungroup() %>%
#       group_by(Genotype) %>%
#       mutate(not_shared_peak = any(peak == T) & any(peak == F),
#              peak_color = ifelse(any(peak == T), id, 33))
#   
#     color_mapping <- setNames(c('black','white', 'grey','black'), c(T, F, 'noGxH', 'GxH'))
#     
#     perturbed_landscape_plot <- ggplot(data = full %>% filter(Species == focal_host), aes(x = num_of_mutations, y = mic, group = step)) +
#       geom_line(aes(linetype = step_effect, color = host_epistasis), size = 0.6)+
#       geom_point(aes(fill = as.factor(peak_color), color = peak), shape = 21, size = 3, stroke = 0.5) +
#       scale_fill_manual(values = fill_mapping)+
#       scale_color_manual(values = color_mapping)+
#       scale_linetype_manual(values = c('solid','longdash','dotted')) +
#       ylab('Fitness')+
#       theme(legend.position = "none")+
#       xlab('Number of mutations') 
#     return(perturbed_landscape_plot)
# }

# mutation_effects_plot <- function(Pid) {
#   landscapes <- random_landscape(seed = Pid, sd_value = sd_value)
#   mutant_steps <- step_effects(landscapes)
#   perturbed_step_effects <- mutant_steps %>%
#     select(Species, step, step_slope) %>%
#     pivot_wider(names_from = Species, values_from = step_slope) %>%
#     mutate(GxH = (Ec > 0 & Kp < 0) | (Ec < 0 & Kp > 0))
#   mutant_effects <- ggplot(perturbed_step_effects, aes(x = Ec, y = Kp)) +
#     geom_point(aes(color = GxH))+
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
#     geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey")+
#     theme(legend.position = "none")+
#     scale_color_manual(values = c('grey', 'black'))+
#     scale_y_continuous(limits = c(-0.5, 1))+
#     scale_x_continuous(limits = c(-0.5, 1))+
#     ylab("mutation effect in transient")+
#     xlab("mutation effect in focal")
#   return(mutant_effects)
# }


# 
# evo_outcome_plot <- function(df) {
#   y_axis_label = expression(paste("fitness"))
#   margins <- c(0.1,0.1,0.1,0.1)
#   point_size = 2
#   landscapes <- random_landscape(seed = df$Partition_ID, sd_value = 0.125)
#   master.landscape.df <- landscapes %>%
#     mutate(num_of_mutations = g + A + E + G + M) 
#   steps_file <- read.csv('data/random/mutational_steps.csv', header = T)
#   steps <- steps_file %>%
#     mutate(step = paste(focal, mutant, sep = "->"))
#   focal_node <- steps %>%
#     select(-mutant) %>%
#     rename(Genotype = focal)
#   mutant_step <- steps  %>%
#     select(-focal) %>%
#     rename(Genotype = mutant)
#   focal <- left_join(master.landscape.df, focal_node, multiple = "all")
#   focal$identity <- 'focal'
#   mutant <- left_join(master.landscape.df, mutant_step, multiple = "all")
#   mutant$identity <- 'mutant'
#   full_data <- full_join(focal, mutant) 
#   full <- full_data %>%
#     separate(col = step, sep = '->', into = c('focal', 'mutant'), remove = F) %>%
#     filter(!(is.na(focal))) %>%
#     group_by(Species, step) %>%
#     mutate(focal_MIC = mic[Genotype == focal], 
#            mutant_MIC = mic[Genotype == mutant],
#            step_slope = (mutant_MIC) - (focal_MIC), 
#            step_effect = case_when((mutant_MIC) > (focal_MIC) ~ 'beneficial',
#                                    (mutant_MIC) < (focal_MIC) ~ 'deleterious',
#                                    TRUE ~ 'neutral')) %>%
#     group_by(step) %>%
#     mutate(host_epistasis = ifelse(length(unique(step_effect)) == 1, 'noGxH', 'GxH')) %>%
#     group_by(Genotype, Species) %>%
#     mutate(peak = all((identity == 'focal' & step_effect == 'deleterious')|
#                         (identity == 'mutant' & step_effect == 'beneficial')),
#            peak_color = ifelse(peak == T, id, 33)) %>%
#     ungroup() %>%
#     group_by(Genotype) %>%
#     mutate(not_shared_peak = any(peak == T) & any(peak == F),
#            peak_color = ifelse(any(peak == T), id, 33))
#   my_comparisons <- list(c('Focal only', 'Focal with HGT to transient'))
#   peak_info <- full %>%
#     ungroup() %>%
#     select(id, peak_color) %>%
#     rename(Genotype = id) %>%
#     distinct()
#   df$data[[1]] <- left_join(df$data[[1]], peak_info)
#     
#   # wrapped_label <- str_wrap(my_comparisons[[1]][2], width = 15)  # Adjust the width as needed
#   plot <- ggplot(data = df$data[[1]], aes(x = Treatment_ID, y = Resistance)) +
#     geom_violin(width = 1)+
#     geom_jitter(aes(color = as.factor(peak_color)), width = 0.15, size = point_size, alpha = 0.1)+
#     stat_summary(fun = "mean", geom = "point", color = "black", size = 3)+
#     #stat_summary(aes(color = Treatment_ID), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size, color = 'dark grey', width = 0.15)+
#     geom_hline(yintercept = df$Ec_N_mean, linetype = 'dashed', color = 'black')+
#     geom_bracket(xmin = 'Focal only', xmax = 'Focal with HGT to transient', y = max(df$data[[1]]$Resistance, na.rm = T) * 1.1, label = df$outcome[1])+
#     theme(legend.position = "none") +
#     ylab(y_axis_label)+
#     xlab(" ")+
#     scale_color_manual(values = fill_mapping)+
#     scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, by = 0.25))+
#     scale_x_discrete(breaks = c('Focal only', 'Focal with HGT to transient'),
#                      labels = c('F', 'F-T'))
#   plot
#   return(plot)
# }

# 
# landscape_focal <- landscape_plot(6622, 'Ec')
# landscape_transient <- landscape_plot(6622, 'Kp')  
# mutant_effects <- mutation_effects_plot(6622)
# evo_outcome <- evo_outcome_plot(df)
#   
# 



