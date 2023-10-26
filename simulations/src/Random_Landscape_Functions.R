# General landsape functions ########
random_landscape <- function(seed, sd_value){
  Pid <- seed
  # Construct master landscapes #####
  set.seed(Pid)
  master.landscape.df <- read.csv("data/random/code_format.csv")
  mutations <- c('g', 'A', 'E', 'G', 'M')
  num_of_mutations <- length(mutations)
  num_of_hosts <- length(unique(master.landscape.df$Species))
  mutation_effects <- runif(5, min = 0, max = 0.2)
  for(mutation in mutations){
    replacement_value <- mutation_effects[mutations==mutation]
    master.landscape.df[paste0(mutation, "_effect")] <- 0
    master.landscape.df[paste0(mutation, "_effect")][master.landscape.df[mutation] == 1] <- replacement_value
  }
  master.landscape.df <- master.landscape.df %>%
    mutate(mic = g_effect + A_effect + E_effect + G_effect + M_effect)
  max_fitness <- max(master.landscape.df$mic)
  master.landscape.df$mic <- master.landscape.df$mic/max_fitness
  master.landscape.df$g_effect <- master.landscape.df$g_effect/max_fitness
  master.landscape.df$A_effect <- master.landscape.df$A_effect/max_fitness
  master.landscape.df$E_effect <- master.landscape.df$E_effect/max_fitness
  master.landscape.df$G_effect <- master.landscape.df$G_effect/max_fitness
  master.landscape.df$M_effect <- master.landscape.df$M_effect/max_fitness
  master.landscape.df$gen_id <- seq(nrow(master.landscape.df))
  # Perturb master landscapes ####
  master.landscape.df$master_mic <- master.landscape.df$mic
  num_gen_perturb <- sample(1:(((2^num_of_mutations)*num_of_hosts)-4), 1)
  potential_gens <- master.landscape.df$gen_id[!master.landscape.df$gen_id %in% c(1, 32, 33, 64)]
  genotypes_to_perturb <- sample(potential_gens, num_gen_perturb)
  master.landscape.df$gen_perturbed <- master.landscape.df$gen_id %in% genotypes_to_perturb
  for (gen in genotypes_to_perturb) {
    old_value <- master.landscape.df$mic[gen]
    valid_perturbed_value <- FALSE
    while (!valid_perturbed_value) {
      perturbed_value <- rnorm(1, mean = old_value, sd = sd_value)
      if (perturbed_value >= 0 && perturbed_value <= 1) {
        valid_perturbed_value <- TRUE
      }
    }
    master.landscape.df$mic[gen] <- perturbed_value
  }
  return(master.landscape.df)
}
calculate_alignment_score <- function(seed_random_landscape) {
  # Calculate mutational steps on perturbed landscapes ####
  master.landscape.df <- seed_random_landscape %>%
    mutate(num_of_mutations = g + A + E + G + M) 
  
  steps_file <- read.csv('data/random/mutational_steps.csv', header = T)
  steps <- steps_file %>%
    mutate(step = paste(focal, mutant, sep = "->"))
  
  focal_node <- steps %>%
    select(-mutant) %>%
    rename(Genotype = focal)
  mutant_step <- steps  %>%
    select(-focal) %>%
    rename(Genotype = mutant)
  
  focal <- left_join(master.landscape.df, focal_node, multiple = "all")
  focal$identity <- 'focal'
  mutant <- left_join(master.landscape.df, mutant_step, multiple = "all")
  mutant$identity <- 'mutant'
  
  full_data <- full_join(focal, mutant) 
  
  full <- full_data %>%
    separate(col = step, sep = '->', into = c('focal', 'mutant'), remove = F) %>%
    filter(!(is.na(focal))) %>%
    group_by(Species, step) %>%
    mutate(focal_MIC = mic[Genotype == focal], 
           mutant_MIC = mic[Genotype == mutant],
           step_slope = (mutant_MIC) - (focal_MIC), 
           step_effect = case_when((mutant_MIC) > (focal_MIC) ~ 'beneficial',
                                   (mutant_MIC) < (focal_MIC) ~ 'deleterious',
                                   TRUE ~ 'neutral')) %>%
    group_by(step) %>%
    mutate(host_epistasis = ifelse(length(unique(step_effect)) == 1, F, T)) 
  
  perturbed_step_effects <- full %>%
    filter(identity == 'focal') %>%
    select(Species, step, step_slope) %>%
    pivot_wider(names_from = Species, values_from = step_slope)
  
  # Calculate alignment metrics ####
  perturbed_step_effects$perpendicular_distance <- abs(perturbed_step_effects$Ec - perturbed_step_effects$Kp) / sqrt(2) #perpendicular distance from the identity line
  perturbed_step_effects$opposite_values <- ifelse(perturbed_step_effects$Ec * perturbed_step_effects$Kp < 0, TRUE, FALSE)
  perpendicular_distance_sum_all <- sum(perturbed_step_effects$perpendicular_distance)
  perpendicular_distance_sum_signeffects <- sum(perturbed_step_effects$perpendicular_distance[perturbed_step_effects$opposite_values])
  alignment_proportion_signeffects <- sum(perturbed_step_effects$opposite_values)/nrow(perturbed_step_effects)
  # Plots ------
  # mutant_effects <- ggplot(perturbed_step_effects, aes(x = Ec, y = Kp)) +
  #   geom_point()+
  #   geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  #   geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  #   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey")
  # 
  # perturbed_landscape_plot <- ggplot(data = full, aes(x = num_of_mutations, y = mic, group = step)) +
  #   geom_line(aes(linetype = step_effect, color = host_epistasis), size = 0.6)+
  #   geom_point(aes(fill = Genotype, color = gen_perturbed), shape = 21, size = 2, stroke = 1) +
  #   facet_grid(~Species)+
  #   scale_color_manual(values = c('grey','black'))+
  #   scale_linetype_manual(values = c('solid','longdash','dotted')) +
  #   ylab('Fitness')+
  #   xlab('Number of mutations') +
  #   theme(legend.position = "none")
  # 
  # combined_plot <- plot_grid(
  #   perturbed_landscape_plot + labs(title = paste("perturbed landscapes, N =", num_gen_perturb, ", alignment =", round(perpendicular_distance_sum_all, 4))),
  #   mutant_effects + labs(title = "mutant effects"),
  #   labels = "AUTO",  # Automatically label each plot with letters (a, b, c)
  #   ncol = 1)  # Arrange the plots in one column
  # 
  # # Save the combined plot to a PDF file
  # ggsave(paste("results/random/partition_landscapes/", Pid, "_", sd_value, "_landscapes", '.pdf', sep = ''), combined_plot, width = 7, height = 10)
  # Output ------
  return(c(perpendicular_distance_sum_all, perpendicular_distance_sum_signeffects, alignment_proportion_signeffects))
}
step_effects <- function(seed_random_landscape) {
  # Calculate mutational steps on perturbed landscapes ####
  master.landscape.df <- seed_random_landscape %>%
    mutate(num_of_mutations = g + A + E + G + M) 
  
  steps_file <- read.csv('data/random/mutational_steps.csv', header = T)
  steps <- steps_file %>%
    mutate(step = paste(focal, mutant, sep = "->"))
  
  focal_node <- steps %>%
    select(-mutant) %>%
    rename(Genotype = focal)
  mutant_step <- steps  %>%
    select(-focal) %>%
    rename(Genotype = mutant)
  
  focal <- left_join(master.landscape.df, focal_node, multiple = "all")
  focal$identity <- 'focal'
  mutant <- left_join(master.landscape.df, mutant_step, multiple = "all")
  mutant$identity <- 'mutant'
  
  full_data <- full_join(focal, mutant) 
  
  full <- full_data %>%
    separate(col = step, sep = '->', into = c('focal', 'mutant'), remove = F) %>%
    filter(!(is.na(focal))) %>%
    group_by(Species, step) %>%
    mutate(focal_MIC = mic[Genotype == focal], 
           mutant_MIC = mic[Genotype == mutant],
           step_slope = (mutant_MIC) - (focal_MIC), 
           step_effect = case_when((mutant_MIC) > (focal_MIC) ~ 'beneficial',
                                   (mutant_MIC) < (focal_MIC) ~ 'deleterious',
                                   TRUE ~ 'neutral')) %>%
    group_by(step) %>%
    mutate(host_epistasis = ifelse(length(unique(step_effect)) == 1, F, T)) 
  
  perturbed_step_effects <- full %>%
    filter(identity == 'focal') %>%
    select(Species, step, step_slope, step_effect)
  return(perturbed_step_effects)
}
find_bin <- function(landscape){
  range_start <- 0
  range_end <- 15
  bin_size <- 0.5
  bins <- seq(range_start, range_end, by = bin_size)
  bin <- bins[findInterval(calculate_alignment_score(landscape)[1], bins)]
  return(bin)
}

# Assigning significance (step 4) ####
perform_wilcoxon <- function(resistance1, resistance2) { # Wilcoxon rank-sum test, also known as the Mann-Whitney U test
  wilcox.test(resistance1, resistance2, correct = TRUE)$p.value}
calculate_mean_resistance <- function(resistance) { # Define a function to calculate mean resistance for a treatment
  mean(resistance, na.rm = TRUE)}

# Plotting functions (step 4) ####
MinMeanSEMMax <- function(x) {
  v <- c(min(x), 
         mean(x) - sd(x)/sqrt(length(x)), 
         mean(x), 
         mean(x) + sd(x)/sqrt(length(x)), 
         max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v}

evo_outcome_plot <- function(df) {
  y_axis_label = expression(paste("fitness"))
  margins <- c(0.1,0.1,0.1,0.1)
  point_size = 2
  landscapes <- random_landscape(seed = df$Partition_ID, sd_value = 0.125)
  master.landscape.df <- landscapes %>%
    mutate(num_of_mutations = g + A + E + G + M) 
  steps_file <- read.csv('data/random/mutational_steps.csv', header = T)
  steps <- steps_file %>%
    mutate(step = paste(focal, mutant, sep = "->"))
  focal_node <- steps %>%
    select(-mutant) %>%
    rename(Genotype = focal)
  mutant_step <- steps  %>%
    select(-focal) %>%
    rename(Genotype = mutant)
  focal <- left_join(master.landscape.df, focal_node, multiple = "all")
  focal$identity <- 'focal'
  mutant <- left_join(master.landscape.df, mutant_step, multiple = "all")
  mutant$identity <- 'mutant'
  full_data <- full_join(focal, mutant) 
  full <- full_data %>%
    separate(col = step, sep = '->', into = c('focal', 'mutant'), remove = F) %>%
    filter(!(is.na(focal))) %>%
    group_by(Species, step) %>%
    mutate(focal_MIC = mic[Genotype == focal], 
           mutant_MIC = mic[Genotype == mutant],
           step_slope = (mutant_MIC) - (focal_MIC), 
           step_effect = case_when((mutant_MIC) > (focal_MIC) ~ 'beneficial',
                                   (mutant_MIC) < (focal_MIC) ~ 'deleterious',
                                   TRUE ~ 'neutral')) %>%
    group_by(step) %>%
    mutate(host_epistasis = ifelse(length(unique(step_effect)) == 1, 'noGxH', 'GxH')) %>%
    group_by(Genotype, Species) %>%
    mutate(peak = all((identity == 'focal' & step_effect == 'deleterious')|
                        (identity == 'mutant' & step_effect == 'beneficial')),
           peak_color = ifelse(peak == T, id, 33)) %>%
    ungroup() %>%
    group_by(Genotype) %>%
    mutate(not_shared_peak = any(peak == T) & any(peak == F),
           peak_color = ifelse(any(peak == T), id, 33))
  my_comparisons <- list(c('Focal only', 'Focal with HGT to transient'))
  peak_info <- full %>%
    ungroup() %>%
    select(id, peak_color) %>%
    rename(Genotype = id) %>%
    distinct()
  df$data[[1]] <- left_join(df$data[[1]], peak_info)
  
  # wrapped_label <- str_wrap(my_comparisons[[1]][2], width = 15)  # Adjust the width as needed
  plot <- ggplot(data = df$data[[1]], aes(x = Treatment_ID, y = Resistance)) +
    geom_violin(width = 1)+
    geom_jitter(aes(color = as.factor(peak_color)), width = 0.15, size = point_size, alpha = 0.1)+
    stat_summary(fun = "mean", geom = "point", color = "black", size = 3)+
    #stat_summary(aes(color = Treatment_ID), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size, color = 'dark grey', width = 0.15)+
    geom_hline(yintercept = df$Ec_N_mean, linetype = 'dashed', color = 'black')+
    geom_bracket(xmin = 'Focal only', xmax = 'Focal with HGT to transient', y = max(df$data[[1]]$Resistance, na.rm = T) * 1.1, label = df$outcome[1])+
    theme(legend.position = "none") +
    ylab(y_axis_label)+
    xlab(" ")+
    scale_color_manual(values = fill_mapping)+
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, by = 0.25))+
    scale_x_discrete(breaks = c('Focal only', 'Focal with HGT to transient'),
                     labels = c('F', 'F-T'))
  plot
  return(plot)
}
mutation_effects_plot <- function(Pid) {
  sd_value = 0.125
  landscapes <- random_landscape(seed = Pid, sd_value = sd_value)
  mutant_steps <- step_effects(landscapes)
  perturbed_step_effects <- mutant_steps %>%
    select(Species, step, step_slope) %>%
    pivot_wider(names_from = Species, values_from = step_slope) %>%
    mutate(GxH = (Ec > 0 & Kp < 0) | (Ec < 0 & Kp > 0))
  mutant_effects <- ggplot(perturbed_step_effects, aes(x = Ec, y = Kp)) +
    geom_point(aes(color = GxH))+
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey")+
    theme(legend.position = "none")+
    scale_color_manual(values = c('grey', 'black'))+
    scale_y_continuous(limits = c(-0.5, 1))+
    scale_x_continuous(limits = c(-0.5, 1))+
    ylab("mutation effect in transient")+
    xlab("mutation effect in focal")
  return(mutant_effects)
}
landscape_plot <- function(Pid, focal_host) {
  landscapes <- random_landscape(seed = Pid, sd_value = 0.125)
  master.landscape.df <- landscapes %>%
    mutate(num_of_mutations = g + A + E + G + M) 
  steps_file <- read.csv('data/random/mutational_steps.csv', header = T)
  steps <- steps_file %>%
    mutate(step = paste(focal, mutant, sep = "->"))
  focal_node <- steps %>%
    select(-mutant) %>%
    rename(Genotype = focal)
  mutant_step <- steps  %>%
    select(-focal) %>%
    rename(Genotype = mutant)
  focal <- left_join(master.landscape.df, focal_node, multiple = "all")
  focal$identity <- 'focal'
  mutant <- left_join(master.landscape.df, mutant_step, multiple = "all")
  mutant$identity <- 'mutant'
  full_data <- full_join(focal, mutant) 
  full <- full_data %>%
    separate(col = step, sep = '->', into = c('focal', 'mutant'), remove = F) %>%
    filter(!(is.na(focal))) %>%
    group_by(Species, step) %>%
    mutate(focal_MIC = mic[Genotype == focal], 
           mutant_MIC = mic[Genotype == mutant],
           step_slope = (mutant_MIC) - (focal_MIC), 
           step_effect = case_when((mutant_MIC) > (focal_MIC) ~ 'beneficial',
                                   (mutant_MIC) < (focal_MIC) ~ 'deleterious',
                                   TRUE ~ 'neutral')) %>%
    group_by(step) %>%
    mutate(host_epistasis = ifelse(length(unique(step_effect)) == 1, 'noGxH', 'GxH')) %>%
    group_by(Genotype, Species) %>%
    mutate(peak = all((identity == 'focal' & step_effect == 'deleterious')|
                        (identity == 'mutant' & step_effect == 'beneficial')),
           peak_color = ifelse(peak == T, id, 33)) %>%
    ungroup() %>%
    group_by(Genotype) %>%
    mutate(not_shared_peak = any(peak == T) & any(peak == F),
           peak_color = ifelse(any(peak == T), id, 33))
  
  color_mapping <- setNames(c('black','white', 'grey','black'), c(T, F, 'noGxH', 'GxH'))
  
  perturbed_landscape_plot <- ggplot(data = full %>% filter(Species == focal_host), aes(x = num_of_mutations, y = mic, group = step)) +
    geom_line(aes(linetype = step_effect, color = host_epistasis), size = 0.6)+
    geom_point(aes(fill = as.factor(peak_color), color = peak), shape = 21, size = 3, stroke = 0.5) +
    scale_fill_manual(values = fill_mapping)+
    scale_color_manual(values = color_mapping)+
    scale_linetype_manual(values = c('solid','longdash','dotted')) +
    ylab('Fitness')+
    theme(legend.position = "none")+
    xlab('Number of mutations') 
  return(perturbed_landscape_plot)
}
create_partition_plots <- function(df) {
  Figure_width = 8 # in width for one column figure in PNAS
  Figure_height = 8
  landscape_focal <- landscape_plot(df$Partition_ID, 'Ec')
  landscape_transient <- landscape_plot(df$Partition_ID, 'Kp')
  mutant_effects <- mutation_effects_plot(df$Partition_ID)
  evo_outcome <- evo_outcome_plot(df)
  output_dir <- "results/random/4_treatment_comparison/plots/"
  Fig <- plot_grid(landscape_focal, landscape_transient, mutant_effects, evo_outcome, ncol = 2, nrow = 2)  
  save_plot(filename = paste0(output_dir, df$bin, "_", df$Partition_ID, ".pdf"), plot = Fig, 
            base_width = Figure_width, 
            base_height = Figure_height)}

