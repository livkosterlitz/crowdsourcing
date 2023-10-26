# Packages--------
library(matrixcalc)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(ggpubr)
master_df <- NA

# Functions--------
# function used in simulate.gradient.selection function
compute.mics <- function(genotypes, 
                         landscape.df) {
  mics<-numeric(0)
  for(i in 1:length(genotypes)) {
    mics<-c(mics,as.vector(landscape.df[landscape.df$id %in% genotypes[i],]$mic))
  }
  return(mics)
}

# function used in pick.next.genotype
find.neighbors <- function(genotypes,
                           neighborhoods) {
  N.loci <- as.integer(log(nrow(neighborhoods),2))
  neighboring.genotypes<-integer(0)
  for(g in genotypes) {
    for(n in 1:N.loci) {
      neighboring.genotypes<-c(neighboring.genotypes,neighborhoods[g,n])
    }
  }
  return(sort(unique(neighboring.genotypes)))
}

# function used in pick.next.genotype
pick.most.resistant.genotypes<-function(genotypes, 
                                        landscape.df) {
  genotype.mics <- compute.mics(genotypes,landscape.df)
  return(genotypes[which(genotype.mics == max(genotype.mics))])
}
# function used in pick.next.genotype
pick.more.resistant.genotypes <- function(genotypes, 
                                          landscape.df, 
                                          mic.value) {
  genotype.mics <- compute.mics(genotypes,landscape.df)
  return(genotypes[which(genotype.mics > mic.value)])
}

# function used in pick.next.genotype
calculate.prob.templates <- function(candidate.template.genotypes,
                                     candidate.mutant.genotypes,
                                     template.genotypes,
                                     neighborhoods,
                                     N.transformants,
                                     mu) {
  EG.size <- length(candidate.template.genotypes) - 1
  G.size <- length(template.genotypes)
  if(length(candidate.mutant.genotypes)==0) {
    return(rep(1/(1+EG.size),EG.size+1))
  } else {
    product <- 1
    for(m in candidate.mutant.genotypes) {
      NG.size <- length(template.genotypes[template.genotypes %in% neighborhoods[m,]])
      product <- product * (exp(-(NG.size*mu*N.transformants)/G.size))
    }
    return(rep(product/(1+EG.size),EG.size+1))
  }
}

# function used in pick.next.genotype
calculate.prob.mutants<- function(candidate.template.genotypes,
                                  candidate.mutant.genotypes,
                                  template.genotypes,
                                  candidate.template.mics,
                                  candidate.mutant.mics,
                                  neighborhoods,
                                  N.transformants,
                                  mu) {
  CM.size <- length(candidate.mutant.genotypes)
  if(CM.size==0) { #no neighbors that have a higher MIC than the template
    return(numeric(0))
  } else {
    G.size <- length(template.genotypes)
    prob <- rep(0.0,CM.size)
    index <- 1
    for(m in candidate.mutant.genotypes) {
      focal.mic <- candidate.mutant.mics[which(candidate.mutant.genotypes == m)]
      NG.size <- length(template.genotypes[template.genotypes %in% neighborhoods[m,]]) # how many templates is this genotype a neighbor of
      first.factor <- (1 - exp(-(NG.size*mu*N.transformants)/G.size))
      more.resistant.mutants <- candidate.mutant.genotypes[which(candidate.mutant.mics > focal.mic)] # which mutational neighbors also being screened have a higher MIC
      if(length(more.resistant.mutants)==0) {
        second.factor <- 1
      } else {
        second.factor <- 1
        for(m.r in more.resistant.mutants) {
          NG.m.r.size <- length(template.genotypes[template.genotypes %in% neighborhoods[m.r,]])
          second.factor <- second.factor * (exp(-(NG.m.r.size*mu*N.transformants)/G.size))
        }
      }
      equally.resistant.mutants <- remove.element(candidate.mutant.genotypes[which(candidate.mutant.mics == focal.mic)],m)
      EM.size <- length(equally.resistant.mutants)
      if(EM.size==0) {
        third.factor <- 1
      } else {
        tt <- generate.subsets(EM.size)
        sum <- 0
        for(s in 1:nrow(tt)) {
          equally.resistant.mutants.present <- equally.resistant.mutants[tt[s,]]
          equally.resistant.mutants.absent <- equally.resistant.mutants[!tt[s,]]
          S.size <- length(equally.resistant.mutants.present)
          if(S.size == 0) {
            first.subfactor <- 1
          } else {
            first.subfactor <- 1
            for(m.p in equally.resistant.mutants.present) {
              NG.m.p.size <- length(template.genotypes[template.genotypes %in% neighborhoods[m.p,]])
              first.subfactor <- first.subfactor * (1-exp(-(NG.m.p.size*mu*N.transformants)/G.size))
            }  
          }
          if(length(equally.resistant.mutants.absent) == 0) {
            second.subfactor <- 1
          } else {
            second.subfactor <- 1
            for(m.a in equally.resistant.mutants.absent) {
              NG.m.a.size <- length(template.genotypes[template.genotypes %in% neighborhoods[m.a,]])
              second.subfactor <- second.subfactor * (exp(-(NG.m.a.size*mu*N.transformants)/G.size))
            }  
          }
          sum <- sum + (first.subfactor * second.subfactor)/(1+S.size)
        }
        third.factor <- sum
      }
      prob[index] <- first.factor * second.factor * third.factor
      index <- index + 1
    }
    return(prob)
  }
}

# function used in calculate.prob.mutants
remove.element<-function(vector,element) {
  new.vector <- integer(0)
  for(i in 1:length(vector)) {
    if(vector[i] != element) {
      new.vector <- c(new.vector,vector[i])
    }
  }
  return(new.vector)
}

# function used in calculate.prob.mutants
generate.subsets <- function(set.size) {
  m <- matrix(TRUE,nrow=2^set.size, ncol=set.size)
  for(j in 1:set.size) {
    m[,j] <- rep(c(rep(FALSE,2^(j-1)),rep(TRUE,2^(j-1))),2^(set.size-j))
  }
  return(m)
}

# function used in to create neighbor matrix
mutational.distance <- function(v1, v2) {
  L <- length(v1)
  md <- 0
  for(i in 1:L) {
    if(v1[i]!=v2[i]) {
      md <- md + 1
    }
  }
  return(md)
}

# function used in simulate.gradient.selection function
pick.next.genotype<-function(template.genotypes,
                             neighborhoods,
                             landscape.df,
                             N.transformants,
                             mu,
                             N.loci) {
  
  #number of templates
  N.template.genotypes <- length(template.genotypes)
  
  # determine mutant genotypes
  mutant.genotypes <- find.neighbors(template.genotypes,neighborhoods)
  
  # determine candidate template genotypes
  candidate.template.genotypes <- pick.most.resistant.genotypes(template.genotypes, landscape.df)
  candidate.template.mics <- compute.mics(candidate.template.genotypes, landscape.df)
  max.mic.templates <- max(candidate.template.mics)
  
  # determine candidate mutant genotypes
  ### which neighbors have a more resistant genotype than either of the templates
  candidate.mutant.genotypes <- pick.more.resistant.genotypes(mutant.genotypes, landscape.df, max.mic.templates)
  candidate.mutant.mics <- compute.mics(candidate.mutant.genotypes, landscape.df)
  
  # calculate the probabilities of the templates
  prob.candidate.templates <- calculate.prob.templates(candidate.template.genotypes,
                                                       candidate.mutant.genotypes,
                                                       template.genotypes,
                                                       neighborhoods,
                                                       N.transformants,
                                                       mu)
  
  # calculate the probabilities of the mutants
  prob.candidate.mutants <- calculate.prob.mutants(candidate.template.genotypes,
                                                   candidate.mutant.genotypes,
                                                   template.genotypes,
                                                   candidate.template.mics,
                                                   candidate.mutant.mics,
                                                   neighborhoods,
                                                   N.transformants,
                                                   mu)
  
  index <- which(rmultinom(1,1,c(prob.candidate.templates,prob.candidate.mutants))[,1] == 1)
  N.candidate.templates <- length(candidate.template.genotypes)
  if(index<=N.candidate.templates) {
    return(candidate.template.genotypes[index])
  } else {
    return(candidate.mutant.genotypes[index-N.candidate.templates])
  }
}




# Function to calculate alignment score ####
calculate_alignment_score <- function(seed, sd_value) {
  Pid <- seed
  # Construct master landscapes #####
  set.seed(Pid)
  landscapes <- read.csv("data/random/code_format.csv")
  mutations <- c('g', 'A', 'E', 'G', 'M')
  num_of_mutations <- length(mutations)
  num_of_hosts <- length(unique(landscapes$Species))
  mutation_effects <- runif(5, min = 0, max = 0.2)
  master.landscape.df <- landscapes 
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
  ## Analyze mutational steps on master landscape ####
  master.landscape.df <- master.landscape.df %>%
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
  
  host_epistasis <- full %>%
    group_by(step) %>%
    filter(host_epistasis == T)
  
  master_step_effects <- full %>%
    filter(identity == 'focal') %>%
    select(Species, step, step_slope) %>%
    pivot_wider(names_from = Species, values_from = step_slope)
  
  master_landscape_plot <- ggplot(data = full, aes(x = num_of_mutations, y = mic, group = step)) +
    geom_line(aes(linetype = step_effect, color = host_epistasis), size = 0.6)+
    geom_point(aes(fill = Genotype), shape = 21, size = 2, stroke = 0) +
    facet_grid(~Species)+
    scale_color_manual(values = c('grey','black'))+
    scale_linetype_manual(values = c('solid','longdash','dotted')) +
    ylab('Fitness')+
    xlab('Number of mutations')+
    theme(legend.position = "none")
  
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
  ## Analyze mutational steps on perturbed landscapes ####
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
  
  host_epistasis <- full %>%
    group_by(step) %>%
    filter(host_epistasis == T)
  
  perturbed_step_effects <- full %>%
    filter(identity == 'focal') %>%
    select(Species, step, step_slope) %>%
    pivot_wider(names_from = Species, values_from = step_slope)
  perturbed_step_effects$landscape <- 'perturbed'
  master_step_effects$landscape <- 'master'
  step_effects <- rbind(master_step_effects, perturbed_step_effects)
  
  master_step_effects$distance_to_identity <- abs(master_step_effects$Ec - master_step_effects$Kp) / sqrt(2)
  sum_distances_master <- sum(master_step_effects$distance_to_identity)
  
  perturbed_step_effects$distance_to_identity <- abs(perturbed_step_effects$Ec - perturbed_step_effects$Kp) / sqrt(2)
  sum_distances_perturbed <- sum(perturbed_step_effects$distance_to_identity)
  
  mutant_effects <- ggplot(step_effects, aes(x = Ec, y = Kp)) +
    geom_point()+
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey")+
    facet_grid(~landscape)
  
  perturbed_landscape_plot <- ggplot(data = full, aes(x = num_of_mutations, y = mic, group = step)) +
    geom_line(aes(linetype = step_effect, color = host_epistasis), size = 0.6)+
    geom_point(aes(fill = Genotype, color = gen_perturbed), shape = 21, size = 2, stroke = 1) +
    facet_grid(~Species)+
    scale_color_manual(values = c('grey','black'))+
    scale_linetype_manual(values = c('solid','longdash','dotted')) +
    ylab('Fitness')+
    xlab('Number of mutations') +
    theme(legend.position = "none")
  
  combined_plot <- plot_grid(
    master_landscape_plot + labs(title = paste("master landscapes, alignment =", sum_distances_master)),
    perturbed_landscape_plot + labs(title = paste("perturbed landscapes, N =", num_gen_perturb, ", alignment =", round(sum_distances_perturbed, 4))),
    mutant_effects + labs(title = "mutant effects"),
    labels = "AUTO",  # Automatically label each plot with letters (a, b, c)
    ncol = 1)  # Arrange the plots in one column
  
  # Save the combined plot to a PDF file
  ggsave(paste("results/random/partition_landscapes/", Pid, "_", sd_value, "_landscapes", '.pdf', sep = ''), combined_plot, width = 7, height = 10)
  
  return(c(num_gen_perturb, sum_distances_perturbed))
}

# Iterating through seeds #####
num_of_seeds = 2000
seed_df <- data.frame(seeds = seq(num_of_seeds), num_of_gen = rep(NA, num_of_seeds), alignment = rep(NA, num_of_seeds))

for (i in seq(num_of_seeds)){
  outputs <- calculate_alignment_score(i, 0.125)
  seed_df$num_of_gen[seed_df$seeds == i] <- outputs[1]
  seed_df$alignment[seed_df$seeds == i] <- outputs[2]
}
write.csv(seed_df, file = "results/random/seeds_df_summary_sd0.125.csv")

ggplot(data = seed_df, aes(x = alignment, y = num_of_gen))+
  geom_point()

seed_hist <- ggplot(data = seed_df, aes(x = alignment))+
  geom_histogram(binwidth = 0.5)
ggsave("results/random/seed_hist_0.0125.pdf", seed_hist)

# Picking the seeds #####
# Set the range and bin size
range_start <- 0
range_end <- 12
bin_size <- 0.5

# Create bins
bins <- seq(range_start, range_end, by = bin_size)

# Initialize an empty list to store seeds for each bin
seeds_per_bin <- vector("list", length(bins))
names(seeds_per_bin) <- bins


# Loop through random seeds
max_seeds <- 2000  # Adjust as needed

for (seed in 1:max_seeds) {
  alignment_score <- seed_df$alignment[seed]
  
  # Find the bin for the alignment score
  bin_index <- findInterval(alignment_score, bins)
  
  # # Print debugging information
  # cat("Seed:", seed, "Alignment Score:", alignment_score, "Bin Index:", bin_index, "\n")
  # Check if the bin has reached the desired number of seeds
  if (length(seeds_per_bin[[bin_index]]) < 15) {
    seeds_per_bin[[bin_index]] <- c(seeds_per_bin[[bin_index]], seed)
    }
  }

seed_df$bin <- NA
for (bin in bins) {
  for (seed in seeds_per_bin[[as.character(bin)]]) {
    seed_df$bin[seed_df$seeds == seed] <- bin
  }
}

seed_df_filtered <- seed_df %>%
  filter(!is.na(bin))
write.csv(seed_df_filtered, file = "results/random/seeds_df_filtered_sd0.125.csv")

ggplot(data = seed_df_filtered, aes(x = alignment, y = num_of_gen))+
  geom_point()

seed_density <- ggplot(data = seed_df_filtered, aes(x = alignment))+
  geom_density()
ggsave("results/random/seed_density_0.0125_filtered.pdf", seed_density)

# sims #####
seed_df_file <- read.csv("results/random/summary_plots/seeds_df_filtered_sd0.125.csv")
seed_df_ordered <- seed_df_file %>% arrange(desc(alignment))
sd_value = 0.125
bin_run = 9

Partition_ID	<- NA
num_of_gen <- NA
alignment <- NA
bin	<- NA
Focal_species <- NA
Resistance.mean_N <- NA
Resistance.mean_HGT <- NA
Resistance.se_N <- NA
Resistance.se_HGT <- NA
Absolute_difference <- NA
SE_difference <- NA
p_value <- NA
df_master <- tibble(Partition_ID,
                    num_of_gen,
                    alignment,
                    bin,
                    Focal_species,
                    Resistance.mean_N,
                    Resistance.mean_HGT,
                    Resistance.se_N,
                    Resistance.se_HGT,
                    Absolute_difference,
                    SE_difference,
                    p_value)

for (i in seed_df_ordered$seeds[seed_df_ordered$bin == bin_run]){
  Pid <- i 
  # Construct master landscapes #####
  set.seed(Pid)
  landscapes <- read.csv("data/random/code_format.csv")
  mutations <- c('g', 'A', 'E', 'G', 'M')
  num_of_mutations <- length(mutations)
  num_of_hosts <- length(unique(landscapes$Species))
  mutation_effects <- runif(5, min = 0, max = 0.2)
  master.landscape.df <- landscapes 
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
  ## Analyze mutational steps on master landscape ####
  master.landscape.df <- master.landscape.df %>%
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
  
  host_epistasis <- full %>%
    group_by(step) %>%
    filter(host_epistasis == T)
  
  master_step_effects <- full %>%
    filter(identity == 'focal') %>%
    select(Species, step, step_slope) %>%
    pivot_wider(names_from = Species, values_from = step_slope)
  
  master_landscape_plot <- ggplot(data = full, aes(x = num_of_mutations, y = mic, group = step)) +
    geom_line(aes(linetype = step_effect, color = host_epistasis), size = 0.6)+
    geom_point(aes(fill = Genotype), shape = 21, size = 2, stroke = 0) +
    facet_grid(~Species)+
    scale_color_manual(values = c('grey','black'))+
    scale_linetype_manual(values = c('solid','longdash','dotted')) +
    ylab('Fitness')+
    xlab('Number of mutations')+
    theme(legend.position = "none")
  
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
  ## Analyze mutational steps on perturbed landscapes ####
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
  
  host_epistasis <- full %>%
    group_by(step) %>%
    filter(host_epistasis == T)
  
  perturbed_step_effects <- full %>%
    filter(identity == 'focal') %>%
    select(Species, step, step_slope) %>%
    pivot_wider(names_from = Species, values_from = step_slope)
  perturbed_step_effects$landscape <- 'perturbed'
  master_step_effects$landscape <- 'master'
  step_effects <- rbind(master_step_effects, perturbed_step_effects)
  
  master_step_effects$distance_to_identity <- abs(master_step_effects$Ec - master_step_effects$Kp) / sqrt(2)
  sum_distances_master <- sum(master_step_effects$distance_to_identity)
  
  perturbed_step_effects$distance_to_identity <- abs(perturbed_step_effects$Ec - perturbed_step_effects$Kp) / sqrt(2)
  sum_distances_perturbed <- sum(perturbed_step_effects$distance_to_identity)
  
  mutant_effects <- ggplot(step_effects, aes(x = Ec, y = Kp)) +
    geom_point()+
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey")+
    facet_grid(~landscape)
  
  perturbed_landscape_plot <- ggplot(data = full, aes(x = num_of_mutations, y = mic, group = step)) +
    geom_line(aes(linetype = step_effect, color = host_epistasis), size = 0.6)+
    geom_point(aes(fill = Genotype, color = gen_perturbed), shape = 21, size = 2, stroke = 1) +
    facet_grid(~Species)+
    scale_color_manual(values = c('grey','black'))+
    scale_linetype_manual(values = c('solid','longdash','dotted')) +
    ylab('Fitness')+
    xlab('Number of mutations') +
    theme(legend.position = "none")
  
  combined_plot <- plot_grid(
    master_landscape_plot + labs(title = paste("master landscapes, alignment =", sum_distances_master)),
    perturbed_landscape_plot + labs(title = paste("perturbed landscapes, N =", num_gen_perturb, ", alignment =", round(sum_distances_perturbed, 4))),
    mutant_effects + labs(title = "mutant effects"),
    labels = "AUTO",  # Automatically label each plot with letters (a, b, c)
    ncol = 1)  # Arrange the plots in one column
  
  # Save the combined plot to a PDF file
  ggsave(paste("results/random/partition_landscapes/", Pid, "_", sd_value, "_landscapes", '.pdf', sep = ''), combined_plot, width = 7, height = 10)
  # Creation of the mutational neighbor matrix -------
  Ecoli.landscape.df <- landscapes %>%
    filter(Species == 'Ec')
  N.g<-length(Ecoli.landscape.df$id) # number of genotypes
  N.l<-as.integer(log(N.g,2)) # number of loci
  neighbor_hoods<-matrix(0,nrow=N.g,ncol=N.l) # neighbor matrix for each genotype created in the for loop below
  for(f.gen in 1:N.g) { #loop through the number of genotypes
    col.counter<-1
    for(o.gen in 1:N.g) { #loop through the number of genotypes
      # compare the two genotypes, if they are one mutation away from eachother proceed
      if(mutational.distance(as.integer(Ecoli.landscape.df[f.gen,2:(1+N.l)]),
                             as.integer(Ecoli.landscape.df[o.gen,2:(1+N.l)]))==1) {
        neighbor_hoods[f.gen,col.counter] <- o.gen # add the genotype to the neighbor matrix in one of the columns
        col.counter <- col.counter + 1}}}
  # Build the partition data frame --------
  landscapes <- master.landscape.df
  treatment_csv <- read.csv("data/random/Treatment_random.csv")
  Partition_ID <- NA
  Treatment_ID <- NA
  Replicate <- NA
  Cumulative_time <- NA
  Species_time <- NA
  Species_order <- NA
  N <- NA
  Resistance <- NA
  Genotype <- NA
  Resistance.mean <- NA
  Resistance.sd <- NA
  Resistance.se <- NA
  partition_endpoint_average <- tibble(Partition_ID,
                                       Treatment_ID,
                                       Cumulative_time,
                                       Species_time,
                                       Species_order,
                                       N,
                                       Resistance.mean,
                                       Resistance.sd,
                                       Resistance.se)
  partition_endpoints <- tibble(Partition_ID,
                                       Treatment_ID,
                                       Replicate,
                                       Cumulative_time,
                                       Species_time,
                                       Species_order,
                                       Genotype,
                                       Resistance)
  # Run the partitions simulations --------
  for(id in treatment_csv$Treatment_ID){
    id<-treatment_csv$Treatment_ID==id
    partition_ID_x <- Pid
    treatment_ID_x <- treatment_csv$Treatment_ID[id]
    Num_of_reps_x <- treatment_csv$Rep[id]
    Cumulative_time_x <- treatment_csv$Cumulative_time[id]
    Species_time_x <- as.numeric(unlist(as.vector(strsplit(treatment_csv$Species_time[id], ','))))
    Species_order_x <- unlist(strsplit(treatment_csv$Species[id], ','))
    N_transformants_x <- treatment_csv$N_transformants[id]
    mu_x <- treatment_csv$mu[id]
    Partition_ID <- NA
    Treatment_ID <- NA
    Replicate <- NA
    Cumulative_time <- NA
    Species_time <- NA
    Species_order <- NA
    Genotype <- NA
    g <- NA
    A <- NA
    E <- NA
    G <- NA
    M <- NA
    Resistance <- NA
    treatment_df <- tibble(Partition_ID,
                           Treatment_ID,
                           Replicate,
                           Cumulative_time,
                           Species_time,
                           Species_order,
                           Genotype,
                           g,
                           A,
                           E,
                           G,
                           M,
                           Resistance)
    if (sum(Species_time_x) != Cumulative_time_x) {
      print(F)
    }
    Cumulative_time <- seq(Cumulative_time_x)
    Species_time <- NULL
    Species_order <- NULL
    # make time vector and species ordering vector according to the treatment
    for (s in seq(length(Species_order_x))) {
      Species_order <- c(Species_order, rep(Species_order_x[s], Species_time_x[s]))
      Species_time <- c(Species_time, seq(Species_time_x[s]))
    }
    # make an empty data frame to store the replicate simulation information for each treatment
    for (r in seq(Num_of_reps_x)) {
      Partition_ID <- rep(partition_ID_x, Cumulative_time_x)
      Treatment_ID <- rep(treatment_ID_x, Cumulative_time_x)
      Replicate <- rep(r, Cumulative_time_x)
      Genotype <- rep(NA, Cumulative_time_x)
      Resistance <- rep(NA, Cumulative_time_x)
      df <- tibble(Partition_ID,
                   Treatment_ID,
                   Replicate,
                   Cumulative_time,
                   Species_time,
                   Species_order,
                   Genotype,
                   g,
                   A,
                   E,
                   G,
                   M,
                   Resistance)
      treatment_df <- rbind(treatment_df, df)
    }
    treatment_df = treatment_df[-1,] #removes the first row with NA values
    # Running the sim --------
    initial.genotype <- 1
    final.genotype <- 32
    N.transformants <- N_transformants_x
    mu <- mu_x
    N.loci <- 5

    # run the number of simulations for the number of replicates specified for the treatment
    for(r in 1:max(treatment_df$Replicate)) {
      row_marker = treatment_df$Replicate == r & treatment_df$Cumulative_time == 1
      treatment_df$Genotype[row_marker] <- initial.genotype
      landscape_row_marker <- landscapes$Species == treatment_df$Species_order[row_marker] & landscapes$id == treatment_df$Genotype[row_marker]
      n.row<-1:nrow(landscapes)
      # mic<-numeric(0)
      # for(n in n.row) {
      #   mic<-c(mic,landscapes[n,sample(7:8,1)])
      # }
      # landscapes$mic <- mic
      treatment_df$Resistance[row_marker] <- landscapes$mic[landscape_row_marker]
      treatment_df$g[row_marker] <- landscapes$g[landscape_row_marker]
      treatment_df$A[row_marker] <- landscapes$A[landscape_row_marker]
      treatment_df$E[row_marker] <- landscapes$E[landscape_row_marker]
      treatment_df$G[row_marker] <- landscapes$G[landscape_row_marker]
      treatment_df$M[row_marker] <- landscapes$M[landscape_row_marker]
      for (t in 2:max(treatment_df$Cumulative_time)) {
        row_marker_previous = row_marker
        row_marker = treatment_df$Replicate == r & treatment_df$Cumulative_time == t
        template.genotype = treatment_df$Genotype[row_marker_previous]
        # mic<-numeric(0)
        # for(n in n.row) {
        #   #### change this! This will just be mic
        #   mic<-c(mic,landscapes[n,sample(7:8,1)])
        # }
        # landscapes$mic <- mic
        landscape.focal = landscapes %>% filter(Species == treatment_df$Species_order[row_marker])
        treatment_df$Genotype[row_marker] <- pick.next.genotype(template.genotype,
                                                                neighbor_hoods,
                                                                landscape.focal,
                                                                N.transformants,
                                                                mu,
                                                                N.loci)
        landscape_row_marker <- landscapes$Species == treatment_df$Species_order[row_marker] & landscapes$id == treatment_df$Genotype[row_marker]
        treatment_df$Resistance[row_marker] <- landscapes$mic[landscape_row_marker]
        treatment_df$g[row_marker] <- landscapes$g[landscape_row_marker]
        treatment_df$A[row_marker] <- landscapes$A[landscape_row_marker]
        treatment_df$E[row_marker] <- landscapes$E[landscape_row_marker]
        treatment_df$G[row_marker] <- landscapes$G[landscape_row_marker]
        treatment_df$M[row_marker] <- landscapes$M[landscape_row_marker]
      }
    }

    # summarise all of the replicate simulations for the treatment by computing averages
    treatment_df_summary <- treatment_df %>%
      group_by(Partition_ID, Treatment_ID, Cumulative_time, Species_time, Species_order) %>%
      summarise(N = n(),
                Resistance.mean = mean(Resistance),
                Resistance.sd = sd(Resistance),
                Resistance.se = Resistance.sd/sqrt(N))
    treatment_df_endpoints <- treatment_df %>%
      select(-g, -A, -E, -G, -M) %>%
      filter(Cumulative_time == max(Cumulative_time))
    partition_endpoints <- rbind(partition_endpoints, treatment_df_endpoints)  
    partition_endpoint_average <- rbind(partition_endpoint_average, treatment_df_summary[nrow(treatment_df_summary),])
    write.csv(x = treatment_df, file = paste("results/random/partition_details/", Pid, "_", Treatment_ID[1], '.csv', sep = ''))
    write.csv(x = treatment_df_summary, file = paste("results/random/partition_details/", Pid, "_",Treatment_ID[1], '_averages.csv', sep = ''))
  }
  partition_endpoints = partition_endpoints[-1,]
  partition_endpoint_average = partition_endpoint_average[-1,]
  partition_endpoints$Treatment_ID <- factor(partition_endpoints$Treatment_ID)
  red_Ec_color <- rgb(255, 48, 48, maxColorValue = 255)
  blue_Kp_color <- rgb(0, 0, 238, maxColorValue = 255)
  purple_Ec.Kp_color <- rgb(178, 58, 238, maxColorValue = 255)
  color_set = c(purple_Ec.Kp_color, red_Ec_color, purple_Ec.Kp_color, blue_Kp_color)
  my_comparisons <- list(c("Ec_N", "Ec_HGT"), c("Kp_N", "Kp_HGT"))
  MinMeanSEMMax <- function(x) {
    v <- c(min(x), 
           mean(x) - sd(x)/sqrt(length(x)), 
           mean(x), 
           mean(x) + sd(x)/sqrt(length(x)), 
           max(x))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
  }
  # Statistics plot -------------------------------------------------------
  text_size_axis_title = 9
  text_size_axis_tick_labels = 7
  axis_line_size = 0.75/(ggplot2::.pt*72.27/96) # pt size for the axis lines
  axis_tick_size = 0.75/(ggplot2::.pt*72.27/96)
  axis_tick_lengths = 0.05
  y_axis_label = expression(paste("resistance"))
  point_jitter_width = 0.15
  point_size = 1
  margins<- c(0.1,0.1,0.1,0.1)
  box_line_size = 0.75/(ggplot2::.pt*72.27/96) 
  Ec.mean <- partition_endpoint_average %>% filter(Cumulative_time == 60 & Treatment_ID == 'Ec_N') %>% select(Resistance.mean)
  Kp.mean <- partition_endpoint_average %>% filter(Cumulative_time == 60 & Treatment_ID == 'Kp_N') %>% select(Resistance.mean)
  E <- ggplot(data = partition_endpoints, 
              aes(x=Treatment_ID, y=Resistance, color = Treatment_ID))+
    stat_summary(data = partition_endpoints, mapping = aes(Treatment_ID, Resistance, color = Treatment_ID), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
    geom_jitter(width = point_jitter_width, size = point_size, alpha = 0.1) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", aes(label = ..p.signif..), method.args = list(p.adjust.method = "bonferroni")) +
    geom_hline(yintercept = Ec.mean$Resistance.mean, linetype = 'dashed', color = red_Ec_color)+
    geom_hline(yintercept = Kp.mean$Resistance.mean, linetype = 'dashed', color = blue_Kp_color)+
    ylab(y_axis_label)+
    theme(axis.title.x = element_blank()) +
    theme(axis.title = element_text(size = text_size_axis_title)) +
    theme(axis.line.y = element_line(size = axis_line_size)) +
    theme(axis.line.x = element_line(size = axis_line_size)) +
    theme(axis.ticks = element_line(size = axis_tick_size))+
    theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
    theme(axis.text.y = element_text(size = text_size_axis_tick_labels))+
    theme(axis.text.x = element_text(size = text_size_axis_tick_labels))+
    theme(legend.position = "none") +
    theme(plot.margin = unit(margins, "in"))+
    scale_color_manual(values = color_set, guide = "none")+
    theme(strip.background = element_blank())+
    theme(strip.text = element_blank())
  ggsave(file = paste("results/random/", Pid, "_statistics", '.pdf', sep = ''), plot = E)
  
  # Partition analysis ####
  write.csv(x = partition_endpoint_average, file = paste("results/random/", Pid, "_endpoints", '.csv', sep = ''))
  partition_analysis <- partition_endpoint_average %>%
    separate(Treatment_ID, into = c("Focal_species", "Treatment"), sep = "_", remove = FALSE) %>%
    select(Partition_ID, Focal_species, Treatment, Resistance.mean, Resistance.se) %>%
    pivot_wider(names_from = Treatment, values_from = c(Resistance.mean, Resistance.se)) %>%
    mutate(Absolute_difference = abs(Resistance.mean_N - Resistance.mean_HGT),
           SE_difference = sqrt(Resistance.se_N^2 + Resistance.se_HGT^2))
  partition_analysis$p_value <- rep(NA, nrow(partition_analysis))
  partition_analysis$p_value[partition_analysis$Focal_species == "Ec"] <- compare_means(Resistance ~ Treatment_ID, data = partition_endpoints %>% filter(Species_order == 'Ec'), ref.group = "Ec_N", method = "wilcox.test", p.adjust.method = "bonferroni")$p
  partition_analysis$p_value[partition_analysis$Focal_species == "Kp"] <- compare_means(Resistance ~ Treatment_ID, data = partition_endpoints %>% filter(Species_order == 'Kp'), ref.group = "Kp_N", method = "wilcox.test", p.adjust.method = "bonferroni")$p
  partition_analysis$num_of_gen <- rep(seed_df_ordered$num_of_gen[seed_df_ordered$seeds==i], nrow(partition_analysis))
  partition_analysis$bin <- rep(seed_df_ordered$bin[seed_df_ordered$seeds==i], nrow(partition_analysis))
  partition_analysis$alignment <- rep(seed_df_ordered$alignment[seed_df_ordered$seeds==i], nrow(partition_analysis))
  df_master <- rbind(df_master, partition_analysis)
  write.csv(x = df_master, file = paste("results/random/", "master", bin_run, '.csv', sep = ''))
}











#master_df <- master_df[-1,]
master_df_plot <- master_df %>%
  separate(Treatment_ID, into = c("Focal_species", "Treatment"), sep = "_", remove = FALSE) %>%
  select(Partition_ID, Focal_species, Treatment, Resistance.mean, Resistance.se, alignment_score) %>%
  pivot_wider(names_from = Treatment, values_from = c(Resistance.mean, Resistance.se)) %>%
  mutate(Absolute_difference = abs(Resistance.mean_N - Resistance.mean_HGT),
         SE_difference = sqrt(Resistance.se_N^2 + Resistance.se_HGT^2))
#write.csv(x = master_df_plot, file = "results/random/random_plots/2023-08-08_masterdf_plot.csv")

# plots ####
# Create density plot
density_plot <- ggplot(data = master_df_plot, aes(x = alignment_score)) +
  geom_density(alpha = 0.5)+
  theme_minimal()

# Create scatter plot
explanatory_text <- "N = 150\nx-axis is the alignment score calculated by summing the distances of data points from the identity line (y = x)\n where the data points are the effects of mutations in each host landscape.\ny-axis is the absolute difference in the average endpoint fitness of the HGT compared to the non-HGT\n treatment from 500 greedy walk simulations for each treatment"
scatter_plot <- ggplot(data = master_df_plot, aes(x = alignment_score, y = Absolute_difference)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  stat_cor(method = "pearson", label.x = 0.1, label.y = 0.04) +
  geom_errorbar(aes(ymin = Absolute_difference - SE_difference, ymax = Absolute_difference + SE_difference), width = 0.2) +
  annotate("text", x = 0, y = max(master_df_plot$Absolute_difference), label = explanatory_text, vjust = 0.1, hjust = 0) +
  theme_minimal()

# Combine plots side by side
combined_plot <- plot_grid(scatter_plot, density_plot, labels = c("Density Plot", "Scatter Plot"), ncol = 1)
ggsave("results/random/random_plots/2023-08-08_alignment_abs_w_density.pdf", combined_plot, width = 10, height = 20)

# Print the combined plot
print(combined_plot)

