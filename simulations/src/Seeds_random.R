# Packages--------
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
source('src/Random_Landscape_Functions.R')

# Initializing the bins #####
range_start <- 0
range_end <- 15
bin_size <- 0.5
bins <- seq(range_start, range_end, by = bin_size)

# Iterating through seeds #####
num_of_seeds = 10
sd_value = 0.125
seed_df <- data.frame(seeds = seq(num_of_seeds))

for (i in seq(num_of_seeds)){
  landscape <- random_landscape(i, sd_value)
  outputs <- calculate_alignment_score(landscape)
  seed_df$num_of_gen[seed_df$seeds == i] <- sum(landscape$gen_perturbed)
  seed_df$alignment_perp_dist_all[seed_df$seeds == i] <- outputs[1]
  seed_df$bin_perp_dist_all[seed_df$seeds == i] <- bins[findInterval(outputs[1], bins)]
  seed_df$alignment_perp_dist_signeffects[seed_df$seeds == i] <- outputs[2]
  seed_df$alignment_proportion_signeffects[seed_df$seeds == i] <- outputs[3]
}

write.csv(seed_df, file = "results/random/1_random_seeds/seed_alignment_scores.csv")