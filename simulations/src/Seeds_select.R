# Packages--------
library(tidyverse)
library(cowplot)
library(stringr)
theme_set(theme_cowplot())

# Filter seeds -----
seed_df <- read.csv("results/random/1_random_seeds/seed_alignment_scores.csv")
bin_count <- 50
selected_seed_df <- seed_df %>%
  arrange(bin_perp_dist_all, seeds) %>%
  group_by(bin_perp_dist_all) %>%
  mutate(row_num = row_number()) %>%
  mutate(selected_seeds = if_else(n() >= bin_count & row_num <= bin_count, T, F)) %>%
  select(-row_num)

# Plot overview --------
bin_overview <- ggplot(data = selected_seed_df, aes(x = bin_perp_dist_all, fill = selected_seeds))+
                geom_bar()+
                geom_hline(yintercept = bin_count, color = 'black')
ggsave("results/random/2_selected_seeds/bin_overview.pdf", bin_overview)

# Generate step 3 shell script -----
# Check if the simulation folder exists
if (file.exists("results/random/3_simulations")) {
  selected_seeds <- selected_seed_df %>%
    filter(selected_seeds) %>%
    pull(seeds)
  existing_files <- list.files("results/random/3_simulations")
  existing_seeds <- str_extract(existing_files, "(?<=_)(\\d+)(?=_)") %>% as.integer()
  new_seeds <- setdiff(selected_seeds, existing_seeds)
  cat("Number of new seeds to be run:", length(new_seeds), "\n")
  cat("Number of seeds that already exist:", length(unique(existing_seeds)), "\n")
  if (length(new_seeds) > 0) {
    cat("Creating step3_runall.sh\n")
    fileConn <- file("results/random/2_selected_seeds/step3_runall.sh", open = "a")  # Open in append mode
    cat("parallel -j 7 <<EOF\n", file = fileConn)
    for (seed in new_seeds) {
      cat("time Rscript src/GradSelSim_random.R -s", seed, "> /dev/null 2>&1 && echo \"seed", seed, "finished\"\n", file = fileConn)
    }
    cat("EOF\n", file = fileConn)
    close(fileConn)
    cat("step3_runall.sh created.\n")
  } else {
    cat("No new seeds to run.\n")
  }
}
