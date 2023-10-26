# Packages ------
library(tidyverse)

# Combine data
combined_data <- data.frame()
folder <- "results/random/3_simulations/"
file_list <- list.files(path = folder, pattern = "*.csv")
for (file in file_list) {
  df <- read.csv(paste(folder, file, sep = ""), header = TRUE)
  df <- df[, -1]
  df <- df %>%
    filter(Cumulative_time == 60)
  bin <- as.numeric(strsplit(file, "_")[[1]][1])
  df$bin <- bin
  combined_data <- bind_rows(combined_data, df)
}

write.csv(combined_data, "results/random/4_treatment_comparison/combined_data.csv", row.names = FALSE)
