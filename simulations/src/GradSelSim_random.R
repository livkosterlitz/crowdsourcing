# Packages--------
library(matrixcalc)
library(tidyverse)
library(cowplot)
library(argparse)
source("src/GradSelSim_Functions.R")
source('src/Random_Landscape_Functions.R')
theme_set(theme_cowplot())

## Command line ------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("-s", "--seed", help="the seed used to generate the landscapes")
args <- parser$parse_args()
Pid <- args$seed

## Troubleshooting variables ---------------------------------------------------------------
#Pid <- 875

# Generate landscapes ------
sd_value = 0.125
landscapes <- random_landscape(Pid, sd_value)
bin <- find_bin(landscapes)

# Creation of the mutational neighbor matrix -------
Ecoli.landscape.df <- landscapes %>%
  filter(Species == 'Ec')
N.g<-length(Ecoli.landscape.df$id) # number of genotypes
N.l<-as.integer(log(N.g,2)) # number of loci
neighbor_hoods <- mutant_neighborhoods(N.g, #number of genotypes
                                       N.l, #number of loci
                                       Ecoli.landscape.df) #this landscape entry must hold the binary genotype information starting in column 2 + column 2+N.l
# Build the partition data frame --------
suppressWarnings({
treatment_csv <- read.csv("data/random/Treatment_random.csv")
})
# Run the partitions simulations for the HGT and no HGT treatments --------
for(id in treatment_csv$Treatment_ID){ #will run both the Ec_HGT and Ec_N treatements
  id<-treatment_csv$Treatment_ID==id
  partition_ID_x <- Pid
  treatment_ID_x <- treatment_csv$Treatment_ID[id]
  Num_of_reps_x <- treatment_csv$Rep[id]
  Cumulative_time_x <- treatment_csv$Cumulative_time[id]
  Species_time_x <- as.numeric(unlist(as.vector(strsplit(treatment_csv$Species_time[id], ','))))
  Species_order_x <- unlist(strsplit(treatment_csv$Species[id], ','))
  N_transformants_x <- treatment_csv$N_transformants[id]
  mu_x <- treatment_csv$mu[id]
  column_names <- c('Partition_ID', # Create a list with column names and initial values
                    'Treatment_ID', 
                    'Replicate', 
                    'Cumulative_time', 
                    'Species_time', 
                    'Species_order',
                    'Genotype',
                    'g',
                    'A',
                    'E',
                    'G',
                    'M',
                    'Resistance')
  column_values <- setNames(rep(list(NA), length(column_names)), column_names) # Create a named list of columns with NA values
  treatment_df <- as_tibble(column_values) # Create the tibble
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
    g <- rep(NA, Cumulative_time_x)
    A <- rep(NA, Cumulative_time_x)
    E <- rep(NA, Cumulative_time_x)
    G <- rep(NA, Cumulative_time_x)
    M <- rep(NA, Cumulative_time_x)
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
  # output the replicate simulations for the treatment
  write.csv(x = treatment_df, file = paste("results/random/3_simulations/", bin, "_", Pid, "_", Treatment_ID[1], '.csv', sep = ''))
}
