# Packages--------
library(matrixcalc)
library(tidyverse)
library(cowplot)
source("src/GradSelSim_Functions.R")
theme_set(theme_cowplot())

# Load landscapes --------
landscape_file <- read.csv("data/bla/MIC_simulation_file_3barcodes.csv")
landscape_file$Species<-recode(landscape_file$Species, `Ec` = "E.coli", `Kp` = "Kleb")
landscape_filter <- landscape_file %>%
  rename(hosts = Species) %>%
  select(-Shorthand, -outlier)
format <- read.csv("data/bla/code_format.csv")
landscapes <- left_join(format, landscape_filter)

id<-1:32
Ecoli.landscape.df<-cbind(id,landscapes[landscapes$hosts=="E.coli",2:8])
Kleb.landscape.df<-cbind(id,landscapes[landscapes$hosts=="Kleb",2:8])
Se.landscape.df<-cbind(id,landscapes[landscapes$hosts=="Se",2:8])
Ecoli.landscape.df$Species <- 'Ec'
Kleb.landscape.df$Species <- 'Kp'
Se.landscape.df$Species <- 'Se'
landscapes <- rbind(Ecoli.landscape.df, Kleb.landscape.df, Se.landscape.df)
# Creation of the mutational neighbor matrix -------
N.g<-length(Ecoli.landscape.df$id) # number of genotypes
N.l<-as.integer(log(N.g,2)) # number of loci
neighbor_hoods <- mutant_neighborhoods(N.g, #number of genotypes
                                       N.l, #number of loci
                                       Ecoli.landscape.df) #this landscape entry must hold the binary genotype information starting in column 2 + column 2+N.l

# Build the treatment data frame --------
treatment_csv <- read.csv("data/bla/Treatment_test.csv")
column_names <- c('Treatment_ID', # Create a list with column names and initial values
                  'Cumulative_time', 
                  'Species_time', 
                  'Species_order', 
                  'N', 
                  'Resistance.mean', 
                  'Resistance.sd', 
                  'Resistance.se')
column_values <- setNames(rep(list(NA), length(column_names)), column_names) # Create a named list of columns with NA values
treatment_average <- as_tibble(column_values) # Create the tibble

for(id in treatment_csv$Treatment_ID){
  id<-treatment_csv$Treatment_ID==id
  treatment_ID_x <- treatment_csv$Treatment_ID[id]
  Num_of_reps_x <- treatment_csv$Rep[id]
  Cumulative_time_x <- treatment_csv$Cumulative_time[id]
  Species_time_x <- as.numeric(unlist(as.vector(strsplit(treatment_csv$Species_time[id], ','))))
  Species_order_x <- unlist(strsplit(treatment_csv$Species[id], ','))
  N_transformants_x <- treatment_csv$N_transformants[id]
  mu_x <- treatment_csv$mu[id]
  column_names <- c('Treatment_ID', # Create a list with column names and initial values
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
    Treatment_ID <- rep(treatment_ID_x, Cumulative_time_x)
    Replicate <- rep(r, Cumulative_time_x)
    Genotype <- rep(NA, Cumulative_time_x)
    Resistance <- rep(NA, Cumulative_time_x)
    g <- rep(NA, Cumulative_time_x)
    A <- rep(NA, Cumulative_time_x)
    E <- rep(NA, Cumulative_time_x)
    G <- rep(NA, Cumulative_time_x)
    M <- rep(NA, Cumulative_time_x)
    df <- tibble(Treatment_ID, 
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
    mic<-numeric(0)
    for(n in n.row) {
      mic<-c(mic,landscapes[n,sample(7:8,1)])
    }
    landscapes$mic <- mic
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
      mic<-numeric(0)
      for(n in n.row) {
        mic<-c(mic,landscapes[n,sample(7:8,1)])
      }
      landscapes$mic <- mic
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
    group_by(Treatment_ID, Cumulative_time, Species_time, Species_order) %>%
    summarise(N = n(),
              Resistance.mean = mean(Resistance), 
              Resistance.sd = sd(Resistance),
              Resistance.se = Resistance.sd/sqrt(N))
  treatment_average <- rbind(treatment_average, treatment_df_summary)
  
  write.csv(x = treatment_df, file = paste("results/bla/", Treatment_ID[1], '.csv', sep = ''))
  write.csv(x = treatment_df_summary, file = paste("results/bla/", Treatment_ID[1], '_averages.csv', sep = ''))
}
