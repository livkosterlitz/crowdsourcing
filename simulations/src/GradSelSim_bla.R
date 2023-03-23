# Packages--------
library(matrixcalc)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

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
neighbor_hoods<-matrix(0,nrow=N.g,ncol=N.l) # neighbor matrix for each genotype created in the for loop below
for(f.gen in 1:N.g) { #loop through the number of genotypes
  col.counter<-1 
  for(o.gen in 1:N.g) { #loop through the number of genotypes
    # compare the two genotypes, if they are one mutation away from eachother proceed
    if(mutational.distance(as.integer(Ecoli.landscape.df[f.gen,2:(1+N.l)]),
                           as.integer(Ecoli.landscape.df[o.gen,2:(1+N.l)]))==1) {
      neighbor_hoods[f.gen,col.counter] <- o.gen # add the genotype to the neighbor matrix in one of the columns
      col.counter <- col.counter + 1}}}

# Build the treatment data frame --------
treatment_csv <- read.csv("data/bla/Treatment_master_bla.csv")
Treatment_ID <- NA
Cumulative_time <- NA
Species_time <- NA
Species_order <- NA
N <- NA
Resistance.mean <- NA
Resistance.sd <- NA
Resistance.se <- NA
treatment_average <- tibble(Treatment_ID,
                            Cumulative_time,
                            Species_time,
                            Species_order,
                            N,
                            Resistance.mean,
                            Resistance.sd,
                            Resistance.se)
for(id in treatment_csv$Treatment_ID){
  id<-treatment_csv$Treatment_ID==id
  treatment_ID_x <- treatment_csv$Treatment_ID[id]
  Num_of_reps_x <- treatment_csv$Rep[id]
  Cumulative_time_x <- treatment_csv$Cumulative_time[id]
  Species_time_x <- as.numeric(unlist(as.vector(strsplit(treatment_csv$Species_time[id], ','))))
  Species_order_x <- unlist(strsplit(treatment_csv$Species[id], ','))
  N_transformants_x <- treatment_csv$N_transformants[id]
  mu_x <- treatment_csv$mu[id]
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
  treatment_df <- tibble(Treatment_ID, 
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
    Treatment_ID <- rep(treatment_ID_x, Cumulative_time_x)
    Replicate <- rep(r, Cumulative_time_x)
    Genotype <- rep(NA, Cumulative_time_x)
    Resistance <- rep(NA, Cumulative_time_x)
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
