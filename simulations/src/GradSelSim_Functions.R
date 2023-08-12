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

# function for mutational neighbor matrix -------
mutant_neighborhoods <- function(N.g,
                                 N.l,
                                 landscape.df){
  neighbor_hoods<-matrix(0,nrow=N.g,ncol=N.l) # neighbor matrix for each genotype created in the for loop below
  for(f.gen in 1:N.g) { #loop through the number of genotypes
    col.counter<-1 
    for(o.gen in 1:N.g) { #loop through the number of genotypes
      # compare the two genotypes, if they are one mutation away from eachother proceed
      if(mutational.distance(as.integer(landscape.df[f.gen,2:(1+N.l)]),
                             as.integer(landscape.df[o.gen,2:(1+N.l)]))==1) {
        neighbor_hoods[f.gen,col.counter] <- o.gen # add the genotype to the neighbor matrix in one of the columns
        col.counter <- col.counter + 1}}}
  return(neighbor_hoods)
}


