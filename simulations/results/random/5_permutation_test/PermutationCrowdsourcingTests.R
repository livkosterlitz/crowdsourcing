#FUNCTIONS

#This function takes a dataframe and 
#puts together a vector for the bins
#(of misalignment scores) and a vector 
#for the proportion of crowdsourcing 
#outcomes in each bin. These vectors
#are returned as a two-element list.
#-------------------------------------
find.vectors<-function(data) {
  bins<-sort(unique(data$bin))
  outs<-rep(0,length(bins))
  index.counter <- 1
  for(b in bins) {
    counter<-0
    for(o in data[data$bin==b,]$outcome) {
      if(o=="crowdsourcing") {
        counter <- counter + 1
      }
    }
    outs[index.counter] <- (counter/length(data[data$bin==b,]$outcome))
    index.counter <- index.counter + 1
  }
  return(list(bins,outs))
}
#-------------------------------------

#This function returns the slope of 
#a best-fit line to a scatterplot
#given by two input vectors.
#-------------------------------------
find.slope<-function(bin.vec, prop.vec) {
  mod <- lm(prop.vec ~ bin.vec)
  return(as.double(mod$coefficients[2]))
}
#-------------------------------------

#This function returns the p-value
#of a permutation test. Specifically,
#it gives the fraction of slopes 
#(calculated from permuted data)
#that are less than the actual slope.
#(Here, the alternative hypothesis to
#our null of zero slope is a negative
#slope-- i.e., the proportion of 
#crowdsourcing decreases with 
#increasing misalignment). The
#"allow.ties" argument is a bool that,
#when TRUE, requires that a slope
#from permuted data is greater than
#*or equal to* the actual slope to
#contribute to the count used for 
#the p-value. When FALSE, the slope 
#from permuted data must be strictly 
#less than the actual slope.
#-------------------------------------
find.p.value<-function(actual.slope, slopes, allow.ties) {
  counter <- 0
  n.perms <- length(slopes)
  for(p in 1:n.perms) {
    if(allow.ties) {
      if(slopes[p] <= actual.slope) {
        counter <- counter + 1
      }
    } else {
      if(slopes[p] < actual.slope) {
        counter <- counter + 1
      }
    }
  }
  return(counter/n.perms)
}
#-------------------------------------

#This function contains the heart of the
#code to perform a permutation test.
#Here the outcomes are randomly permuted
#across all the runs (note there are multiple
#runs per bin). For each permutation, a
#slope of the within-bin proportion of 
#crowdsourcing against the misalignment 
#scores (which define the bins) is computed.
#This is done "n.perms" times. To get
#the status of permutations being executed,
#"print.status" can be set to TRUE. Note, 
#this is written as the difference between
#the permuted slope and the actual slope
#(the slope of the original data) so that
#the user can see whether the actual slope
#is extreme (this is so when most to all of 
#these differences are positive). If the
#user wishes, a histogram giving the slopes
#resulting from the permutations can be
#displayed (display.graphics=TRUE) where the
#actual slope should appear as a red dot
#unless it is so deviant that it doesn't
#appear on the graph. A p-value for the
#permutation is returned (see above function
#description; also see above for explanation
#of "allow.ties" bool)
#-------------------------------------
perform.full.perm.test<-function(data,n.perms=10000,allow.ties=TRUE,
                                 print.status=FALSE,display.graphics=TRUE) {
  vec <- find.vectors(data)
  actual.slope <- find.slope(vec[[1]],vec[[2]])
  new.data <- data
  slopes <- numeric(n.perms)
  for(p in 1:n.perms){
    new.data$outcome <- sample(data$outcome)
    vec <- find.vectors(new.data)
    slopes[p] <- find.slope(vec[[1]],vec[[2]])
    if(print.status & (p%%100==0)) {
      print(paste("Permutation", p, "has a slope difference of", slopes[p]-actual.slope))
    }
  }
  if(display.graphics) {
    hist(slopes,xlim=c(2*min(slopes),-2*min(slopes)))
    points(actual.slope, 0, pch=19, cex=1, col="red")
  }
  return(find.p.value(actual.slope,slopes,n.perms,allow.ties))
}
#-------------------------------------

#This function works essentially by the same
#principle as the above function with one
#exception. Here the values of proportion
#of runs that crowdsource are permuted
#across the bins. The default is to do
#a random permutation, but (by setting
#with.replacement=TRUE) the user can
#assign each bin to a random proportion
#drawn from the total set of possibilities
#with replacement. Otherwise, this function
#works the same as the previous one.
#-------------------------------------
perform.bin.perm.test<-function(data,n.perms=10000,with.replacement=FALSE,
                                allow.ties=TRUE,print.status=FALSE,
                                display.graphics=TRUE) {
  vec <- find.vectors(data)
  b <- vec[[1]]
  p <- vec[[2]]
  actual.slope <- find.slope(b,p)
  slopes <- numeric(n.perms)
  for(i in 1:n.perms){
    p.new <- sample(p, length(p), replace=with.replacement)
    slopes[i] <- find.slope(b,p.new)
    if(print.status & (i%%100==0)) {
      print(paste("Permutation", i, "has a slope difference of", slopes[i]-actual.slope))
    }
  }
  if(display.graphics) {
    hist(slopes,xlim=c(2*min(slopes),-2*min(slopes)))
    points(actual.slope, 0, pch=19, cex=1, col="red")
  }
  return(find.p.value(actual.slope,slopes,n.perms,allow.ties))
}
#-------------------------------------

#-----------------
#RUNNING THE TESTS
#-----------------
#Set the working directory
setwd("/Users/ben/Downloads")
#Read in the data
data <- read.csv("step_3_out.csv",header=TRUE)
#Perform the permutation test (where the bins are permuted)
perform.bin.perm.test(data,n.perms=100000,with.replacement=FALSE,
                      allow.ties=TRUE,print.status=TRUE,display.graphics=TRUE)
#Perform a second permutation test (where all the outcomes are permuted)
perform.full.perm.test(data,n.perms=100000,allow.ties=TRUE,
                       print.status=TRUE,display.graphics=TRUE)
#-----------------
