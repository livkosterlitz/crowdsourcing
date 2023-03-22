# Packages-----
## Download and install the package
#install.packages("igraph")
## Load package
library(igraph)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(gridGraphics)
theme_set(theme_cowplot())
citation(package = "igraph")

# Colors ---------
red_Ec_color <- rgb(255, 48, 48, maxColorValue = 255)
blue_Kp_color <- rgb(0, 0, 238, maxColorValue = 255)
yellow_Se_color <- rgb(242, 184, 0, maxColorValue = 255)
purple_Ec.Kp_color <- rgb(178, 58, 238, maxColorValue = 255)
orange_Ec.Se_color <- rgb(238, 118, 0, maxColorValue = 255)
green_Se.Kp_color <- rgb(70, 178, 0, maxColorValue = 255)
pie_color <- rgb(169, 169, 169, maxColorValue = 255)
pie_arrow <- rgb(226, 210, 195, maxColorValue = 255)

# Functions -----
draw.pie<-function(center.x,
                   center.y,
                   radius,
                   edge.color,
                   edge.width,
                   start.angle,
                   N.slices,
                   label.size,
                   label.pos,
                   colors,
                   labels,
                   text.colors,
                   x.y.ratio=1) {
  for(s in 1:N.slices) {
    angles <- seq(start.angle + (s-1)*(2*pi/N.slices),
                  start.angle + s*(2*pi/N.slices),
                  0.01)
    x.values <- center.x + radius*cos(angles)
    y.values <- center.y + radius*sin(angles)*x.y.ratio
    x.values <- c(center.x, x.values)
    y.values <- c(center.y, y.values)
    polygon(x.values,y.values,border=edge.color,col=colors[N.slices-s+1],lwd=edge.width)
    text(center.x + label.pos*radius*cos(start.angle + (s-.5)*(2*pi/N.slices)),
         center.y + label.pos*radius*sin(start.angle + (s-.5)*(2*pi/N.slices)),
         labels[N.slices-s+1],cex=label.size,col=text.colors[N.slices-s+1],
         adj=c(0.5,0.5),font=2)
  }
}

draw.pie.yaxis<-function(center.x,center.y,radius,edge.color,edge.width,start.angle,N.slices,label.size,label.pos,colors,labels,text.colors,axis.multiplier,xy.ratio) {
  for(s in 1:N.slices) {
    angles <- seq(start.angle + (s-1)*(2*pi/N.slices),
                  start.angle + s*(2*pi/N.slices),
                  0.01)
    x.values <- center.x + radius*cos(angles)
    y.values <- center.y + radius*sin(angles)*axis.multiplier*xy.ratio
    x.values <- c(center.x, x.values)
    y.values <- c(center.y, y.values)
    polygon(x.values,y.values,border=edge.color,col=colors[N.slices-s+1],lwd=edge.width)
    text(center.x + label.pos*radius*cos(start.angle + (s-.5)*(2*pi/N.slices)),
         center.y + label.pos*radius*sin(start.angle + (s-.5)*(2*pi/N.slices)),
         labels[N.slices-s+1],cex=label.size,col=text.colors[N.slices-s+1],
         adj=c(0.5,0.5),font=2)
  }
}

draw.arrow<-function(x.start,
                     y.start,
                     x.end,
                     y.end,
                     prop.drawn,
                     displacement,
                     color,
                     length,
                     angle,
                     width,
                     type=1) {
  if(x.start < x.end) {
    d <- sqrt(((y.end-y.start)^2) + ((x.end-x.start)^2))
    m <- (y.end-y.start)/(x.end-x.start)
    delta.x <- (d*0.5*(1 - prop.drawn))/sqrt(1+(m^2))
    new.x.start <- x.start + delta.x
    new.y.start <- y.start + m * delta.x
    delta.x.2 <- (d*prop.drawn)/sqrt(1+(m^2)) 
    new.x.end <- new.x.start + delta.x.2
    new.y.end <- new.y.start + m * delta.x.2
    if(displacement >= 0) {
      if(y.start < y.end) {
        delta.x.3 <- displacement/sqrt(1+(1/m)^2)
        delta.y.3 <- abs((-1/m)*delta.x.3)
        new.x.start <- new.x.start - delta.x.3
        new.y.start <- new.y.start + delta.y.3
        new.x.end <- new.x.end - delta.x.3
        new.y.end <- new.y.end + delta.y.3 
      } else {
        if(y.start > y.end) {
          delta.x.3 <- displacement/sqrt(1+(1/m)^2)
          delta.y.3 <- abs((-1/m)*delta.x.3)
          new.x.start <- new.x.start + delta.x.3
          new.y.start <- new.y.start + delta.y.3
          new.x.end <- new.x.end + delta.x.3
          new.y.end <- new.y.end + delta.y.3 
        } else {
          new.y.start <- new.y.start + displacement
          new.y.end <- new.y.end + displacement 
        }
      }
    } else {
      if(y.start < y.end) {
        delta.x.3 <- -displacement/sqrt(1+(1/m)^2)
        delta.y.3 <- abs((-1/m)*delta.x.3)
        new.x.start <- new.x.start + delta.x.3
        new.y.start <- new.y.start - delta.y.3
        new.x.end <- new.x.end + delta.x.3
        new.y.end <- new.y.end - delta.y.3 
      } else {
        if(y.start > y.end) {
          delta.x.3 <- -displacement/sqrt(1+(1/m)^2)
          delta.y.3 <- abs((-1/m)*delta.x.3)
          new.x.start <- new.x.start - delta.x.3
          new.y.start <- new.y.start - delta.y.3
          new.x.end <- new.x.end - delta.x.3
          new.y.end <- new.y.end - delta.y.3 
        } else {
          new.y.start <- new.y.start + displacement
          new.y.end <- new.y.end + displacement 
        }
      }
    }
    arrows(new.x.start,new.y.start,
           new.x.end,new.y.end,
           length, angle,
           col=color,
           lwd=width,
           lty=type)
  } else {
    if(x.start > x.end) {
      x.s <- x.end
      y.s <- y.end
      x.e <- x.start
      y.e <- y.start
      d <- sqrt(((y.e-y.s)^2) + ((x.e-x.s)^2))
      m <- (y.e-y.s)/(x.e-x.s)
      delta.x <- (d*0.5*(1 - prop.drawn))/sqrt(1+(m^2))
      new.x.start <- x.s + delta.x
      new.y.start <- y.s + m * delta.x
      delta.x.2 <- (d*prop.drawn)/sqrt(1+(m^2)) 
      new.x.end <- new.x.start + delta.x.2
      new.y.end <- new.y.start + m * delta.x.2
      if(displacement >= 0) {
        if(y.s < y.e) {
          delta.x.3 <- displacement/sqrt(1+(1/m)^2)
          delta.y.3 <- abs((-1/m)*delta.x.3)
          new.x.start <- new.x.start - delta.x.3
          new.y.start <- new.y.start + delta.y.3
          new.x.end <- new.x.end - delta.x.3
          new.y.end <- new.y.end + delta.y.3 
        } else {
          if(y.s > y.e) {
            delta.x.3 <- displacement/sqrt(1+(1/m)^2)
            delta.y.3 <- abs((-1/m)*delta.x.3)
            new.x.start <- new.x.start + delta.x.3
            new.y.start <- new.y.start + delta.y.3
            new.x.end <- new.x.end + delta.x.3
            new.y.end <- new.y.end + delta.y.3 
          } else {
            new.y.start <- new.y.start + displacement
            new.y.end <- new.y.end + displacement 
          }
        }
      } else {
        if(y.s < y.e) {
          delta.x.3 <- -displacement/sqrt(1+(1/m)^2)
          delta.y.3 <- abs((-1/m)*delta.x.3)
          new.x.start <- new.x.start + delta.x.3
          new.y.start <- new.y.start - delta.y.3
          new.x.end <- new.x.end + delta.x.3
          new.y.end <- new.y.end - delta.y.3 
        } else {
          if(y.s > y.e) {
            delta.x.3 <- -displacement/sqrt(1+(1/m)^2)
            delta.y.3 <- abs((-1/m)*delta.x.3)
            new.x.start <- new.x.start - delta.x.3
            new.y.start <- new.y.start - delta.y.3
            new.x.end <- new.x.end - delta.x.3
            new.y.end <- new.y.end - delta.y.3 
          } else {
            new.y.start <- new.y.start + displacement
            new.y.end <- new.y.end + displacement 
          }
        }
      }
      arrows(new.x.end,new.y.end,
             new.x.start,new.y.start,
             length, angle,
             col=color,
             lwd=width,
             lty=type)
    } else {
      if(y.start < y.end) {
        d <- sqrt(((y.end-y.start)^2) + ((x.end-x.start)^2))
        delta.y <- (d*0.5*(1 - prop.drawn))
        new.x.start <- x.start
        new.y.start <- y.start + delta.y
        new.x.end <- x.end
        new.y.end <- y.end - delta.y
        if(displacement >=0) {
          new.x.start <- new.x.start + displacement
          new.x.end <- new.x.end + displacement
        } else {
          new.x.start <- new.x.start + displacement
          new.x.end <- new.x.end + displacement
        }
        arrows(new.x.start,new.y.start,
               new.x.end,new.y.end,
               length, angle,
               col=color,
               lwd=width,
               lty=type)
      } else {
        d <- sqrt(((y.end-y.start)^2) + ((x.end-x.start)^2))
        delta.y <- (d*0.5*(1 - prop.drawn))
        new.x.start <- x.start
        new.y.start <- y.start - delta.y
        new.x.end <- x.end
        new.y.end <- y.end + delta.y
        if(displacement >=0) {
          new.x.start <- new.x.start + displacement
          new.x.end <- new.x.end + displacement
        } else {
          new.x.start <- new.x.start + displacement
          new.x.end <- new.x.end + displacement
        }
        arrows(new.x.start,new.y.start,
               new.x.end,new.y.end,
               length, angle,
               col=color,
               lwd=width,
               lty=type)
      }
    }
  }
}
# create draw.arrow function for the network where we specify the shaft using  
# segment() and segment for the arrows with a solid line type

draw.circle.asp<-function(center.x,
                          center.y,
                          radius,
                          fill.color,
                          edge.color,
                          edge.width) {
  angles <- seq(0,2*pi,0.01)
  x.values <- center.x + radius*cos(angles)
  y.values <- center.y + radius*sin(angles)
  polygon(x.values,y.values,border=edge.color,col=fill.color,lwd=edge.width)
}

is.apart <- function(df,N.loci,g1,g2,m) {
  N.differences<-0
  for(p in 1:N.loci) {
    if(df[g1,1+p]!=df[g2,1+p]) {
      N.differences <- N.differences + 1
    }
  }
  if(N.differences==m) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.distance <- function(df,N.loci,g,m) {
  N.mutations<-0
  for(p in 1:N.loci) {
    if(df[g,1+p]==1) {
      N.mutations <- N.mutations + 1
    }
  }
  if(N.mutations==m) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

draw.landscape.graph<-function(df,
                               N.loci,
                               starting.position.labels,
                               ending.position.labels,
                               N.fit.estimates,
                               buffer.x,
                               buffers.y,
                               node.radius,
                               node.edge.width,
                               node.label.size,
                               label.pos,
                               node.edge.color="gray",
                               starting.color="white",
                               ending.color="black",
                               starting.text.color="black",
                               ending.text.color="white",
                               arrow.host.colors,
                               arrow.host.displacements,
                               arrow.prop.drawn,
                               arrow.width,
                               arrow.head,
                               gamma=1) {
  N.hosts<-length(unique(as.vector(df$hosts)))
  hosts<-unique(as.vector(df$hosts))
  df$x <- 0
  df$y <- 0
  for(h in hosts) {
    for(m in 0:N.loci) {
      if(m<=N.loci/2) {
        xpos <- buffer.x + ((m/N.loci)^gamma)*(1-2*buffer.x)
      } else {
        xpos <- 1-buffer.x - (((N.loci-m)/N.loci)^gamma)*(1-2*buffer.x)
      }
      count <- 1
      for(g in 1:length(df$hosts)) {
        if(df[g,1]==h & is.distance(df,N.loci,g,m)) {
          df[g,]$x <- xpos
          if(m==0 | m==N.loci) {
            df[g,]$y <- 0.5
          } else {
            df[g,]$y <- 1-buffers.y[m]-(1-2*buffers.y[m])*((count-1)/(choose(N.loci,m)-1)) 
          }
          count <- count + 1
        }
      }
    }  
  }
  
  for(m in 0:(N.loci-1)) {
    for(g in 1:length(df[df$hosts==hosts[1],]$hosts)) {
      for(h in hosts) {
        for(i in 1:length(df$hosts)) {
          if(df[i,1]==h & is.apart(df,N.loci,g,i,0)) {
            focal.genotype <- i
            break
          }
        }
        neighbor.genotypes<-numeric(0)
        for(i in 1:length(df$hosts)) {
          if(df[i,1]==h & is.apart(df,N.loci,focal.genotype,i,1)) {
            neighbor.genotypes<-c(neighbor.genotypes,i)
          }
        }
        
        for(n in neighbor.genotypes) {
          arrow.type <- determine.direction(df,N.loci,N.fit.estimates,focal.genotype,n)
          if(arrow.type == "forward") {
            xs <- df[focal.genotype,]$x
            ys <- df[focal.genotype,]$y
            xe <- df[n,]$x
            ye <- df[n,]$y
            head.length <- arrow.head
            arrow.type <- 1
          } else {
            if(arrow.type == "backward") {
              xe <- df[focal.genotype,]$x
              ye <- df[focal.genotype,]$y
              xs <- df[n,]$x
              ys <- df[n,]$y
              head.length <- arrow.head
              arrow.type <- 1
            } else {
              xs <- df[focal.genotype,]$x
              ys <- df[focal.genotype,]$y
              xe <- df[n,]$x
              ye <- df[n,]$y
              head.length <- 0
              arrow.type <- 2
            }
          }
          
          p<-arrow.prop.drawn
          delx <- (1-p)/2
          lenx <- abs(xs-xe)
          leny <- abs(ys-ye)
          slope <- leny / lenx
          len <- sqrt(1 + slope^2)
          #dely <- slope*delx
          #pp <- 1 - 2*(sqrt(delx^2 + dely^2))
          if(slope==0) {
            pp<-p
          } else {
            pp <- 1 - (2*delx)/len
          }
          draw.arrow(x.start=xs,
                     y.start=ys,
                     x.end=xe,
                     y.end=ye,
                     prop.drawn=pp,
                     displacement=arrow.host.displacements[which(hosts==h)],
                     color=arrow.host.colors[which(hosts==h)],
                     length=head.length,
                     angle=30,
                     width = arrow.width,
                     type=arrow.type)
        }
      }
    }
  }
  
  new.df<-df[df$hosts==hosts[1],]
  for(g in 1:length(new.df$hosts)) {
    cols<-character(0)
    t.cols<-character(0)
    labs<-character(0)
    for(p in 1:N.loci) {
      if(new.df[g,1+p]==0) {
        cols<-c(cols,starting.color)
        t.cols<-c(t.cols,starting.text.color)
        labs<-c(labs,starting.position.labels[p])
      } else {
        cols<-c(cols,ending.color)
        t.cols<-c(t.cols,ending.text.color)
        labs<-c(labs,ending.position.labels[p])
      }
    }
    draw.pie(center.x = new.df[g,]$x,
             center.y = new.df[g,]$y,
             radius=node.radius,
             edge.color=node.edge.color,
             edge.width=node.edge.width,
             start.angle=pi/2+pi/N.loci,
             N.slices=N.loci,
             label.size=node.label.size,
             label.pos=label.pos,
             colors=cols,
             labels=labs,
             text.colors=t.cols)
  }
}

generate.mutation.effects <- function(df,N.loci,N.fit.estimates,used.hosts) {
  h1.effects <- double(0)
  h2.effects <- double(0)
  h1.df <- df[df$hosts==used.hosts[1],]
  h2.df <- df[df$hosts==used.hosts[2],]
  for(i in 1:(length(h1.df$hosts)-1)) {
    for(j in (i+1):length(h1.df$hosts)) {
      if(is.apart(h1.df, N.loci, i, j, 1) & is.apart(h2.df, N.loci, i, j, 1)) {
        h1.first.node.fit <- mean(as.double(h1.df[i,(2+N.loci):(1+N.loci+N.fit.estimates)]))
        h1.second.node.fit <- mean(as.double(h1.df[j,(2+N.loci):(1+N.loci+N.fit.estimates)]))
        h2.first.node.fit <- mean(as.double(h2.df[i,(2+N.loci):(1+N.loci+N.fit.estimates)]))
        h2.second.node.fit <- mean(as.double(h2.df[j,(2+N.loci):(1+N.loci+N.fit.estimates)]))
        h1.effects <- c(h1.effects,log2(h1.second.node.fit)-log2(h1.first.node.fit))
        h2.effects <- c(h2.effects,log2(h2.second.node.fit)-log2(h2.first.node.fit))
      }
    }
  }
  return(list(h1.effects,h2.effects))
}

generate.network.list <- function(df,N.loci,N.fit.estimates,used.hosts,
                                  wt.alleles, mut.alleles) {
  network.list <- character(0)
  for(h in used.hosts) {
    new.df <- df[df$hosts==h,]
    for(i in 1:(length(new.df$hosts)-1)) {
      for(j in (i+1):length(new.df$hosts)) {
        if(is.apart(new.df, N.loci, i, j, 1)) {
          first.node.fit <- mean(as.double(new.df[i,(2+N.loci):(1+N.loci+N.fit.estimates)]))
          second.node.fit <- mean(as.double(new.df[j,(2+N.loci):(1+N.loci+N.fit.estimates)]))
          if(first.node.fit != second.node.fit) {
            first.node.mutations <- as.integer(new.df[i,2:(1+N.loci)])
            second.node.mutations <- as.integer(new.df[j,2:(1+N.loci)])
            if(first.node.fit < second.node.fit) {
              network.list <- c(network.list,
                                name.node(first.node.mutations,wt.alleles, mut.alleles),
                                name.node(second.node.mutations,wt.alleles, mut.alleles))
            } else {
              network.list <- c(network.list,
                                name.node(second.node.mutations,wt.alleles, mut.alleles),
                                name.node(first.node.mutations,wt.alleles, mut.alleles))
            }
          }
        }
      }
    }
  }
  return(network.list)
}

name.node<-function(mutations,wt.alleles,mut.alleles) {
  alleles<-character(0)
  for(i in 1:length(mutations)) {
    if(mutations[i]==0) {
      alleles<-c(alleles,wt.alleles[i])
    } else {
      alleles<-c(alleles,mut.alleles[i])
    }
  }
  return(paste(alleles,collapse=""))
}

FindCycles = function(g) {
  Cycles = NULL
  for(v1 in V(g)) {
    if(degree(g, v1, mode="in") == 0) { next }
    GoodNeighbors = neighbors(g, v1, mode="out")
    GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
    for(v2 in GoodNeighbors) {
      TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
      TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
      TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles  = c(Cycles, TempCyc)
    }
  }
  Cycles
}

draw.landscape.merge.graph.se<-function(df,
                                        N.loci,
                                        starting.position.labels,
                                        ending.position.labels,
                                        N.fit.estimates,
                                        buffer.x,
                                        buffers.y,
                                        node.radius,
                                        node.edge.width,
                                        node.label.size,
                                        label.pos,
                                        node.edge.color="gray",
                                        starting.color="white",
                                        ending.color="black",
                                        starting.text.color="black",
                                        ending.text.color="white",
                                        merge_arrow_color="grey",
                                        arrow.host.colors,
                                        arrow.host.displacements,
                                        arrow.prop.drawn,
                                        arrow.width,
                                        arrow.head,
                                        gamma=1,
                                        x.y.ratio=1) {
  N.hosts<-length(unique(as.vector(df$hosts)))
  hosts<-unique(as.vector(df$hosts))
  df$x <- 0
  df$y <- 0
  ### This loop assigns and x and y value position for each genotype
  for(h in hosts) { # loop through species (e.g., Ec, Kp, Se)
    for(m in 0:N.loci) { # loop through the number of loci (e.g., 1, 2, 3, 4, 5)
      if(m<=N.loci/2) { 
        xpos <- buffer.x + ((m/N.loci)^gamma)*(1-2*buffer.x) 
      } else {
        xpos <- 1-buffer.x - (((N.loci-m)/N.loci)^gamma)*(1-2*buffer.x)
      }
      count <- 1
      for(g in 1:length(df$hosts)) { # loop through all rows in the dataframe
        if(df[g,1]==h & is.distance(df,N.loci,g,m)) { # if the row is the right species and the rows genotype has the right number of mutations
          df[g,]$x <- xpos
          if(m==0 | m==N.loci) { # if the genotype is the ancestor (i.e., no mutations) or the genotype is the evolved (i.e., all of the mutations)
            df[g,]$y <- 0.5 # put the ancestor and evolved genotypes in the middle of the y-axis
          } else {
            df[g,]$y <- 1-buffers.y[m]-(1-2*buffers.y[m])*((count-1)/(choose(N.loci,m)-1)) 
          }
          count <- count + 1
        }
      }
    }  
  }
  for(g in 1:length(df[df$hosts==hosts[1],]$hosts)) { # loop through the number of genotypes (1 to 32)
    for(h in hosts) { # loop through hosts (E.coli, Kleb, Se)
      for(i in 1:length(df$hosts)) { # loop through rows (1 to 96)
        if(df[i,1]==h & is.apart(df,N.loci,g,i,0)) { # if the row is for the particular host and matches the genotype 'g' then this is the focal genotype
          focal.genotype <- i
          break
        }
      }
      # g: 1 through 32
      # h: Ec, Kp, Se
      # i: the line number for g in the data frame
      neighbor.genotypes<-numeric(0)
      for(i in 1:length(df$hosts)) { # loop through the rows 
        if(df[i,1]==h & is.apart(df,N.loci,focal.genotype,i,1)) { # if the host matches and the genotype differs by one mutation then store the row number for the neighbor genotype 
          neighbor.genotypes<-c(neighbor.genotypes,i)
        }
      }
      merge.neighbor.genotypes <- vector(mode = "list", length = length(hosts))
      names(merge.neighbor.genotypes) <- hosts
      for(i in 1:length(df$hosts)) {  
        if(is.apart(df,N.loci,focal.genotype,i,1)) {
          merge.neighbor.genotypes[[df[i,1]]] <- append(merge.neighbor.genotypes[[df[i,1]]], i)
        }
      }
      merge.focal.genotypes <- vector(mode = "list", length = length(hosts))
      names(merge.focal.genotypes) <- hosts
      for(i in 1:length(df$hosts)) {  
        if(is.apart(df,N.loci,focal.genotype,i,0)) {
          merge.focal.genotypes[[df[i,1]]] <- append(merge.focal.genotypes[[df[i,1]]], i)
        }
      }
      for(n in neighbor.genotypes) { # loop through the neighbor genotypes
        if (df$Num_of_mutations[focal.genotype] < df$Num_of_mutations[n]) {
          arrow.type <- determine.direction.se(df,N.loci,N.fit.estimates,focal.genotype,n)
          index = which(neighbor.genotypes == n)
          merge.check <- vector()
          for(s in hosts){
            merge.check <- c(merge.check,determine.direction.se(df,N.loci,N.fit.estimates,merge.focal.genotypes[[s]],merge.neighbor.genotypes[[s]][index]))
          }
          if(arrow.type == "forward") {
            xs <- df[focal.genotype,]$x
            ys <- df[focal.genotype,]$y
            xe <- df[n,]$x
            ye <- df[n,]$y
            head.length <- arrow.head
            arrow.type.num <- 1
          } else {
            if(arrow.type == "backward") {
              xe <- df[focal.genotype,]$x
              ye <- df[focal.genotype,]$y
              xs <- df[n,]$x
              ys <- df[n,]$y
              head.length <- arrow.head
              arrow.type.num <- 5
            } else {
              xs <- df[focal.genotype,]$x
              ys <- df[focal.genotype,]$y
              xe <- df[n,]$x
              ye <- df[n,]$y
              head.length <- 0
              arrow.type.num <- 3
            }
          }
          p<-arrow.prop.drawn
          delx <- (1-p)/2
          lenx <- abs(xs-xe)
          leny <- abs(ys-ye)
          slope <- leny / lenx
          len <- sqrt(1 + slope^2)
          if(slope==0) {
            pp<-p
          } else {
            pp <- 1 - (2*delx)/len
          }
          if(length(unique(merge.check)) == 1) {
            draw.arrow(x.start=xs,
                       y.start=ys,
                       x.end=xe,
                       y.end=ye,
                       prop.drawn=pp,
                       displacement=0,
                       color=merge_arrow_color,
                       length=head.length,
                       angle=30,
                       width = arrow.width,
                       type=arrow.type.num)
          }   
        }
      }
    }
  }
  for(g in 1:length(df[df$hosts==hosts[1],]$hosts)) { # loop through the number of genotypes (1 to 32)
    for(h in hosts) { # loop through hosts (E.coli, Kleb, Se)
      for(i in 1:length(df$hosts)) { # loop through rows (1 to 96)
        if(df[i,1]==h & is.apart(df,N.loci,g,i,0)) { # if the row is for the particular host and matches the genotype 'g' then this is the focal genotype
          focal.genotype <- i
          break
        }
      }
      # g: 1 through 32
      # h: Ec, Kp, Se
      # i: the line number for g in the data frame
      neighbor.genotypes<-numeric(0)
      for(i in 1:length(df$hosts)) { # loop through the rows 
        if(df[i,1]==h & is.apart(df,N.loci,focal.genotype,i,1)) { # if the host matches and the genotype differs by one mutation then store the row number for the neighbor genotype 
          neighbor.genotypes<-c(neighbor.genotypes,i)
        }
      }
      merge.neighbor.genotypes <- vector(mode = "list", length = length(hosts))
      names(merge.neighbor.genotypes) <- hosts
      for(i in 1:length(df$hosts)) {  
        if(is.apart(df,N.loci,focal.genotype,i,1)) {
          merge.neighbor.genotypes[[df[i,1]]] <- append(merge.neighbor.genotypes[[df[i,1]]], i)
        }
      }
      merge.focal.genotypes <- vector(mode = "list", length = length(hosts))
      names(merge.focal.genotypes) <- hosts
      for(i in 1:length(df$hosts)) {  
        if(is.apart(df,N.loci,focal.genotype,i,0)) {
          merge.focal.genotypes[[df[i,1]]] <- append(merge.focal.genotypes[[df[i,1]]], i)
        }
      }
      for(n in neighbor.genotypes) { # loop through the neighbor genotypes
        if (df$Num_of_mutations[focal.genotype] < df$Num_of_mutations[n]) {
          arrow.type <- determine.direction.se(df,N.loci,N.fit.estimates,focal.genotype,n)
          index = which(neighbor.genotypes == n)
          merge.check <- vector()
          for(s in hosts){
            merge.check <- c(merge.check,determine.direction.se(df,N.loci,N.fit.estimates,merge.focal.genotypes[[s]],merge.neighbor.genotypes[[s]][index]))
          }
          if(arrow.type == "forward") {
            xs <- df[focal.genotype,]$x
            ys <- df[focal.genotype,]$y
            xe <- df[n,]$x
            ye <- df[n,]$y
            head.length <- arrow.head
            arrow.type.num <- 1
          } else {
            if(arrow.type == "backward") {
              xe <- df[focal.genotype,]$x
              ye <- df[focal.genotype,]$y
              xs <- df[n,]$x
              ys <- df[n,]$y
              head.length <- arrow.head
              arrow.type.num <- 5
            } else {
              xs <- df[focal.genotype,]$x
              ys <- df[focal.genotype,]$y
              xe <- df[n,]$x
              ye <- df[n,]$y
              head.length <- 0
              arrow.type.num <- 3
            }
          }
          p<-arrow.prop.drawn
          delx <- (1-p)/2
          lenx <- abs(xs-xe)
          leny <- abs(ys-ye)
          slope <- leny / lenx
          len <- sqrt(1 + slope^2)
          if(slope==0) {
            pp<-p
          } else {
            pp <- 1 - (2*delx)/len
          }
          if(length(unique(merge.check)) != 1) {
            draw.arrow(x.start=xs,
                       y.start=ys,
                       x.end=xe,
                       y.end=ye,
                       prop.drawn=pp,
                       displacement=arrow.host.displacements[which(hosts==h)],
                       color=arrow.host.colors[which(hosts==h)],
                       length=head.length,
                       angle=30,
                       width = arrow.width,
                       type=arrow.type.num)
          }  
        }
      }
    }
  }
  new.df<-df[df$hosts==hosts[1],]
  for(g in 1:length(new.df$hosts)) {
    cols<-character(0)
    t.cols<-character(0)
    labs<-character(0)
    for(p in 1:N.loci) {
      if(new.df[g,1+p]==0) {
        cols<-c(cols,starting.color)
        t.cols<-c(t.cols,starting.text.color)
        labs<-c(labs,starting.position.labels[p])
      } else {
        cols<-c(cols,ending.color)
        t.cols<-c(t.cols,ending.text.color)
        labs<-c(labs,ending.position.labels[p])
      }
    }
    draw.pie(center.x = new.df[g,]$x,
             center.y = new.df[g,]$y,
             radius=node.radius,
             edge.color=node.edge.color,
             edge.width=node.edge.width,
             start.angle=pi/2+pi/N.loci,
             N.slices=N.loci,
             label.size=node.label.size,
             label.pos=label.pos,
             colors=cols,
             labels=labs,
             text.colors=t.cols,
             x.y.ratio=x.y.ratio)
  }
}

determine.direction <- function(df,N.loci,N.fit.estimates,focal.genotype,n) {
  fit.origin <- 0.0
  fit.terminus <- 0.0
  for(e in 1:N.fit.estimates) {
    fit.origin<-fit.origin + as.numeric(df[focal.genotype,1+N.loci+e])
    fit.terminus<-fit.terminus + as.numeric(df[n,1+N.loci+e])
  }
  fit.origin <- fit.origin/N.fit.estimates
  fit.terminus <- fit.terminus/N.fit.estimates
  if(fit.origin < fit.terminus) {
    return("forward")
  } else {
    if(fit.origin > fit.terminus) {
      return("backward")
    } else {
      return("neutral")
    }
  }
}

determine.direction.se <- function(df,N.loci,N.fit.estimates,focal.genotype,n) {
  fit.origin <- 0.0
  fit.terminus <- 0.0
  for(e in 1:N.fit.estimates) {
    fit.origin<-fit.origin + as.numeric(df[focal.genotype,1+N.loci+e])
    fit.terminus<-fit.terminus + as.numeric(df[n,1+N.loci+e])
  }
  fit.origin <- fit.origin/N.fit.estimates
  fit.terminus <- fit.terminus/N.fit.estimates
  fit.origin.se <- df$MIC_se[focal.genotype]
  fit.terminus.se <- df$MIC_se[n]
  if(fit.origin < fit.terminus & ((fit.origin+fit.origin.se)<(fit.terminus-fit.terminus.se))) {
    return("forward")
  } else {
    if(fit.origin > fit.terminus & ((fit.origin-fit.origin.se)>(fit.terminus+fit.terminus.se))) {
      return("backward")
    } else {
      return("neutral")
    }
  }
}

# TEM landscape key --------------
draw.figure.TEM.ancestor.key<-function(){
  df<-read.csv("inputs/ancestor.csv",header=TRUE)
  plot.new()
  draw.landscape.graph(df,
                       N.loci=5,
                       starting.position.labels=c("g","A","E","M","G"),
                       ending.position.labels=c("a","G","K","T","S"),
                       N.fit.estimates=1,
                       buffer.x=0.5,
                       buffers.y=c(0,0,0,0),
                       node.radius=0.4,
                       node.edge.width=1,
                       node.label.size=1,
                       label.pos=.66,
                       node.edge.color='black',
                       starting.color="white",
                       ending.color=pie_color,
                       starting.text.color=pie_color,
                       ending.text.color="black",
                       arrow.host.colors=c("black"),
                       arrow.host.displacements=c(0),
                       arrow.prop.drawn=0.6,
                       arrow.width=1,
                       arrow.head=.2,
                       gamma=1)
}
draw.figure.TEM.evolved.key<-function(){
  df<-read.csv("inputs/evolved.csv",header=TRUE)
  plot.new()
  draw.landscape.graph(df,
                       N.loci=5,
                       starting.position.labels=c("g","A","E","M","G"),
                       ending.position.labels=c("a","G","K","T","S"),
                       N.fit.estimates=1,
                       buffer.x=0.5,
                       buffers.y=c(0,0,0,0),
                       node.radius=0.4,
                       node.edge.width=1,
                       node.label.size=1,
                       label.pos=.66,
                       node.edge.color='black',
                       starting.color="white",
                       ending.color=pie_color,
                       starting.text.color=pie_color,
                       ending.text.color="black",
                       arrow.host.colors=c("black"),
                       arrow.host.displacements=c(0),
                       arrow.prop.drawn=0.6,
                       arrow.width=1,
                       arrow.head=.2,
                       gamma=1)
}

# TEM landscape network graph --------------
draw.figure.TEM.network <-function(network.x.y.ratio){
  df<-read.csv("inputs/resistance_levels.csv",header=TRUE)
  plot.new()
  draw.landscape.merge.graph.se(df,
                                N.loci=5,
                                starting.position.labels=c("","","","",""),
                                ending.position.labels=c("","","","",""),
                                N.fit.estimates=1,
                                buffer.x=0,
                                buffers.y=c(0.2,0,0,0.2),
                                node.radius=0.025,
                                node.edge.width=0.75,
                                node.label.size=0,
                                label.pos=.65,
                                node.edge.color='black',
                                starting.color="white",
                                ending.color=pie_color,
                                starting.text.color=pie_color,
                                ending.text.color="black",
                                merge_arrow_color=pie_arrow,
                                arrow.host.colors=c(red_Ec_color, blue_Kp_color, yellow_Se_color),
                                arrow.host.displacements=c(-0.012,0,0.012),
                                arrow.prop.drawn=0.6,
                                arrow.width=1.75,
                                arrow.head=.08, 
                                gamma=1,
                                x.y.ratio=network.x.y.ratio)}

draw.landscape.blowups <-function(df,
                                  N.loci,
                                  starting.position.labels,
                                  ending.position.labels,
                                  starting.blowup.labels,
                                  ending.blowup.labels,
                                  N.fit.estimates,
                                  buffer.x,
                                  buffers.y,
                                  node.radius,
                                  node.edge.width,
                                  blowup.node.radius,
                                  blowup.distance,
                                  node.label.size,
                                  label.pos,
                                  node.edge.color="gray",
                                  starting.color="white",
                                  ending.color="black",
                                  starting.text.color="black",
                                  ending.text.color="white",
                                  merge_arrow_color="grey",
                                  arrow.host.colors,
                                  arrow.host.displacements,
                                  arrow.prop.drawn,
                                  arrow.width,
                                  arrow.head,
                                  gamma=1,
                                  x.y.ratio=1) {
  N.hosts<-length(unique(as.vector(df$hosts)))
  hosts<-unique(as.vector(df$hosts))
  df$x <- 0
  df$y <- 0
  ### This loop assigns and x and y value position for each genotype
  for(h in hosts) { # loop through species (e.g., Ec, Kp, Se)
    for(m in 0:N.loci) { # loop through the number of loci (e.g., 1, 2, 3, 4, 5)
      if(m<=N.loci/2) { 
        xpos <- buffer.x + ((m/N.loci)^gamma)*(1-2*buffer.x) 
      } else {
        xpos <- 1-buffer.x - (((N.loci-m)/N.loci)^gamma)*(1-2*buffer.x)
      }
      count <- 1
      for(g in 1:length(df$hosts)) { # loop through all rows in the dataframe
        if(df[g,1]==h & is.distance(df,N.loci,g,m)) { # if the row is the right species and the rows genotype has the right number of mutations
          df[g,]$x <- xpos
          if(m==0 | m==N.loci) { # if the genotype is the ancestor (i.e., no mutations) or the genotype is the evolved (i.e., all of the mutations)
            df[g,]$y <- 0.5 # put the ancestor and evolved genotypes in the middle of the y-axis
          } else {
            df[g,]$y <- 1-buffers.y[m]-(1-2*buffers.y[m])*((count-1)/(choose(N.loci,m)-1)) 
          }
          count <- count + 1
        }
      }
    }  
  }
  for(g in 1:length(df[df$hosts==hosts[1],]$hosts)) { # loop through the number of genotypes (1 to 32)
    for(h in hosts) { # loop through hosts (E.coli, Kleb, Se)
      for(i in 1:length(df$hosts)) { # loop through rows (1 to 96)
        if(df[i,1]==h & is.apart(df,N.loci,g,i,0)) { # if the row is for the particular host and matches the genotype 'g' then this is the focal genotype
          focal.genotype <- i
          break
        }
      }
      # g: 1 through 32
      # h: Ec, Kp, Se
      # i: the line number for g in the data frame
      neighbor.genotypes<-numeric(0)
      for(i in 1:length(df$hosts)) { # loop through the rows 
        if(df[i,1]==h & is.apart(df,N.loci,focal.genotype,i,1)) { # if the host matches and the genotype differs by one mutation then store the row number for the neighbor genotype 
          neighbor.genotypes<-c(neighbor.genotypes,i)
        }
      }
      merge.neighbor.genotypes <- vector(mode = "list", length = length(hosts))
      names(merge.neighbor.genotypes) <- hosts
      for(i in 1:length(df$hosts)) {  
        if(is.apart(df,N.loci,focal.genotype,i,1)) {
          merge.neighbor.genotypes[[df[i,1]]] <- append(merge.neighbor.genotypes[[df[i,1]]], i)
        }
      }
      merge.focal.genotypes <- vector(mode = "list", length = length(hosts))
      names(merge.focal.genotypes) <- hosts
      for(i in 1:length(df$hosts)) {  
        if(is.apart(df,N.loci,focal.genotype,i,0)) {
          merge.focal.genotypes[[df[i,1]]] <- append(merge.focal.genotypes[[df[i,1]]], i)
        }
      }
      for(n in neighbor.genotypes) { # loop through the neighbor genotypes
        if (df$Num_of_mutations[focal.genotype] < df$Num_of_mutations[n]) {
          arrow.type <- determine.direction.se(df,N.loci,N.fit.estimates,focal.genotype,n)
          index = which(neighbor.genotypes == n)
          merge.check <- vector()
          for(s in hosts){
            merge.check <- c(merge.check,determine.direction.se(df,N.loci,N.fit.estimates,merge.focal.genotypes[[s]],merge.neighbor.genotypes[[s]][index]))
          }
          if(arrow.type == "forward") {
            xs <- df[focal.genotype,]$x
            ys <- df[focal.genotype,]$y
            xe <- df[n,]$x
            ye <- df[n,]$y
            head.length <- arrow.head
            arrow.type.num <- 1
          } else {
            if(arrow.type == "backward") {
              xe <- df[focal.genotype,]$x
              ye <- df[focal.genotype,]$y
              xs <- df[n,]$x
              ys <- df[n,]$y
              head.length <- arrow.head
              arrow.type.num <- 5
            } else {
              xs <- df[focal.genotype,]$x
              ys <- df[focal.genotype,]$y
              xe <- df[n,]$x
              ye <- df[n,]$y
              head.length <- 0
              arrow.type.num <- 3
            }
          }
          p<-arrow.prop.drawn
          delx <- (1-p)/2
          lenx <- abs(xs-xe)
          leny <- abs(ys-ye)
          slope <- leny / lenx
          len <- sqrt(1 + slope^2)
          if(slope==0) {
            pp<-p
          } else {
            pp <- 1 - (2*delx)/len
          }
          if(length(unique(merge.check)) == 1) {
            draw.arrow(x.start=xs,
                       y.start=ys,
                       x.end=xe,
                       y.end=ye,
                       prop.drawn=pp,
                       displacement=0,
                       color=merge_arrow_color,
                       length=head.length,
                       angle=30,
                       width = arrow.width,
                       type=arrow.type.num)
          }   
        }
      }
    }
  }
  for(g in 1:length(df[df$hosts==hosts[1],]$hosts)) { # loop through the number of genotypes (1 to 32)
    for(h in hosts) { # loop through hosts (E.coli, Kleb, Se)
      for(i in 1:length(df$hosts)) { # loop through rows (1 to 96)
        if(df[i,1]==h & is.apart(df,N.loci,g,i,0)) { # if the row is for the particular host and matches the genotype 'g' then this is the focal genotype
          focal.genotype <- i
          break
        }
      }
      # g: 1 through 32
      # h: Ec, Kp, Se
      # i: the line number for g in the data frame
      neighbor.genotypes<-numeric(0)
      for(i in 1:length(df$hosts)) { # loop through the rows 
        if(df[i,1]==h & is.apart(df,N.loci,focal.genotype,i,1)) { # if the host matches and the genotype differs by one mutation then store the row number for the neighbor genotype 
          neighbor.genotypes<-c(neighbor.genotypes,i)
        }
      }
      merge.neighbor.genotypes <- vector(mode = "list", length = length(hosts))
      names(merge.neighbor.genotypes) <- hosts
      for(i in 1:length(df$hosts)) {  
        if(is.apart(df,N.loci,focal.genotype,i,1)) {
          merge.neighbor.genotypes[[df[i,1]]] <- append(merge.neighbor.genotypes[[df[i,1]]], i)
        }
      }
      merge.focal.genotypes <- vector(mode = "list", length = length(hosts))
      names(merge.focal.genotypes) <- hosts
      for(i in 1:length(df$hosts)) {  
        if(is.apart(df,N.loci,focal.genotype,i,0)) {
          merge.focal.genotypes[[df[i,1]]] <- append(merge.focal.genotypes[[df[i,1]]], i)
        }
      }
      for(n in neighbor.genotypes) { # loop through the neighbor genotypes
        if (df$Num_of_mutations[focal.genotype] < df$Num_of_mutations[n]) {
          arrow.type <- determine.direction.se(df,N.loci,N.fit.estimates,focal.genotype,n)
          index = which(neighbor.genotypes == n)
          merge.check <- vector()
          for(s in hosts){
            merge.check <- c(merge.check,determine.direction.se(df,N.loci,N.fit.estimates,merge.focal.genotypes[[s]],merge.neighbor.genotypes[[s]][index]))
          }
          if(arrow.type == "forward") {
            xs <- df[focal.genotype,]$x
            ys <- df[focal.genotype,]$y
            xe <- df[n,]$x
            ye <- df[n,]$y
            head.length <- arrow.head
            arrow.type.num <- 1
          } else {
            if(arrow.type == "backward") {
              xe <- df[focal.genotype,]$x
              ye <- df[focal.genotype,]$y
              xs <- df[n,]$x
              ys <- df[n,]$y
              head.length <- arrow.head
              arrow.type.num <- 5
            } else {
              xs <- df[focal.genotype,]$x
              ys <- df[focal.genotype,]$y
              xe <- df[n,]$x
              ye <- df[n,]$y
              head.length <- 0
              arrow.type.num <- 3
            }
          }
          p<-arrow.prop.drawn
          delx <- (1-p)/2
          lenx <- abs(xs-xe)
          leny <- abs(ys-ye)
          slope <- leny / lenx
          len <- sqrt(1 + slope^2)
          if(slope==0) {
            pp<-p
          } else {
            pp <- 1 - (2*delx)/len
          }
          if(length(unique(merge.check)) != 1) {
            draw.arrow(x.start=xs,
                       y.start=ys,
                       x.end=xe,
                       y.end=ye,
                       prop.drawn=pp,
                       displacement=arrow.host.displacements[which(hosts==h)],
                       color=arrow.host.colors[which(hosts==h)],
                       length=head.length,
                       angle=30,
                       width = arrow.width,
                       type=arrow.type.num)
          }  
        }
      }
    }
  }
  new.df<-df[df$hosts==hosts[1],]
  for(g in 1:length(new.df$hosts)) {
    cols<-character(0)
    t.cols<-character(0)
    labs<-character(0)
    for(p in 1:N.loci) {
      if(new.df[g,1+p]==0) {
        cols<-c(cols,starting.color)
        t.cols<-c(t.cols,starting.text.color)
        labs<-c(labs,starting.position.labels[p])
      } else {
        cols<-c(cols,ending.color)
        t.cols<-c(t.cols,ending.text.color)
        labs<-c(labs,ending.position.labels[p])
      }
    }
    draw.pie(center.x = new.df[g,]$x,
             center.y = new.df[g,]$y,
             radius=node.radius,
             edge.color=node.edge.color,
             edge.width=node.edge.width,
             start.angle=pi/2+pi/N.loci,
             N.slices=N.loci,
             label.size=node.label.size,
             label.pos=label.pos,
             colors=cols,
             labels=labs,
             text.colors=t.cols,
             x.y.ratio=x.y.ratio)
  }
  ## draw blowups
  for(g in 1:length(new.df$hosts)) {
    num_of_mutations = sum(new.df[g,2:(1+N.loci)])
    if (num_of_mutations == 0 | num_of_mutations == N.loci){
      if (num_of_mutations == 0){
        xx <- c(new.df[g,]$x-node.radius, new.df[g,]$x+node.radius, new.df[g,]$x+blowup.node.radius, new.df[g,]$x-blowup.node.radius) #bottom left, bottom right, top right, top left
        yy <- c(new.df[g,]$y, new.df[g,]$y, new.df[g,]$y+blowup.distance, new.df[g,]$y+blowup.distance)
        y_store <- yy[3]
        polygon(xx, yy, col=adjustcolor(pie_color,alpha.f = 0.25) , border=NA)
      }
      if (num_of_mutations == N.loci){
        xx <- c(new.df[g,]$x-node.radius, new.df[g,]$x+node.radius, new.df[g,]$x+blowup.node.radius, new.df[g,]$x-blowup.node.radius)
        yy <- c(new.df[g,]$y, new.df[g,]$y, new.df[g,]$y-blowup.distance, new.df[g,]$y-blowup.distance)
        y_store <- yy[3]
        polygon(xx, yy, col=adjustcolor(pie_color,alpha.f = 0.25), border=NA)
      }
      cols<-character(0)
      t.cols<-character(0)
      labs<-character(0)
      blowuplabs<-character(0)
      for(p in 1:N.loci) {
        if(new.df[g,1+p]==0) {
          cols<-c(cols,starting.color)
          t.cols<-c(t.cols,starting.text.color)
          labs<-c(labs,starting.position.labels[p])
          blowuplabs<-c(blowuplabs,starting.blowup.labels[p])
        } else {
          cols<-c(cols,ending.color)
          t.cols<-c(t.cols,ending.text.color)
          labs<-c(labs,ending.position.labels[p])
          blowuplabs<-c(blowuplabs,ending.blowup.labels[p])
        }
      }
      draw.pie(center.x = new.df[g,]$x,
               center.y = new.df[g,]$y,
               radius=node.radius,
               edge.color=node.edge.color,
               edge.width=node.edge.width,
               start.angle=pi/2+pi/N.loci,
               N.slices=N.loci,
               label.size=node.label.size,
               label.pos=label.pos,
               colors=cols,
               labels=labs,
               text.colors=t.cols,
               x.y.ratio=x.y.ratio)
      draw.pie(center.x = new.df[g,]$x,
               center.y = y_store,
               radius=blowup.node.radius,
               edge.color=node.edge.color,
               edge.width=node.edge.width,
               start.angle=pi/2+pi/N.loci,
               N.slices=N.loci,
               label.size=node.label.size,
               label.pos=label.pos,
               colors=cols,
               labels=blowuplabs,
               text.colors=t.cols,
               x.y.ratio=x.y.ratio)
    }
  }
}

draw.figure.TEM.network.blowups <-function(network.x.y.ratio,
                                           buffer.x = 0.05,
                                           buffer.y = c(0.2,0,0,0.2)){
  df<-read.csv("inputs/resistance_levels.csv",header=TRUE)
  plot.new()
  draw.landscape.blowups(df,
                         N.loci=5,
                         starting.position.labels=c("","","","",""),
                         ending.position.labels=c("","","","",""),
                         starting.blowup.labels=c("g","A","E","M","G"),
                         ending.blowup.labels=c("a","G","K","T","S"),
                         N.fit.estimates=1,
                         buffer.x=buffer.x,
                         buffers.y=buffer.y,
                         node.radius=0.025,
                         blowup.node.radius=0.08,
                         blowup.distance=0.3,
                         node.edge.width=0.75,
                         node.label.size=1,
                         label.pos=.65,
                         node.edge.color='black',
                         starting.color="white",
                         ending.color=pie_color,
                         starting.text.color=pie_color,
                         ending.text.color="black",
                         merge_arrow_color=pie_arrow,
                         arrow.host.colors=c(red_Ec_color, blue_Kp_color, yellow_Se_color),
                         arrow.host.displacements=c(-0.012,0,0.012),
                         arrow.prop.drawn=0.6,
                         arrow.width=1.75,
                         arrow.head=.08,
                         gamma=1,
                         x.y.ratio=network.x.y.ratio)}

# Effect of mutations species correlations-----------
sp.slopes.df <- read.csv("inputs/mutational_step_slopes_only.csv")
sp.slopes <- sp.slopes.df %>%
  select(-X)

Ec.main.label = expression('mutation effect in ' * phantom('Ec'))
Ec.sp.label = expression(phantom('mutation effect in ') * 'Ec')
Ec.label = c(Ec.main.label, Ec.sp.label)

Kp.main.label = expression('mutation effect in ' * phantom('Kp'))
Kp.sp.label = expression(phantom('mutation effect in ') * 'Kp')
Kp.label = c(Kp.main.label,Kp.sp.label)

Se.main.label = expression('mutation effect in ' * phantom('Se'))
Se.sp.label = expression(phantom('mutation effect in ') * 'Se')
Se.label = c(Se.main.label,Se.sp.label)

Ec.Kp.no.epistasis <- sp.slopes %>%
  filter(Species %in% c('Ec','Kp')) %>%
  pivot_wider(names_from = Species, values_from = c(step_slope,step_effect)) %>%
  filter(step_effect_Ec == step_effect_Kp) %>%
  rename(Ec = step_slope_Ec, Kp = step_slope_Kp)

Ec.Kp.yes.epistasis <- sp.slopes %>%
  filter(Species %in% c('Ec','Kp')) %>%
  pivot_wider(names_from = Species, values_from = c(step_slope,step_effect)) %>%
  filter(step_effect_Ec != step_effect_Kp) %>%
  rename(Ec = step_slope_Ec, Kp = step_slope_Kp)

Ec.Se.no.epistasis <- sp.slopes %>%
  filter(Species %in% c('Ec','Se')) %>%
  pivot_wider(names_from = Species, values_from = c(step_slope,step_effect)) %>%
  filter(step_effect_Ec == step_effect_Se) %>%
  rename(Ec = step_slope_Ec, Se = step_slope_Se)

Ec.Se.yes.epistasis <- sp.slopes %>%
  filter(Species %in% c('Ec','Se')) %>%
  pivot_wider(names_from = Species, values_from = c(step_slope,step_effect)) %>%
  filter(step_effect_Ec != step_effect_Se) %>%
  rename(Ec = step_slope_Ec, Se = step_slope_Se)

Kp.Se.no.epistasis <- sp.slopes %>%
  filter(Species %in% c('Kp','Se')) %>%
  pivot_wider(names_from = Species, values_from = c(step_slope,step_effect)) %>%
  filter(step_effect_Kp == step_effect_Se) %>%
  rename(Kp = step_slope_Kp, Se = step_slope_Se)

Kp.Se.yes.epistasis <- sp.slopes %>%
  filter(Species %in% c('Kp','Se')) %>%
  pivot_wider(names_from = Species, values_from = c(step_slope,step_effect)) %>%
  filter(step_effect_Kp != step_effect_Se) %>%
  rename(Kp = step_slope_Kp, Se = step_slope_Se)

correlation_plot <- function(df.no, df.yes, x.host, y.host, x.label, x.color, y.label, y.color, no.color, yes.color, global.cex, tick.cex = 1) {
  plot((-2:22), (-2:22), type="n", xlab='', ylab='', frame.plot=F, axes = F)
  mtext(side=1, line=0.6, text = c(x.label[1], x.label[2]), col=c("black", x.color), cex=global.cex)
  mtext(side=2, line=0.6, text = c(y.label[1], y.label[2]), col=c("black", y.color), cex=global.cex)
  points(as.list(df.no[x.host])[[1]], as.list(df.no[y.host])[[1]], pch=16, col=no.color, cex = 0.75, lwd=0.5)
  points(as.list(df.yes[x.host])[[1]], as.list(df.yes[y.host])[[1]],pch=16, col=yes.color, cex = 0.75, lwd=0.5)
  axis(side = 1, pos = 0, at = c(5,10,15,20), col.ticks = 'black', col = 'white', cex.axis=tick.cex, padj = -1.5, tck=-(0.05/b_square_in))
  axis(side = 2, pos = 0, at = c(5,10,15,20), col.ticks = 'black', col = 'white', cex.axis=tick.cex, padj = 1.1, tck=-(0.05/b_square_in))
  arrows(x0 = 0, y0 = -2, x1 = 0, y1 = 22, lwd = 1, length = 0.05)
  arrows(x0 = -2, y0 = 0, x1 = 22, y1 = 0, lwd = 1, length = 0.05)
  lines(-2:22, -2:22, lwd = 1, col = "black")
}

# Fitness landscapes ----------
fitness.landscape.graph <- function(host, 
                                    host.color, 
                                    xyratio, 
                                    axis.label = F,
                                    radius.buffer,
                                    node.radius,
                                    global.cex) {
  df<-read.csv("inputs/resistance_levels.csv",header=TRUE)
  df.rel <- df %>%
    group_by(hosts) %>%
    mutate(MIC_WT = fit.1[Num_of_mutations == 0],
           fit.1 = fit.1 - MIC_WT)
  df.segments <- read.csv('inputs/mutational_step_slopes_rel.csv')
  df.host <- df.rel %>% filter(hosts == host)
  df.host.segments <- df.segments %>% filter(Species == host)
  n_mut <- max(df.rel$Num_of_mutations)
  plot(0:n_mut, seq((min(df.rel$fit.1)-0.5),(max(df.rel$fit.1)+0.5),length.out = 6), type="n", xaxt='n', xlab='', ylab='', axes=FALSE, frame.plot=TRUE)
  axis(1, at = 0:n_mut, cex.axis=tick.cex, padj = -1.5, tck=-(0.05/c_square_w_in), lwd = 0, lwd.ticks = 1)
  mtext(side=1, line=1.5, text = 'number of mutations', cex = global.cex)
  if (axis.label == T) {
    axis(2, at = seq(round(min(df.rel$fit.1)),round(max(df.rel$fit.1)),by = 5), cex = global.cex, cex.axis=tick.cex, padj = 1.1, tck=-(0.05/c_square_w_in), lwd = 0, lwd.ticks = 1)
    mtext(side=2, line=1.5, text = 'relative resistance level', cex = global.cex)
  } else {
    axis(2, at = seq(round(min(df.rel$fit.1)),round(max(df.rel$fit.1)),by = 5), labels = FALSE, cex = global.cex, cex.axis=tick.cex, padj = 1.1, tck=-(0.05/c_square_w_in), lwd = 0, lwd.ticks = 1)
  }
  starting.color="white"
  ending.color=pie_color
  N.loci = 5
  axis.multiplier = ((max(df.rel$fit.1)+0.5)-(min(df.rel$fit.1)-0.5))/N.loci
  node.edge.color = 'black'
  node.edge.width = 0.5
  label.pos = 0.5
  node.label.size = 0
  x.y.ratio = xyratio
  radius_filter = node.radius*radius.buffer*axis.multiplier*x.y.ratio
  for (s in 1:length(df.host.segments$step)){
    if (df.host.segments$step_effect[s] == 'neutral') {
      arrows(x0 = df.host.segments$focal_num[s], y0 = df.host.segments$focal_MIC[s], x1 = df.host.segments$mutant_num[s], y1 = df.host.segments$mutant_MIC[s],
             col = host.color, length=0, angle=30, code=1, lty = 3)
    }
    if (df.host.segments$step_effect[s] == 'deleterious') {
      arrows(x0 = df.host.segments$focal_num[s], y0 = df.host.segments$focal_MIC[s], x1 = df.host.segments$mutant_num[s], y1 = df.host.segments$mutant_MIC[s],
             col = host.color, length=0, angle=30, code=1, lty = 5)
    }
    if (df.host.segments$step_effect[s] == 'beneficial') {
      arrows(x0 = df.host.segments$focal_num[s], y0 = df.host.segments$focal_MIC[s], x1 = df.host.segments$mutant_num[s], y1 = df.host.segments$mutant_MIC[s],
             col = host.color, length=0, angle=30, code=1, lty = 1)
    }
  } 
  df.host <- df.host %>%
    group_by(Num_of_mutations) %>%
    arrange(fit.1)%>%
    mutate(N = n(),
           closest_below = fit.1 - lag(fit.1), 
           closest_above = lead(fit.1) - fit.1,
           filter_below = ifelse(is.na(closest_below) | closest_below > radius_filter, T, F),
           filter_above = ifelse(is.na(closest_above) | closest_above > radius_filter, T, F),
           pie_filter = ifelse(filter_below == F | filter_above == F, F, T))
  df.host.filter <- df.host %>%
    filter(pie_filter == T)
  offsets <- read.csv(file = "inputs/Offsets.csv")
  offsets <- offsets %>%
    select(hosts, g, A, E, M, G, xC, yC)
  df.host$xC <- NULL
  df.host$yC <- NULL
  df.host <- left_join(df.host, offsets)
  for(g in 1:length(df.host$hosts)) {
    if (df.host[g,]$pie_filter == F) {
      mA = (df.host[g,]$yC-df.host[g,]$fit.1)/(df.host[g,]$xC-df.host[g,]$Num_of_mutations)
      mB = -1/mA
      deltaX = node.radius*cos(atan(mB))
      deltaY = node.radius*sin(atan(mB))*axis.multiplier*xyratio
      x.values <- c(df.host[g,]$Num_of_mutations, df.host[g,]$xC+deltaX, df.host[g,]$xC-deltaX)
      y.values <- c (y0 = df.host[g,]$fit.1, df.host[g,]$yC+deltaY, df.host[g,]$yC-deltaY)
      polygon(x.values,y.values,border=NA,col=adjustcolor(pie_color,alpha.f = 0.45),lwd=0.1)
    }
  }
  for(g in 1:length(df.host$hosts)) {
    if (df.host[g,]$pie_filter == F) {
      cols<-character(0)
      for(p in 1:N.loci) {
        if(df.host[g,1+p]==0) {
          cols<-c(cols,starting.color)
        } else {
          cols<-c(cols,ending.color)
        }
      }
      draw.pie.yaxis(center.x = df.host[g,]$xC,
                     center.y = df.host[g,]$yC,
                     radius=node.radius,
                     edge.color=node.edge.color,
                     edge.width=node.edge.width,
                     start.angle=pi/2+pi/N.loci,
                     N.slices=N.loci,
                     label.size=node.label.size,
                     label.pos=label.pos,
                     colors=cols,
                     labels='',
                     text.colors='black',
                     axis.multiplier=((max(df.rel$fit.1)+0.5)-(min(df.rel$fit.1)-0.5))/N.loci,
                     xy.ratio = xyratio)
    }
  }
  points(df.host$Num_of_mutations, df.host$fit.1, pch = 21, col="black", bg=pie_color,cex = 0.5,lwd=0.4)
  for(g in 1:length(df.host$hosts)) {
    if (df.host[g,]$pie_filter == T) {
      cols<-character(0)
      for(p in 1:N.loci) {
        if(df.host[g,1+p]==0) {
          cols<-c(cols,starting.color)
        } else {
          cols<-c(cols,ending.color)
        }
      }
      draw.pie.yaxis(center.x = df.host[g,]$Num_of_mutations,
                     center.y = df.host[g,]$fit.1,
                     radius=node.radius,
                     edge.color=node.edge.color,
                     edge.width=node.edge.width,
                     start.angle=pi/2+pi/N.loci,
                     N.slices=N.loci,
                     label.size=node.label.size,
                     label.pos=label.pos,
                     colors=cols,
                     labels='',
                     text.colors='black',
                     axis.multiplier=((max(df.rel$fit.1)+0.5)-(min(df.rel$fit.1)-0.5))/N.loci,
                     xy.ratio = xyratio)
    }
  }
}

# Evolutionary simulations -----------
simulation_data_load <- function(filename, reference_filename){
  df_all_only <- read.csv(filename)
  df_ref <- read.csv(reference_filename)
  df_ref <- df_ref %>%
    select(-Shorthand, - fit.1, -fit.2, -outlier, -fit.log) %>%
    rename(Species_order = Species)
  df_only <- left_join(df_all_only, df_ref)
  df_only <- df_only %>%
    select(-X, -Resistance) %>%
    rename(Resistance = fit_average) %>%
    group_by(Treatment_ID, Cumulative_time, Species_time, Species_order) %>%
    summarise(N = n(),
              Resistance.mean = mean(Resistance), 
              Resistance.sd = sd(Resistance),
              Resistance.se = Resistance.sd/sqrt(N))
  return(df_only)
}


simulation_plot <- function(df,
                            x.main.label,
                            global.cex,
                            y.main.label,
                            focal_color,
                            x_max,
                            y_max) {
  
  plot((0:x_max), seq(0, y_max, y_max/x_max), type="n", xlab='', ylab='', frame.plot=T, axes = F)  
  mtext(side=1, line=1.5, text =x.main.label, col='black', cex=global.cex)
  mtext(side=2, line=1.5, text =y.main.label, col='black',  cex=global.cex)
  axis(2, cex.axis=tick.cex, padj = 1.1, tck=-(0.05/d_square_h_in), lwd = 0, lwd.ticks = 1)
  axis(1, cex.axis=tick.cex, padj = -1.5, tck=-(0.05/d_square_h_in), lwd = 0, lwd.ticks = 1)
  
  polygon(x = c(df$Cumulative_time,rev(df$Cumulative_time)), y = c(df$Resistance.mean-df$Resistance.se, rev(df$Resistance.mean+df$Resistance.se)), col=adjustcolor(focal_color,alpha.f = 0.25), border=NA)
  lines(df$Cumulative_time, df$Resistance.mean, col = focal_color, lwd = 1.5)
}

simulation_plot_split <- function(df,
                                  df.mean,
                                  x.main.label,
                                  global.cex,
                                  y.main.label,
                                  combine_color,
                                  focal_color,
                                  transient_color,
                                  x_max,
                                  y_max,
                                  arrow_proportion, 
                                  line_proportion_multiplier,
                                  buffer_proportion) {
  
  df_early <- df %>% filter(Cumulative_time <= 20)
  df_late <- df %>% filter(Cumulative_time > 40)
  df_middle <- df %>% filter(Cumulative_time <= 40 & Cumulative_time > 20)
  plot((0:x_max), seq(0, y_max, y_max/x_max), type="n", xlab='', ylab='', frame.plot=T, axes = F)
  axis(2, labels = FALSE, cex.axis=tick.cex, padj = -1.5, tck=-(0.05/d_square_h_in), lwd = 0, lwd.ticks = 1)
  axis(1, cex.axis=tick.cex, padj = -1.5, tck=-(0.05/d_square_h_in), lwd = 0, lwd.ticks = 1)
  mtext(side=1, line=1.5, text =x.main.label, col='black', cex=global.cex)
  mtext(side=2, line=1.5, text =y.main.label, col='black',  cex=global.cex)
  y_arrow_length = y_max/arrow_proportion
  x_arrow_length = x_max/arrow_proportion
  line_proportion = arrow_proportion*line_proportion_multiplier
  y_distance_from_arrowpoint = y_max/line_proportion
  x_distance_from_arrowpoint = x_max/line_proportion
  hgt_times = c(20, 40)
  polygon(x = c(df_early$Cumulative_time,rev(df_early$Cumulative_time)), y = c(df_early$Resistance.mean-df_early$Resistance.se, rev(df_early$Resistance.mean+df_early$Resistance.se)), col=adjustcolor(focal_color,alpha.f = 0.25), border=NA)
  lines(df_early$Cumulative_time, df_early$Resistance.mean, col = focal_color, lwd = 1.5)
  polygon(x = c(df_middle$Cumulative_time,rev(df_middle$Cumulative_time)), y = c(df_middle$Resistance.mean-df_middle$Resistance.se, rev(df_middle$Resistance.mean+df_middle$Resistance.se)), col=adjustcolor(transient_color,alpha.f = 0.25), border=NA)
  lines(df_middle$Cumulative_time, df_middle$Resistance.mean, col = transient_color, lwd = 1.5)
  polygon(x = c(df_late$Cumulative_time,rev(df_late$Cumulative_time)), y = c(df_late$Resistance.mean-df_late$Resistance.se, rev(df_late$Resistance.mean+df_late$Resistance.se)), col=adjustcolor(focal_color,alpha.f = 0.25), border=NA)
  lines(df_late$Cumulative_time, df_late$Resistance.mean, col = focal_color, lwd = 1.5)
  df_mean <- df.mean %>% filter(Cumulative_time == x_max) %>% select(Resistance.mean)
  arrows(x0 = 0, y0 = df_mean$Resistance.mean, x1 = x_max, y1 = df_mean$Resistance.mean, length = 0, col = focal_color, lty = 'dashed')
  
  for (t in c(1:length(hgt_times))) { 
    distance_from_endpoint = y_max/buffer_proportion[t]
    if (df$Resistance.mean[df$Cumulative_time==hgt_times[t]]>df$Resistance.mean[df$Cumulative_time==hgt_times[t]+1])  {
      arrow_delta_x = x_arrow_length*cos(45)
      arrow_delta_y = y_arrow_length*sin(45)
      y0 = df$Resistance.mean[df$Cumulative_time==hgt_times[t]+1]+distance_from_endpoint
      arrows(x0 = hgt_times[t]+0.5, y0 = y0, x1 = hgt_times[t]+0.5+arrow_delta_x, y1 = y0+arrow_delta_y, length = 0, col = combine_color)
      arrows(x0 = hgt_times[t]+0.5, y0 = y0, x1 = hgt_times[t]+0.5-arrow_delta_x, y1 = y0+arrow_delta_y, length = 0, col = combine_color)
      line_delta_x = x_distance_from_arrowpoint*cos(45)
      line_delta_y = y_distance_from_arrowpoint*sin(45)
      y1 = df$Resistance.mean[df$Cumulative_time==hgt_times[t]]-distance_from_endpoint
      arrows(x0 = hgt_times[t]+0.5+line_delta_x, y0 = y0+line_delta_y, x1 = hgt_times[t]+0.5+line_delta_x, y1 = y1, length = 0, col = combine_color)
      arrows(x0 = hgt_times[t]+0.5-line_delta_x, y0 = y0+line_delta_y, x1 = hgt_times[t]+0.5-line_delta_x, y1 = y1, length = 0, col = combine_color)
    } else {
      arrow_delta_x = x_arrow_length*cos(45)
      arrow_delta_y = y_arrow_length*sin(45)
      y0 = df$Resistance.mean[df$Cumulative_time==hgt_times[t]+1]-distance_from_endpoint
      arrows(x0 = hgt_times[t]+0.5, y0 = y0, x1 = hgt_times[t]+0.5+arrow_delta_x, y1 = y0-arrow_delta_y, length = 0, col = combine_color)
      arrows(x0 = hgt_times[t]+0.5, y0 = y0, x1 = hgt_times[t]+0.5-arrow_delta_x, y1 = y0-arrow_delta_y, length = 0, col = combine_color)
      line_delta_x = x_distance_from_arrowpoint*cos(45)
      line_delta_y = y_distance_from_arrowpoint*sin(45)
      y1 = df$Resistance.mean[df$Cumulative_time==hgt_times[t]]+distance_from_endpoint
      arrows(x0 = hgt_times[t]+0.5+line_delta_x, y0 = y0-line_delta_y, x1 = hgt_times[t]+0.5+line_delta_x, y1 = y1, length = 0, col = combine_color)
      arrows(x0 = hgt_times[t]+0.5-line_delta_x, y0 = y0-line_delta_y, x1 = hgt_times[t]+0.5-line_delta_x, y1 = y1, length = 0, col = combine_color)
    }
  }
}


# Arrange plot --------------
### Figure 3  --------------
matrix_layout <- as.matrix(read.csv('inputs/matrix.csv', header = F))
n.row = nrow(matrix_layout)
n.col = ncol(matrix_layout)
outer_margin = 0.01
tick_marks = 0.05
# columns
c1_left_margin = 0.45 + tick_marks
c1_right_margin = 0.25 
c2_left_margin = 0 + tick_marks
c2_right_margin = 0.25
c3_left_margin = 0 + tick_marks
c3_right_margin = 0.13
# rows
r1_top_margin = 0
r1_bottom_margin = 0.35
r2_top_margin = 0
r2_bottom_margin = 0.35
r3_top_margin = 0
r3_bottom_margin = 0.35
r4_top_margin = 0
r4_bottom_margin = 0.4 + tick_marks
r5_top_margin = 0
r5_bottom_margin = 0.4 + tick_marks
# b part dimensions
b_square_in = 1.25
b_left_margin = 0.6
b_right_margin = 1.93-0.6-1.25
# c part dimensions
c_square_h_in = 2
c_square_w_in = 1.75
# d part dimensions
d_square_h_in = 1.5
d_square_w_in = 1.75
c1 = c1_left_margin + c1_right_margin + c_square_w_in
c2 = c2_left_margin + c2_right_margin + c_square_w_in  # 2.05
c3 = c3_left_margin + c3_right_margin + c_square_w_in # 1.94
r1 = r1_bottom_margin + r1_top_margin + b_square_in  # 1.6
r2 = r2_bottom_margin + r2_top_margin + b_square_in  # 1.6
r3 = r3_bottom_margin + r3_top_margin + b_square_in # 1.6
r4 = r4_bottom_margin + r4_top_margin + c_square_h_in  # 2.45
r5 = r5_bottom_margin + r5_top_margin + d_square_h_in  # 1.95
pdf.width = c1 + c2 + c3 + (outer_margin*2) # 6.5 max
pdf.height = r1 + r2 + r3 + r4 + r5 + (outer_margin*2) # 9.25 max
### Initiate  --------
pointsize = 12
globalpointsize = 9
tickpointsize = 7
pdf(file = 'Fig3.pdf', width = pdf.width, height = pdf.height, pointsize = pointsize)
layout(matrix(matrix_layout, nrow = n.row, ncol = n.col),
       widths = c(lcm(c1*2.54),lcm(c2*2.54),lcm(c3*2.54)), 
       heights = c(lcm(r1*2.54),lcm(r2*2.54),lcm(r3*2.54),lcm(r4*2.54),lcm(r5*2.54)))
par(omi=c(outer_margin,outer_margin,outer_margin,outer_margin)) # outer margins (b, l, t, r)
global.cex = globalpointsize / pointsize
tick.cex = tickpointsize / globalpointsize
### Panel A----------
par(mai=c(0,0,0,0), cex = global.cex) # margins
panel.A.x.y.ratio = (c1 + c2) / (r1 + r2 + r3)
draw.figure.TEM.network.blowups(panel.A.x.y.ratio, buffer.y = c(0.21,0.01,0.01,0.21))
mtext("a", line = -1.1, side=3, at = -0.02) 

### Panel B---------
par(mai=c(r1_bottom_margin,b_left_margin,r1_top_margin,b_right_margin), cex = global.cex) # margins
correlation_plot(df.no = Ec.Kp.no.epistasis, df.yes = Ec.Kp.yes.epistasis, x.host = 'Ec', y.host = 'Kp', x.label = Ec.label, y.label = Kp.label, x.color = red_Ec_color, y.color = blue_Kp_color, no.color = pie_arrow, yes.color = purple_Ec.Kp_color, global.cex = global.cex, tick.cex = tick.cex)
mtext("b", line = -1.1, side=3, at = -13)  
correlation_plot(df.no = Kp.Se.no.epistasis, df.yes = Kp.Se.yes.epistasis, x.host = 'Kp', y.host = 'Se', x.label = Kp.label, y.label = Se.label, x.color = blue_Kp_color, y.color = yellow_Se_color, no.color = pie_arrow, yes.color = green_Se.Kp_color, global.cex = global.cex, tick.cex = tick.cex)
correlation_plot(df.no = Ec.Se.no.epistasis, df.yes = Ec.Se.yes.epistasis, x.host = 'Se', y.host = 'Ec', x.label =  Se.label, y.label = Ec.label, x.color = yellow_Se_color, y.color = red_Ec_color, no.color = pie_arrow, yes.color = orange_Ec.Se_color, global.cex = global.cex, tick.cex = tick.cex)

### Panel C------------
par(mai=c(r4_bottom_margin,c1_left_margin,r4_top_margin,c1_right_margin), cex = global.cex)
xyratio.pdf <- c_square_w_in/c_square_h_in
radiusbuffer = 1.2 
noderadius = 0.125 #0.125
fitness.landscape.graph(host = 'Ec',
                        host.color = red_Ec_color,
                        xyratio = xyratio.pdf, 
                        axis.label = T,
                        radius.buffer = radiusbuffer,
                        node.radius = noderadius,
                        global.cex = global.cex)
mtext("c", line = -0.5, side=3, at = -1.5, cex = 1)  
par(mai=c(r4_bottom_margin,c2_left_margin,r4_top_margin,c2_right_margin), cex = global.cex)
fitness.landscape.graph(host = 'Kp',
                        host.color = blue_Kp_color,
                        xyratio = xyratio.pdf,
                        axis.label = F,
                        radius.buffer = radiusbuffer,
                        node.radius = noderadius,
                        global.cex = global.cex)
par(mai=c(r4_bottom_margin,c3_left_margin,r4_top_margin,c3_right_margin), cex = global.cex)
fitness.landscape.graph(host = 'Se',
                        host.color = yellow_Se_color,
                        xyratio = xyratio.pdf,
                        axis.label = F,
                        radius.buffer = radiusbuffer,
                        node.radius = noderadius,
                        global.cex = global.cex)


### Panel D, E, G---------
sim_ymax = 875
par(mai=c(r5_bottom_margin,c1_left_margin,r5_top_margin,c1_right_margin), cex = global.cex)

df_Ec_only <- simulation_data_load(filename = "inputs/76.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Ec_Kp <- simulation_data_load(filename = "inputs/77.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Ec_Se <- simulation_data_load(filename = "inputs/78.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
simulation_plot(df = df_Ec_only, 
                x.main.label = 'time', 
                global.cex = global.cex, 
                y.main.label = "resistance", 
                focal_color = red_Ec_color, 
                x_max = 60, 
                y_max = sim_ymax)
mtext("d", line = -0.5, side=3, at = -18.5)  
par(mai=c(r5_bottom_margin,c2_left_margin,r5_top_margin,c2_right_margin), cex = global.cex)
arrow_prop = 17.5
line_prop_mult = 4
buffer_prop = c(1000, 40) # was 80
simulation_plot_split(df = df_Ec_Kp,
                      df.mean = df_Ec_only, 
                      x.main.label = 'time', 
                      global.cex = global.cex, 
                      y.main.label = "", 
                      combine_color = purple_Ec.Kp_color, 
                      focal_color = red_Ec_color, 
                      transient_color = blue_Kp_color, 
                      x_max = 60, 
                      y_max = sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
mtext("e", line = -0.5, side=3, at = -7.5)  
par(mai=c(r5_bottom_margin,c3_left_margin,r5_top_margin,c3_right_margin), cex = global.cex)
buffer_prop = c(80, 20)
simulation_plot_split(df = df_Ec_Se, 
                      df.mean = df_Ec_only,
                      x.main.label = 'time', 
                      global.cex = global.cex, 
                      y.main.label = "", 
                      combine_color = orange_Ec.Se_color, 
                      focal_color = red_Ec_color, 
                      transient_color = yellow_Se_color, 
                      x_max = 60, 
                      y_max = sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
mtext("f", line = -0.5, side=3, at = -7.5)  
dev.off()
