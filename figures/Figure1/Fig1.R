# Packages-----
library(tidyverse)
library(cowplot)
library(broom)
library(drc)
library(modelr)
library(pBrackets)
theme_set(theme_cowplot())
# Colors ---------
red_Ec_color <- rgb(255, 48, 48, maxColorValue = 255)
blue_Kp_color <- rgb(0, 0, 238, maxColorValue = 255)
yellow_Se_color <- rgb(242, 184, 0, maxColorValue = 255)
purple_Ec.Kp_color <- rgb(178, 58, 238, maxColorValue = 255)
orange_Ec.Se_color <- rgb(238, 118, 0, maxColorValue = 255)
green_Se.Kp_color <- rgb(70, 178, 0, maxColorValue = 255)
pie_color <- rgb(169, 169, 169, maxColorValue = 255)
pie_arrow <- rgb(226, 210, 195, maxColorValue = 255)
gene_color <- rgb(245, 245, 245, maxColorValue = 255)
# Functions -----
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

# Landscape plot ------
fitness.landscape.graph <- function(df_filename,
                                    offsets_filename,
                                    host, 
                                    host.color, 
                                    xyratio, 
                                    axis.label = F,
                                    radius.buffer,
                                    node.edge.width,
                                    global.cex,
                                    node.radius) {
  df_format<-read.csv(df_filename,header=TRUE)
  data <- read.csv(df_filename, header = T)
  steps_file <- read.csv('inputs/mutational_steps.csv', header = T)
  steps <- steps_file %>%
    mutate(step = paste(focal, mutant, sep = "->"))
  focal_node <- steps %>%
    dplyr::select(-mutant) %>%
    rename(Genotype = focal)
  mutant_step <- steps  %>%
    dplyr::select(-focal) %>%
    rename(Genotype = mutant)
  focal <- left_join(data, focal_node)
  focal$identity <- 'focal'
  mutant <- left_join(data, mutant_step)
  mutant$identity <- 'mutant'
  full_data <- full_join(focal, mutant) 
  full <- full_data %>%
    separate(col = step, sep = '->', into = c('focal', 'mutant'), remove = F) %>%
    filter(!(is.na(focal))) %>%
    group_by(Species, step) %>%
    mutate(focal_MIC = MIC_mean[Genotype == focal], 
           mutant_MIC = MIC_mean[Genotype == mutant],
           focal_se = MIC_se[Genotype == focal],
           mutant_se = MIC_se[Genotype == mutant],
           step_slope = (mutant_MIC) - (focal_MIC), 
           step_effect = case_when((mutant_MIC-mutant_se) > (focal_MIC+focal_se) ~ 'beneficial',
                                   (mutant_MIC+mutant_se) < (focal_MIC-focal_se) ~ 'deleterious',
                                   TRUE ~ 'neutral')) %>%
    group_by(step) %>%
    mutate(host_epistasis = ifelse(length(unique(step_effect)) == 1, F, T))
  step_data <- full %>%
    group_by(Species)%>%
    filter(!(is.na(step)),
           identity == 'focal')
  df.segments <- step_data %>%
    mutate(focal_num = nchar(focal),
           mutant_num = nchar(mutant),
           focal_num = ifelse(focal == 'WT', 0, focal_num)) %>%
    dplyr::select(Species, step, step_slope, step_effect, focal_MIC, mutant_MIC, focal_num, mutant_num)
  df<-df_format %>%
    rename(hosts = Species,
           fit.1 = MIC_mean) %>%
    dplyr::select(-Genotype)
  df.rel <- df %>%
    group_by(hosts) %>%
    mutate(MIC_WT = fit.1[Num_of_mutations == 0],
           fit.1 = fit.1)
  df.host <- df.rel %>% filter(hosts == host)
  df.host.segments <- df.segments %>% filter(Species == host)
  n_mut <- max(df.rel$Num_of_mutations)
  plot(seq(-0.25,n_mut+0.25,length.out = 2), seq((min(df.rel$fit.1)-1),(max(df.rel$fit.1)+1),length.out = 2), type="n", xaxt='n', xlab='', ylab='', axes=FALSE, frame.plot=TRUE)
  axis(1, at = 0:n_mut, cex.axis = 0.9, tck=-0.04, padj = -1)
  axis(2, labels = FALSE, tick = F, cex = global.cex)
  if (axis.label == T) {
    mtext(side=2, line=0.5, text = 'level of resistance', cex = global.cex)
  }
  mtext(side=1, line=1.55, text = 'number of mutations', cex = global.cex)
  
  starting.color="white"
  ending.color=pie_color
  N.loci = 3
  axis.multiplier = ((max(df.rel$fit.1)+0.1)-(min(df.rel$fit.1)-0.1))/N.loci
  node.edge.color = 'black'
  node.edge.width = node.edge.width
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
  offsets <- read.csv(file = offsets_filename)
  offsets <- offsets %>%
    dplyr::select(hosts, A, B, C, xC, yC)
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
                     axis.multiplier=axis.multiplier,
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
                     axis.multiplier=axis.multiplier,
                     xy.ratio = xyratio)
    }
  }
}

# Evolution plot -----
mutation_plot <- function(df_filename,
                          df_trajectories_filename,
                          treatment_filter,
                          x.main.label,
                          focal_color,
                          line_thickness,
                          node.edge.width,
                          xyratio,
                          radius.buffer,
                          global.cex,
                          node.radius,
                          arrow_length) {
  node.label.size = 1
  label.pos = 0.5
  node.edge.color = 'black'
  N.loci = 3
  df_MIC = read.csv(df_filename)
  df_trajectories = read.csv(df_trajectories_filename)
  df_trajectories_filter = df_trajectories %>%
    filter(treatment == treatment_filter)
  
  df_traj_resistance <- left_join(df_trajectories_filter, df_MIC)
  x_max = max(df_trajectories$time)+1
  x_min = -0.25
  y_max = (max(df_MIC$MIC_mean))+1
  y_min = (min(df_MIC$MIC_mean))-1
  
  plot(seq(x_min, x_max, length=2), seq(y_min, y_max, length=2), type="n", xlab='', ylab='', frame.plot=T, axes = F)
  arrows(x0 = 0, x1 = x_max, y0 = y_min - 1.5, y1 = y_min - 1.5, lwd = line_thickness, length = arrow_length, xpd = TRUE) # x-axis arrow
  mtext(side=1, line=1.55, text =x.main.label, col='black', cex=global.cex)

  starting.color="white"
  ending.color=pie_color
  axis.multiplier = (y_max-y_min)/(x_max - x_min)
  
  for(g in 2:length(df_traj_resistance$Species)) {
    if ((df_traj_resistance[g,]$MIC_mean - df_traj_resistance[g-1,]$MIC_mean) < 0) {
      current_time = df_traj_resistance$time[g]
      lines(x = c(current_time-1+node.radius,max(df_traj_resistance$time)+1), y = rep(df_traj_resistance$MIC_mean[g-1], 2), col = focal_color, lwd = line_thickness)
      cols<-character(0)
      for(p in 5:(5+N.loci-1)) {
        if(df_traj_resistance[g,1+p]==0) {
          cols<-c(cols,starting.color)
        } else {
          cols<-c(cols,ending.color)
        }
      }
      draw.pie.yaxis(center.x = df_traj_resistance[g-1,]$time,
                     center.y = df_traj_resistance[g-1,]$MIC_mean,
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
                     axis.multiplier=axis.multiplier,
                     xy.ratio = xyratio)
      break
    } else {
      current_time = df_traj_resistance$time[g]
      lines(x = c(current_time-1+node.radius,current_time), y = rep(df_traj_resistance$MIC_mean[g-1], 2), col = focal_color, lwd = line_thickness)
      if (current_time == max(df_traj_resistance$time)){
        lines(x = c(current_time+node.radius,current_time+1), y = rep(df_traj_resistance$MIC_mean[g], 2), col = focal_color, lwd = line_thickness)
      }
      if (df_traj_resistance[g,]$event == 'mutation') {
        arrows(x0 = df_traj_resistance[g,]$time, x1 = df_traj_resistance[g,]$time,
               y0 = df_traj_resistance[g-1,]$MIC_mean,
               y1 = df_traj_resistance[g,]$MIC_mean-(node.radius*axis.multiplier),
               length = arrow_length, lwd = line_thickness, col = focal_color)
      }
      cols<-character(0)
      for(p in 5:(5+N.loci-1)) {
        if(df_traj_resistance[g,1+p]==0) {
          cols<-c(cols,starting.color)
        } else {
          cols<-c(cols,ending.color)
        }
      }
      draw.pie.yaxis(center.x = df_traj_resistance[g,]$time,
                     center.y = df_traj_resistance[g,]$MIC_mean,
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
                     axis.multiplier=axis.multiplier,
                     xy.ratio = xyratio)
      if (current_time == max(df_traj_resistance$time)) {
        cols<-character(0)
        for(p in 5:(5+N.loci-1)) {
          if(df_traj_resistance[g,1+p]==0) {
            cols<-c(cols,starting.color)
          } else {
            cols<-c(cols,ending.color)
          }
        }
        draw.pie.yaxis(center.x = df_traj_resistance[g,]$time,
                       center.y = df_traj_resistance[g,]$MIC_mean,
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
                       axis.multiplier=axis.multiplier,
                       xy.ratio = xyratio)
      }
    }
  }
  cols<-character(0)
  for(p in 5:(5+N.loci-1)) {
    if(df_traj_resistance[1,1+p]==0) {
      cols<-c(cols,starting.color)
    } else {
      cols<-c(cols,ending.color)
    }
  }
  draw.pie.yaxis(center.x = df_traj_resistance[1,]$time,
                 center.y = df_traj_resistance[1,]$MIC_mean,
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
                 axis.multiplier=axis.multiplier,
                 xy.ratio = xyratio)
}

# df_filename = 'inputs/aligned.csv'
# df_trajectories_filename = 'inputs/Trajectory.csv'
# treatment_filter = 'HGT'
# x.main.label = 'time'
# focal_color = red_Ec_color
# transient_color = blue_Kp_color
# combine_color = purple_Ec.Kp_color
# line_thickness = 1
# node.edge.width = 0.5
# xyratio = 1
# radius.buffer = 1
# global.cex = global.cex
# node.radius = 0.15
# arrow_length = 0.075
# arrow_proportion = 35
# line_proportion_multiplier = 4
# buffer_proportion = 80

HGT_plot <- function(df_filename,
                          df_trajectories_filename,
                          treatment_filter,
                          x.main.label,
                          focal_color,
                          transient_color,
                          combine_color,
                          line_thickness,
                          node.edge.width,
                          xyratio,
                          radius.buffer,
                          global.cex,
                          node.radius,
                          arrow_proportion,
                          line_proportion_multiplier,
                          buffer_proportion,
                          arrow_length) {
  node.label.size = 1
  label.pos = 0.5
  node.edge.color = 'black'
  N.loci = 3
  df_MIC = read.csv(df_filename)
  df_trajectories = read.csv(df_trajectories_filename)
  df_trajectories_filter = df_trajectories %>%
    filter(treatment == treatment_filter)
  df_traj_resistance <- left_join(df_trajectories_filter, df_MIC)
  x_max = max(df_trajectories$time)+1
  x_min = -0.25
  y_max = (max(df_MIC$MIC_mean))+1
  y_min = (min(df_MIC$MIC_mean))-1
  plot(seq(x_min, x_max, length=2), seq(y_min, y_max, length=2), type="n", xlab='', ylab='', frame.plot=T, axes = F)
  arrows(x0 = 0, x1 = x_max, y0 = y_min - 1.5, y1 = y_min - 1.5, lwd = line_thickness, length = arrow_length, xpd = TRUE) # x-axis arrow
  mtext(side=1, line=1.55, text =x.main.label, col='black', cex=global.cex)
  starting.color="white"
  ending.color=pie_color
  axis.multiplier = (y_max-y_min)/(x_max - x_min)
  
  for(g in 2:length(df_traj_resistance$Species)) {
      if (df_traj_resistance$Species[g-1] == 'Ec'){
        lines(x = c(df_traj_resistance$time[g-1],df_traj_resistance$time[g]), y = rep(df_traj_resistance$MIC_mean[g-1], 2), col = focal_color, lwd = line_thickness)
      } else {
        lines(x = c(df_traj_resistance$time[g-1],df_traj_resistance$time[g]), y = rep(df_traj_resistance$MIC_mean[g-1], 2), col = transient_color, lwd = line_thickness)
      }
      if (df_traj_resistance$time[g] == max(df_traj_resistance$time)){
        lines(x = c(df_traj_resistance$time[g],df_traj_resistance$time[g]+1), y = rep(df_traj_resistance$MIC_mean[g], 2), col = focal_color, lwd = line_thickness)
      }
      if (df_traj_resistance[g,]$event == 'mutation') {
        if (df_traj_resistance$Species[g-1] == 'Ec') {
          arrows(x0 = df_traj_resistance[g,]$time, x1 = df_traj_resistance[g,]$time,
                 y0 = df_traj_resistance[g-1,]$MIC_mean,
                 y1 = df_traj_resistance[g,]$MIC_mean-(node.radius*axis.multiplier),
                 length = arrow_length, lwd = line_thickness, col = focal_color)
        } else {
          arrows(x0 = df_traj_resistance[g,]$time, x1 = df_traj_resistance[g,]$time,
                 y0 = df_traj_resistance[g-1,]$MIC_mean,
                 y1 = df_traj_resistance[g,]$MIC_mean-(node.radius*axis.multiplier),
                 length = arrow_length, lwd = line_thickness, col = transient_color)
        }
      }
      if (df_traj_resistance[g-1,]$event == 'mutation') {    
        cols<-character(0)
        for(p in 5:(5+N.loci-1)) {
          if(df_traj_resistance[g,1+p]==0) {
            cols<-c(cols,starting.color)
          } else {
            cols<-c(cols,ending.color)
          }
        }
        draw.pie.yaxis(center.x = df_traj_resistance[g-1,]$time,
                       center.y = df_traj_resistance[g-1,]$MIC_mean,
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
                       axis.multiplier=axis.multiplier,
                       xy.ratio = xyratio)
      }
    if (df_traj_resistance[g,]$event == 'HGT') {
      y_arrow_length = y_max/arrow_proportion
      x_arrow_length = x_max/arrow_proportion
      line_proportion = arrow_proportion*line_proportion_multiplier
      y_distance_from_arrowpoint = y_max/line_proportion
      x_distance_from_arrowpoint = x_max/line_proportion
      distance_from_endpoint = y_max/buffer_proportion
      if (df_traj_resistance$MIC_mean[g-1] > df_traj_resistance$MIC_mean[g]) {
        arrow_delta_x = x_arrow_length*cos(45)
        arrow_delta_y = y_arrow_length*sin(45)
        y0 = df_traj_resistance[g,]$MIC_mean+distance_from_endpoint
        arrows(x0 = df_traj_resistance[g,]$time, y0 = y0, x1 = df_traj_resistance[g,]$time+arrow_delta_x, y1 = y0+arrow_delta_y, length = 0, col = combine_color)
        arrows(x0 = df_traj_resistance[g,]$time, y0 = y0, x1 = df_traj_resistance[g,]$time-arrow_delta_x, y1 = y0+arrow_delta_y, length = 0, col = combine_color)
        line_delta_x = x_distance_from_arrowpoint*cos(45)
        line_delta_y = y_distance_from_arrowpoint*sin(45)
        y1 = df_traj_resistance[g-1,]$MIC_mean-distance_from_endpoint
        arrows(x0 = df_traj_resistance[g,]$time+line_delta_x, y0 = y0+line_delta_y, x1 = df_traj_resistance[g,]$time+line_delta_x, y1 = y1, length = 0, col = combine_color)
        arrows(x0 = df_traj_resistance[g,]$time-line_delta_x, y0 = y0+line_delta_y, x1 = df_traj_resistance[g,]$time-line_delta_x, y1 = y1, length = 0, col = combine_color)
      } else {
        arrow_delta_x = x_arrow_length*cos(45)
        arrow_delta_y = y_arrow_length*sin(45)
        y0 = df_traj_resistance[g,]$MIC_mean-distance_from_endpoint
        arrows(x0 = df_traj_resistance[g,]$time, y0 = y0, x1 = df_traj_resistance[g,]$time+arrow_delta_x, y1 = y0-arrow_delta_y, length = 0, col = combine_color)
        arrows(x0 = df_traj_resistance[g,]$time, y0 = y0, x1 = df_traj_resistance[g,]$time-arrow_delta_x, y1 = y0-arrow_delta_y, length = 0, col = combine_color)
        line_delta_x = x_distance_from_arrowpoint*cos(45)
        line_delta_y = y_distance_from_arrowpoint*sin(45)
        y1 = df_traj_resistance[g-1,]$MIC_mean+distance_from_endpoint
        arrows(x0 = df_traj_resistance[g,]$time+line_delta_x, y0 = y0-line_delta_y, x1 = df_traj_resistance[g,]$time+line_delta_x, y1 = y1, length = 0, col = combine_color)
        arrows(x0 = df_traj_resistance[g,]$time-line_delta_x, y0 = y0-line_delta_y, x1 = df_traj_resistance[g,]$time-line_delta_x, y1 = y1, length = 0, col = combine_color)      }
    }
    if (df_traj_resistance$time[g] == max(df_traj_resistance$time)) {
      cols<-character(0)
      for(p in 5:(5+N.loci-1)) {
        if(df_traj_resistance[g,1+p]==0) {
          cols<-c(cols,starting.color)
        } else {
          cols<-c(cols,ending.color)
        }
      }
      draw.pie.yaxis(center.x = df_traj_resistance[g,]$time,
                     center.y = df_traj_resistance[g,]$MIC_mean,
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
                     axis.multiplier=axis.multiplier,
                     xy.ratio = xyratio)
      }
  }
  cols<-character(0)
  for(p in 5:(5+N.loci-1)) {
    if(df_traj_resistance[1,1+p]==0) {
      cols<-c(cols,starting.color)
    } else {
      cols<-c(cols,ending.color)
    }
  }
  draw.pie.yaxis(center.x = df_traj_resistance[1,]$time,
                 center.y = df_traj_resistance[1,]$MIC_mean,
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
                 axis.multiplier=axis.multiplier,
                 xy.ratio = xyratio)
}
# HGT_plot(df_filename = 'inputs/aligned.csv',
#          df_trajectories_filename = 'inputs/Trajectory.csv',
#          treatment_filter = 'HGT',
#          x.main.label = 'time',
#          focal_color = red_Ec_color,
#          transient_color = blue_Kp_color,
#          combine_color = purple_Ec.Kp_color,
#          line_thickness = 1,
#          node.edge.width = 0.5,
#          xyratio = 1,
#          radius.buffer = 1,
#          global.cex = global.cex,
#          node.radius = 0.15,
#          arrow_length = 0.075,
#          arrow_proportion = 35,
#          line_proportion_multiplier = 4,
#          buffer_proportion = 50)

###################### ------------
# Arrange Fig1 -------
matrix_layout <- as.matrix(read.csv('inputs/matrix.csv', header = F))
n.row = nrow(matrix_layout)
n.col = ncol(matrix_layout)
pdf.width = 3.25*2
outer_margin = 0.01
square_in = 1.4
L_left_in = 0.2
L_right_in = 0.15
R_left_in = 0
R_right_in = 0
R_column = square_in + R_left_in + R_right_in
L_column = square_in + L_left_in + L_right_in
if (pdf.width < ((R_column*3) + L_column + outer_margin + outer_margin)){
  print('check plot widths')}
top_in = 0
T_top_in = 0.1
bottom_in = 0.4
row = square_in + top_in + bottom_in
pdf.height = (row*3)+T_top_in+outer_margin+outer_margin #max 9.25 inches

# Initiate Fig1 ----
pdf(file = 'Fig1.pdf', width = pdf.width, height = pdf.height)
layout(matrix(matrix_layout, nrow = n.row, ncol = n.col), 
       widths = c(lcm(L_column*2.54),lcm(R_column*2.54),lcm(L_column*2.54),lcm(R_column*2.54)), 
       heights = c(lcm((row+T_top_in)*2.54),lcm(row*2.54),lcm(row*2.54),lcm(row*2.54)))
par(omi=c(outer_margin,outer_margin,outer_margin,outer_margin)) # outer margins (b, l, t, r)
global.cex = 0.75

# landscape parameters -----
node_radius = 0.15
node_edge_width =  0.5
line_thickness = 1
# evolution parameters ------
x_landscape_length = 0.25+3+0.25
x_evolution_length = 0.25+4
radius_multiplier = x_evolution_length/x_landscape_length
arrow_len = 0.075
arrow_prop = 17.5
line_prop_mult = 3
buffer_prop = 50
#### A panel ##### -----
# left ------
par(mai=c(bottom_in,L_left_in,T_top_in,L_right_in))
fitness.landscape.graph(df_filename = 'inputs/aligned.csv',
                        offsets_filename = "inputs/aligned_offsets.csv",
                        host = 'Ec',
                        host.color = red_Ec_color, 
                        node.edge.width = node_edge_width,
                        xyratio = 1,
                        axis.label = T,
                        radius.buffer = 1,
                        node.radius = node_radius,
                        global.cex = global.cex)
mtext("a", side=2, line = 1, at = 11, cex = global.cex, las = 1)
# right ------
par(mai=c(bottom_in,R_left_in,T_top_in,R_right_in))
fitness.landscape.graph(df_filename = 'inputs/aligned.csv',
                        offsets_filename = "inputs/aligned_offsets.csv",
                        host = 'Kp',
                        host.color = blue_Kp_color,
                        node.edge.width = node_edge_width,
                        xyratio = 1,
                        axis.label = F,
                        radius.buffer = 1,
                        node.radius = node_radius,
                        global.cex = global.cex)

#### C panel ##### -----
# left ------
par(mai=c(bottom_in,L_left_in,top_in,L_right_in))
fitness.landscape.graph(df_filename = 'inputs/anticorrelated.csv',
                        offsets_filename = "inputs/aligned_offsets.csv",
                        host = 'Ec',
                        host.color = red_Ec_color, 
                        node.edge.width = node_edge_width,
                        xyratio = 1,
                        axis.label = T,
                        radius.buffer = 1,
                        node.radius = node_radius,
                        global.cex = global.cex)
mtext("c", side=2, line = 1, at = 11, cex = global.cex, las = 1) 
# right ------
par(mai=c(bottom_in,R_left_in,top_in,R_right_in))
fitness.landscape.graph(df_filename = 'inputs/anticorrelated.csv',
                        offsets_filename = "inputs/aligned_offsets.csv",
                        host = 'Kp',
                        host.color = blue_Kp_color,
                        node.edge.width = node_edge_width,
                        xyratio = 1,
                        axis.label = F,
                        radius.buffer = 1,
                        node.radius = node_radius,
                        global.cex = global.cex)
#### E panel ##### -----
# left ------
par(mai=c(bottom_in,L_left_in,top_in,L_right_in))
fitness.landscape.graph(df_filename = 'inputs/suboptimal.csv',
                        offsets_filename = "inputs/aligned_offsets.csv",
                        host = 'Ec',
                        host.color = red_Ec_color, 
                        node.edge.width = node_edge_width,
                        xyratio = 1,
                        axis.label = T,
                        radius.buffer = 1,
                        node.radius = node_radius,
                        global.cex = global.cex)
mtext("e", side=2, line = 1, at = 11, cex = global.cex, las = 1) 
# right  ------
par(mai=c(bottom_in,R_left_in,top_in,R_right_in))
fitness.landscape.graph(df_filename = 'inputs/suboptimal.csv',
                        offsets_filename = "inputs/aligned_offsets.csv",
                        host = 'Kp',
                        host.color = blue_Kp_color,
                        node.edge.width = node_edge_width,
                        xyratio = 1,
                        axis.label = F,
                        radius.buffer = 1,
                        node.radius = node_radius,
                        global.cex = global.cex)
#### B panel ##### -----
# left ------
par(mai=c(bottom_in,L_left_in,T_top_in,L_right_in))
mutation_plot(df_filename = 'inputs/aligned.csv',
              df_trajectories_filename = 'inputs/Trajectory.csv',
              treatment_filter = 'noHGT',
              x.main.label = 'time',
              focal_color = red_Ec_color,
              line_thickness = line_thickness,
              node.edge.width = node_edge_width,
              xyratio = 1,
              radius.buffer = 1,
              global.cex = global.cex,
              node.radius = node_radius * radius_multiplier,
              arrow_length = arrow_len)

mtext("b", side=2, line = 0.5, at = 11, cex = global.cex, las = 1) 
# right ------
par(mai=c(bottom_in,R_left_in,T_top_in,R_right_in))
HGT_plot(df_filename = 'inputs/aligned.csv',
         df_trajectories_filename = 'inputs/Trajectory.csv',
         treatment_filter = 'HGT',
         x.main.label = 'time',
         focal_color = red_Ec_color,
         transient_color = blue_Kp_color,
         combine_color = purple_Ec.Kp_color,
         line_thickness = line_thickness,
         node.edge.width = node_edge_width,
         xyratio = 1,
         radius.buffer = 1,
         global.cex = global.cex,
         node.radius = node_radius * radius_multiplier,
         arrow_length = arrow_len,
         arrow_proportion = arrow_prop,
         line_proportion_multiplier = line_prop_mult,
         buffer_proportion = buffer_prop)
#### D panel ##### -----
# left ------
par(mai=c(bottom_in,L_left_in,top_in,L_right_in))
mutation_plot(df_filename = 'inputs/anticorrelated.csv',
              df_trajectories_filename = 'inputs/Trajectory.csv',
              treatment_filter = 'noHGT',
              x.main.label = 'time',
              focal_color = red_Ec_color,
              line_thickness = line_thickness,
              node.edge.width = node_edge_width,
              xyratio = 1,
              radius.buffer = 1,
              global.cex = global.cex,
              node.radius = node_radius * radius_multiplier,
              arrow_length = arrow_len)

mtext("d", side=2, line = 0.5, at = 11, cex = global.cex, las = 1) 
# right ------
par(mai=c(bottom_in,R_left_in,top_in,R_right_in))
HGT_plot(df_filename = 'inputs/anticorrelated.csv',
         df_trajectories_filename = 'inputs/Trajectory_anticorrelated.csv',
         treatment_filter = 'HGT',
         x.main.label = 'time',
         focal_color = red_Ec_color,
         transient_color = blue_Kp_color,
         combine_color = purple_Ec.Kp_color,
         line_thickness = line_thickness,
         node.edge.width = node_edge_width,
         xyratio = 1,
         radius.buffer = 1,
         global.cex = global.cex,
         node.radius = node_radius * radius_multiplier,
         arrow_length = arrow_len,
         arrow_proportion = arrow_prop,
         line_proportion_multiplier = line_prop_mult,
         buffer_proportion = buffer_prop)
#### F panel ##### -----
# left ------
par(mai=c(bottom_in,L_left_in,top_in,L_right_in))
mutation_plot(df_filename = 'inputs/suboptimal.csv',
              df_trajectories_filename = 'inputs/Trajectory.csv',
              treatment_filter = 'noHGT',
              x.main.label = 'time',
              focal_color = red_Ec_color,
              line_thickness = line_thickness,
              node.edge.width = node_edge_width,
              xyratio = 1,
              radius.buffer = 1,
              global.cex = global.cex,
              node.radius = node_radius * radius_multiplier,
              arrow_length = arrow_len)

mtext("f", side=2, line = 0.5, at = 11, cex = global.cex, las = 1) 
# right  ------
par(mai=c(bottom_in,R_left_in,top_in,R_right_in))
HGT_plot(df_filename = 'inputs/suboptimal.csv',
         df_trajectories_filename = 'inputs/Trajectory.csv',
         treatment_filter = 'HGT',
         x.main.label = 'time',
         focal_color = red_Ec_color,
         transient_color = blue_Kp_color,
         combine_color = purple_Ec.Kp_color,
         line_thickness = line_thickness,
         node.edge.width = node_edge_width,
         xyratio = 1,
         radius.buffer = 1,
         global.cex = global.cex,
         node.radius = node_radius * radius_multiplier,
         arrow_length = arrow_len,
         arrow_proportion = arrow_prop,
         line_proportion_multiplier = line_prop_mult,
         buffer_proportion = buffer_prop)
dev.off()

