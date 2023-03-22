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

# 4 parameter model fitting and plotting  -------------------------------------------------------
drm.func <- function(x) {
  drm(GrowthRate ~ Concentration, 
      fct = L.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                 fixed = c(1, 0, 10, NA)), 
      data = x)
}
predict.fun <- function(x) {
  add_predictions(data.frame(Concentration = c(seq(-3, 11, 0.01))), x)}

# Build dose-response data frame---------
Data_for_fit <- data.frame('hosts' = NA, 'Genotype' = NA, 'Concentration' = NA, 'GrowthRate' = NA)
df_reformat <- read.csv('inputs/suboptimal.csv')
for(g in 1:length(df_reformat$Species)) {
  IC50 = df_reformat[g,]$MIC_mean
  df_g <- data.frame('hosts' = rep(df_reformat[g,]$Species, 3), 
                     'Genotype' = rep(df_reformat[g,]$Genotype, 3), 
                     'Concentration' = c(IC50-1, IC50, IC50+1), 
                     'GrowthRate' = c(10,5,0))
  Data_for_fit <- rbind(Data_for_fit, df_g)
}
Data_for_fit <- Data_for_fit[2:length(Data_for_fit$hosts),]

# Model fitting ----------
Data_fit <- Data_for_fit %>% 
  group_by(hosts, Genotype) %>% 
  nest() %>%
  mutate(drmod = map(data, drm.func), 
         pred = map(drmod, predict.fun))

# Plot dose-response -------
Dose_response_plot <- function(df_reformat,
                               Data_fit,
                               host,
                               focal_color,
                               node.radius, 
                               node.edge.width,
                               axis.label = F, 
                               global.cex,
                               xyratio){
  max_concentration = 15
  N.loci = 3
  y_min = -2.25
  y_max = 10
  x_min = -3
  x_max = 11
  x_label = 'drug concentration'
  y_label = 'growth rate'
  pie_y_center = -1.5
  node.edge.color = 'black'
  node.edge.width = node.edge.width
  starting.color = 'white'
  ending.color = pie_color
  axis.multiplier = (y_max-(-abs(y_min)))/(x_max + abs(x_min))
  label.pos = 0.1
  node.label.size = 0.1
  color_genotypes = c('WT', 'A', 'B')
  df<-df_reformat %>%
    rename(hosts = Species,
           fit.1 = MIC_mean) %>%
    dplyr::select(hosts, A, B, C, fit.1, MIC_se, Num_of_mutations, Genotype)
  df.rel <- df %>%
    group_by(hosts) %>%
    mutate(MIC_WT = fit.1[Num_of_mutations == 0],
           fit.1 = fit.1)
  df.host <- df.rel %>% filter(hosts == host)
  fit.host <- Data_fit %>% filter(hosts == host)
  plot(seq(x_min,x_max,length.out = 2), seq(y_min,y_max,length.out = 2), type="n", xaxt='n', xlab='', ylab='', axes=FALSE, frame.plot=F, xaxs="i")
  for(g in 1:length(df.host$hosts)) {
    focal_g = df.host[g,]$Genotype
    if (focal_g %in% color_genotypes){
      lines(fit.host$pred[fit.host$Genotype == focal_g][[1]]$Concentration, fit.host$pred[fit.host$Genotype == focal_g][[1]]$pred, col = focal_color) # colored curves
      arrows(x0 = fit.host$data[fit.host$Genotype == focal_g][[1]]$Concentration[2], y0 = 0, x1 = fit.host$data[fit.host$Genotype == focal_g][[1]]$Concentration[2], y1 = y_max/2, length = 0, lty = 2, col = focal_color)
      points(fit.host$data[fit.host$Genotype == focal_g][[1]]$Concentration[2], y_max/2, pch = 16, col=focal_color, cex = 0.5, lwd=1)
      cols<-character(0)
      for(p in 1:N.loci) {
        if(df.host[g,1+p]==0) {
          cols<-c(cols,starting.color)
        } else {
          cols<-c(cols,ending.color)
        }
      }
      draw.pie.yaxis(center.x = df.host[g,]$fit.1,
                     center.y = pie_y_center,
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
    } else {
      lines(fit.host$pred[fit.host$Genotype == focal_g][[1]]$Concentration, fit.host$pred[fit.host$Genotype == focal_g][[1]]$pred, col = pie_color) # grey curves
    }
  }
  if (axis.label == T) {
    mtext(side = 2, line = 0.5, text = y_label, cex = global.cex)
  }
  axis(side = 1, pos = 0, labels = NA, lwd.ticks = 0, at = c(x_min,x_max), cex = global.cex)
  axis(side = 3, pos = y_max+0.3, labels = NA, lwd.ticks = 0, at = c(x_min,x_max), cex = global.cex)
  axis(side = 4, pos = x_max, labels = NA, lwd.ticks = 0, at = c(0,y_max+0.3), cex = global.cex)
  axis(side = 2, pos = x_min, labels = NA, lwd.ticks = 0, at = c(0,y_max+0.3), cex = global.cex)
  mtext(side = 1, line = 0.25, text = x_label, cex = global.cex)
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

# Plasmid functions ---------
convert.pol.to.car<-function(r,a) {
  x<-r*cos(a)
  y<-r*sin(a)
  return(c(x,y))
}

convert.pol.to.ang<-function(p) {
  return((360*p)/(2*pi))
}


draw.plasmid.lab<-function(x, y, r, w, r.l,
                       genes.pos, genes.col,
                       genes.lab, genes.lab.col,
                       genes.lab.cex,
                       genes.lwd,
                       genes.lab.sp,
                       genes.lab.family="Helvetica",
                       genes.lab.face=2,
                       backbone.col="black",
                       backbone.width=1,
                       res=100) {
  cx<-convert.pol.to.car(r,0)[1]+x
  cy<-convert.pol.to.car(r,0)[2]+y
  for(z in 1:res) {
    a<-0+z*(2*pi)/res
    cx<-c(cx,convert.pol.to.car(r,a)[1]+x)
    cy<-c(cy,convert.pol.to.car(r,a)[2]+y)
  }
  polygon(cx,cy,col='white',border=backbone.col,lwd=backbone.width)
  for(i in 1:length(genes.lab)) {
    ox<-convert.pol.to.car(r+w/2,genes.pos[[i]][1])[1]+x
    oy<-convert.pol.to.car(r+w/2,genes.pos[[i]][1])[2]+y
    for(j in 1:res) {
      a<-genes.pos[[i]][1]-j*(genes.pos[[i]][1]-genes.pos[[i]][2])/res
      ox<-c(ox,convert.pol.to.car(r+w/2,a)[1]+x)
      oy<-c(oy,convert.pol.to.car(r+w/2,a)[2]+y)
    }
    ix<-convert.pol.to.car(r-w/2,genes.pos[[i]][2])[1]+x
    iy<-convert.pol.to.car(r-w/2,genes.pos[[i]][2])[2]+y
    for(j in 1:res) {
      a<-genes.pos[[i]][2]+j*(genes.pos[[i]][1]-genes.pos[[i]][2])/res
      ix<-c(ix,convert.pol.to.car(r-w/2,a)[1]+x)
      iy<-c(iy,convert.pol.to.car(r-w/2,a)[2]+y)
    }
    polygon(c(ox,ix),c(oy,iy),col=genes.col[i], lwd = genes.lwd[i])
    name<-strsplit(genes.lab[i],"")[[1]]
    le<-length(name)
    ce<-(genes.pos[[i]][1]+genes.pos[[i]][2])/2.0
    if(le>0) {
      for(k in 1:le) {
        a <- (ce+(le/2)*genes.lab.sp)-(k-1)*((le)/(le-1))*genes.lab.sp
        text(convert.pol.to.car(r.l,a)[1]+x,convert.pol.to.car(r.l,a)[2]+y,
             name[k], srt=(-90+convert.pol.to.ang(a)), adj=c(0.5,0.5), 
             col=genes.lab.col[i], 
             font=genes.lab.face, family=genes.lab.family,
             cex=genes.lab.cex[i])
      }
    }
  } 
}
draw.plasmid<-function(x, 
                       y, 
                       r, 
                       w,
                       genes.pos, 
                       genes.col,
                       genes.lwd,
                       backbone.col="black",
                       backbone.width=1,
                       res=100) {
  cx<-convert.pol.to.car(r,0)[1]+x
  cy<-convert.pol.to.car(r,0)[2]+y
  for(z in 1:res) {
    a<-0+z*(2*pi)/res
    cx<-c(cx,convert.pol.to.car(r,a)[1]+x)
    cy<-c(cy,convert.pol.to.car(r,a)[2]+y)
  }
  polygon(cx,cy,col='white',border=backbone.col,lwd=backbone.width)
  for(i in 1:length(genes.pos)) {
    ox<-convert.pol.to.car(r+w/2,genes.pos[[i]][1])[1]+x
    oy<-convert.pol.to.car(r+w/2,genes.pos[[i]][1])[2]+y
    for(j in 1:res) {
      a<-genes.pos[[i]][1]-j*(genes.pos[[i]][1]-genes.pos[[i]][2])/res
      ox<-c(ox,convert.pol.to.car(r+w/2,a)[1]+x)
      oy<-c(oy,convert.pol.to.car(r+w/2,a)[2]+y)
    }
    ix<-convert.pol.to.car(r-w/2,genes.pos[[i]][2])[1]+x
    iy<-convert.pol.to.car(r-w/2,genes.pos[[i]][2])[2]+y
    for(j in 1:res) {
      a<-genes.pos[[i]][2]+j*(genes.pos[[i]][1]-genes.pos[[i]][2])/res
      ix<-c(ix,convert.pol.to.car(r-w/2,a)[1]+x)
      iy<-c(iy,convert.pol.to.car(r-w/2,a)[2]+y)
    }
    polygon(c(ox,ix),c(oy,iy),col=genes.col[i], lwd = genes.lwd[i])
  }
}
draw.plasmid.back<-function(x, 
                       y, 
                       r, 
                       w,
                       genes.pos, 
                       genes.col,
                       genes.lwd,
                       backbone.col="black",
                       border.col,
                       backbone.width=1,
                       res=100) {
  cx<-convert.pol.to.car(r,0)[1]+x
  cy<-convert.pol.to.car(r,0)[2]+y
  for(z in 1:res) {
    a<-0+z*(2*pi)/res
    cx<-c(cx,convert.pol.to.car(r,a)[1]+x)
    cy<-c(cy,convert.pol.to.car(r,a)[2]+y)
  }
  polygon(cx,cy,col='white',border=backbone.col,lwd=backbone.width)
  for(i in 1:length(genes.pos)) {
    ox<-convert.pol.to.car(r+w/2,genes.pos[[i]][1])[1]+x
    oy<-convert.pol.to.car(r+w/2,genes.pos[[i]][1])[2]+y
    for(j in 1:res) {
      a<-genes.pos[[i]][1]-j*(genes.pos[[i]][1]-genes.pos[[i]][2])/res
      ox<-c(ox,convert.pol.to.car(r+w/2,a)[1]+x)
      oy<-c(oy,convert.pol.to.car(r+w/2,a)[2]+y)
    }
    ix<-convert.pol.to.car(r-w/2,genes.pos[[i]][2])[1]+x
    iy<-convert.pol.to.car(r-w/2,genes.pos[[i]][2])[2]+y
    for(j in 1:res) {
      a<-genes.pos[[i]][2]+j*(genes.pos[[i]][1]-genes.pos[[i]][2])/res
      ix<-c(ix,convert.pol.to.car(r-w/2,a)[1]+x)
      iy<-c(iy,convert.pol.to.car(r-w/2,a)[2]+y)
    }
    polygon(c(ox,ix),c(oy,iy),col=genes.col[i], border = border.col[i], lwd = genes.lwd[i])
  }
}
draw.segment<-function(x, y, r, w,
                       genes.pos, 
                       genes.rip,
                       rip.w, 
                       rip.n, 
                       genes.col,
                       genes.lwd, 
                       backbone.col="black",
                       border.col,
                       backbone.width=1) {
  lines(c(x-r,x+r), c(y,y), col=backbone.col, lwd=backbone.width)
  for(i in 1:length(genes.pos)) {
    rect(genes.pos[[i]][1],y-w,genes.pos[[i]][2],y+w,col=genes.col[i], border = border.col[i], lwd = genes.lwd[i])
  }
}

# Tube functions --------
draw.tube<-function(center,width,height,liquid.level,side.view,resolution=100,angle=0,liquid="lightblue",meniscus=c("black",2),glass=c("black",2)) {
  rad<-width/2
  bottom.x<-numeric(0)
  bottom.y<-numeric(0)
  rad.seq<-rev(seq(-pi,0,pi/resolution))
  for(r in rad.seq) {
    bottom.x<-c(bottom.x,rad*cos(r))
    bottom.y<-c(bottom.y,rad*sin(r))
  }
  liq.x<-c(rad,rad,bottom.x,-rad,-rad)
  liq.y<-c((liquid.level-0.5)*height,(-0.5)*height,
           bottom.y-(0.5)*height,(-0.5)*height,
           (liquid.level-0.5)*height)
  rotate.x<-liq.x*cos(angle) - liq.y*sin(angle)
  rotate.y<-liq.x*sin(angle) + liq.y*cos(angle)
  liquid.x<-center[1]+rotate.x
  liquid.y<-center[2]+rotate.y
  polygon(liquid.x, liquid.y, col=liquid, border=NA)
  rad.seq<-seq(0,2*pi,pi/resolution)
  meni.x<-numeric(0)
  meni.y<-numeric(0)
  a<-rad
  b<-(1-side.view)*rad
  for(r in rad.seq) {
    meni.r<-(a*b)/(sqrt((a*sin(r))^2 +(b*cos(r))^2)) 
    meni.x<-c(meni.x,meni.r*cos(r))
    meni.y<-c(meni.y,meni.r*sin(r))
  }
  base.x<-meni.x
  base.y<-meni.y+(liquid.level-0.5)*height
  rotate.x<-base.x*cos(angle) - base.y*sin(angle)
  rotate.y<-base.x*sin(angle) + base.y*cos(angle)
  meniscus.x<-center[1]+rotate.x
  meniscus.y<-center[2]+rotate.y
  polygon(meniscus.x, meniscus.y, col=liquid, 
          border=meniscus[1],lwd=meniscus[2])
  meni.half.x<-numeric(0)
  meni.half.y<-numeric(0)
  rad.seq.half<-rev(seq(0,pi,pi/resolution))
  for(r in rad.seq.half) {
    meni.r<-(a*b)/(sqrt((a*sin(r))^2 +(b*cos(r))^2)) 
    meni.half.x<-c(meni.half.x,meni.r*cos(r))
    meni.half.y<-c(meni.half.y,meni.r*sin(r))
  }
  base.x<-c(meni.x,rad,rad,bottom.x,-rad,-rad,meni.half.x)
  base.y<-c(height/2+meni.y,height/2,-height/2,
            bottom.y-(0.5)*height,-height/2,
            height/2,height/2+meni.half.y)
  rotate.x<-base.x*cos(angle) - base.y*sin(angle)
  rotate.y<-base.x*sin(angle) + base.y*cos(angle)
  glass.x<-center[1]+rotate.x
  glass.y<-center[2]+rotate.y
  polygon(glass.x, glass.y, col=NA, 
          border=glass[1],lwd=glass[2])
}
draw.petri.dish<-function(center.x,
                          center.y,
                          radius,
                          agar.color,
                          edge.color,
                          edge.width,
                          N.colonies,
                          colony.radius,
                          colony.color,
                          critical.spacing) {
  angles<-seq(0,2*pi,0.01)
  polygon(center.x+radius*cos(angles),center.y+radius*sin(angles),
          col=agar.color,border=NA)
  lines(center.x+radius*cos(angles),center.y+radius*sin(angles),
        col=edge.color,lwd=edge.width)
  x.list<-numeric(0)
  y.list<-numeric(0)
  for(c in 1:N.colonies) {
    is.looking.for.colony<-TRUE
    while(is.looking.for.colony) {
      is.isolated <- TRUE
      random.angle <- runif(1,min=0,max=2*pi)
      random.radius <- runif(1,min=0,max=radius-2*colony.radius)
      delta.x <- random.radius*cos(random.angle)
      delta.y <- random.radius*sin(random.angle)
      if(c>=2) {
        for(oc in 1:(c-1)) {
          inter.colony.distace <- distance(center.x+delta.x,
                                           center.y+delta.y,
                                           x.list[oc],
                                           y.list[oc])
          if(inter.colony.distace < critical.spacing) {
            is.isolated <- FALSE
          }
        }  
      }
      if(is.isolated) {
        is.looking.for.colony<-FALSE
        polygon(center.x+delta.x+colony.radius*cos(angles),
                center.y+delta.y+colony.radius*sin(angles),
                col=colony.color,border=NA)
        x.list <- c(x.list,center.x+delta.x)
        y.list <- c(y.list,center.y+delta.y)
      }
    }
  }
}
draw.cell<-function(center,length,width,rounded,resolution=100,angle=0,cytoplasm="white",membrane=c("black",2),...) {
  short.rad<-width/2
  long.rad<-length/2
  right.cap.x<-numeric(0)
  right.cap.y<-numeric(0)
  left.cap.x<-numeric(0)
  left.cap.x<-numeric(0)
  rad.seq<-rev(seq(-pi/2,pi/2,pi/resolution))
  for(r in rad.seq) {
    if(r >= pi/4) {
      right.cap.x<-c(right.cap.x,(rounded*(short.rad)+(1-rounded)*(short.rad/cos((pi/2)-r)))*cos(r))
      right.cap.y<-c(right.cap.y,(rounded*(short.rad)+(1-rounded)*(short.rad/cos((pi/2)-r)))*sin(r))
    } else if(r>=-pi/4) {
      right.cap.x<-c(right.cap.x,(rounded*(short.rad)+(1-rounded)*(short.rad/cos(r)))*cos(r))
      right.cap.y<-c(right.cap.y,(rounded*(short.rad)+(1-rounded)*(short.rad/cos(r)))*sin(r))
    } else {
      right.cap.x<-c(right.cap.x,(rounded*(short.rad)+(1-rounded)*(short.rad/cos((pi/2)+r)))*cos(r))
      right.cap.y<-c(right.cap.y,(rounded*(short.rad)+(1-rounded)*(short.rad/cos((pi/2)+r)))*sin(r))
    }
  }
  left.cap.x<-rev(-right.cap.x)
  left.cap.y<-rev(right.cap.y)
  base.x<-c(-(long.rad-short.rad),(long.rad-short.rad),
            right.cap.x+(long.rad-short.rad),(long.rad-short.rad),
            -(long.rad-short.rad),left.cap.x-(long.rad-short.rad))
  base.y<-c(short.rad,short.rad,right.cap.y,-short.rad,-short.rad,left.cap.y)
  rotate.x<-base.x*cos(angle) - base.y*sin(angle)
  rotate.y<-base.x*sin(angle) + base.y*cos(angle)
  cell.x<-center[1]+rotate.x
  cell.y<-center[2]+rotate.y
  polygon(cell.x, cell.y, col=cytoplasm, border=membrane[1], lwd=membrane[2])
  plasmids<-list(...)
  Nplasmids<-length(plasmids)
  if(length(plasmids)>0) {
    for(p in 1:Nplasmids) {
      rad.seq<-seq(0,2*pi,pi/resolution)
      p.rad<-short.rad*plasmids[[p]][3]
      base.p.x<-numeric(0)
      base.p.y<-numeric(0)
      for(r in rad.seq) {
        base.p.x<-c(base.p.x,p.rad*cos(r)) 
        base.p.y<-c(base.p.y,p.rad*sin(r)) 
      }
      shift.p.x<-long.rad*plasmids[[p]][1]+base.p.x
      shift.p.y<-short.rad*plasmids[[p]][2]+base.p.y
      rotate.p.x<-shift.p.x*cos(angle) - shift.p.y*sin(angle)
      rotate.p.y<-shift.p.x*sin(angle) + shift.p.y*cos(angle)
      plasmid.x<-center[1]+rotate.p.x
      plasmid.y<-center[2]+rotate.p.y
      polygon(plasmid.x, plasmid.y, col=NA, border=rgb(plasmids[[p]][4],plasmids[[p]][5],plasmids[[p]][6]), lwd=plasmids[[p]][7])
    }
  }
}
# Dish functions -------
getDistance <- function(point.a, point.b, resolution=1){
                        distance.ab <- 0
                        #get distances (grid)
                        dx <- point.a[1]-point.b[1]
                        dy <- point.a[2]-point.b[2]
                        #get distances (degree)
                        xdist <- dx * resolution
                        ydist <- dy * resolution
                        #pythagoras
                        distance.ab <- sqrt(xdist^2 + ydist^2)
                        return(distance.ab)}
draw.dish<-function(center,
                    width,
                    height,
                    liquid.level,
                    side.view,
                    resolution=100,
                    angle=0,
                    liquid="lightblue",
                    meniscus=c("black",2),
                    glass=c("black",2)) {
  rad<-width/2
  rad.seq<-rev(seq(-pi,0,pi/resolution)) # bottom of the tube
  bottom.x<-numeric(0)
  bottom.y<-numeric(0)
  a<-rad
  b<-(1-side.view)*rad
  for(r in rad.seq) {
    bottom.r<-(a*b)/(sqrt((a*sin(r))^2 +(b*cos(r))^2)) 
    bottom.x<-c(bottom.x,bottom.r*cos(r))
    bottom.y<-c(bottom.y,bottom.r*sin(r))}
  liq.x<-c(rad,
           rad,
           bottom.x,
           -rad,
           -rad)
  liq.y<-c((liquid.level-0.5)*height,
           (-0.5)*height,
           bottom.y-(0.5)*height,
           (-0.5)*height,
           (liquid.level-0.5)*height)
  rotate.x<-liq.x*cos(angle) - liq.y*sin(angle)
  rotate.y<-liq.x*sin(angle) + liq.y*cos(angle)
  liquid.x<-center[1]+rotate.x
  liquid.y<-center[2]+rotate.y
  polygon(liquid.x, liquid.y, col=liquid, border=NA)  
  rad.seq<-seq(0,2*pi,pi/resolution)
  meni.x<-numeric(0)
  meni.y<-numeric(0)
  a<-rad
  b<-(1-side.view)*rad
  for(r in rad.seq) {
    meni.r<-(a*b)/(sqrt((a*sin(r))^2 +(b*cos(r))^2)) 
    meni.x<-c(meni.x,meni.r*cos(r))
    meni.y<-c(meni.y,meni.r*sin(r))}
  base.x<-meni.x
  base.y<-meni.y+(liquid.level-0.5)*height
  rotate.x<-base.x*cos(angle) - base.y*sin(angle)
  rotate.y<-base.x*sin(angle) + base.y*cos(angle)
  meniscus.x<-center[1]+rotate.x
  meniscus.y<-center[2]+rotate.y
  polygon(meniscus.x, meniscus.y, col=liquid, border=meniscus[1],lwd=meniscus[2])
  meni.half.x<-numeric(0)
  meni.half.y<-numeric(0)
  rad.seq.half<-rev(seq(0,pi,pi/resolution))
  for(r in rad.seq.half) {
    meni.r<-(a*b)/(sqrt((a*sin(r))^2 +(b*cos(r))^2)) 
    meni.half.x<-c(meni.half.x,meni.r*cos(r))
    meni.half.y<-c(meni.half.y,meni.r*sin(r))}
  base.x<-c(meni.x,rad,rad,bottom.x,-rad,-rad,meni.half.x)
  base.y<-c(height/2+meni.y,height/2,-height/2,
            bottom.y-(0.5)*height,-height/2,
            height/2,height/2+meni.half.y)
  rotate.x<-base.x*cos(angle) - base.y*sin(angle)
  rotate.y<-base.x*sin(angle) + base.y*cos(angle)
  glass.x<-center[1]+rotate.x
  glass.y<-center[2]+rotate.y
  polygon(glass.x, glass.y, col=NA, border=glass[1],lwd=glass[2])
}

draw.dish.colonies<-function(center,
                    width,
                    height,
                    liquid.level,
                    side.view,
                    N.colonies,
                    colony.radius,
                    colony.color,
                    critical.spacing,
                    edge_spacing,
                    resolution=100,
                    angle=0,
                    liquid="lightblue",
                    meniscus=c("black",2),
                    glass=c("black",2)) {
  rad<-width/2
  ### Fix this to be the miniscus ###
  rad.seq<-rev(seq(-pi,0,pi/resolution)) # bottom of the tube
  bottom.x<-numeric(0)
  bottom.y<-numeric(0)
  a<-rad
  b<-(1-side.view)*rad
  for(r in rad.seq) {
    bottom.r<-(a*b)/(sqrt((a*sin(r))^2 +(b*cos(r))^2)) 
    bottom.x<-c(bottom.x,bottom.r*cos(r))
    bottom.y<-c(bottom.y,bottom.r*sin(r))}
  liq.x<-c(rad,
           rad,
           bottom.x,
           -rad,
           -rad)
  liq.y<-c((liquid.level-0.5)*height,
           (-0.5)*height,
           bottom.y-(0.5)*height,
           (-0.5)*height,
           (liquid.level-0.5)*height)
  rotate.x<-liq.x*cos(angle) - liq.y*sin(angle)
  rotate.y<-liq.x*sin(angle) + liq.y*cos(angle)
  liquid.x<-center[1]+rotate.x
  liquid.y<-center[2]+rotate.y
  polygon(liquid.x, liquid.y, col=liquid, border=NA)  
  rad.seq<-seq(0,2*pi,pi/resolution)
  meni.x<-numeric(0)
  meni.y<-numeric(0)
  a<-rad
  b<-(1-side.view)*rad
  for(r in rad.seq) {
    meni.r<-(a*b)/(sqrt((a*sin(r))^2 +(b*cos(r))^2)) 
    meni.x<-c(meni.x,meni.r*cos(r))
    meni.y<-c(meni.y,meni.r*sin(r))}
  base.x<-meni.x
  base.y<-meni.y+(liquid.level-0.5)*height
  rotate.x<-base.x*cos(angle) - base.y*sin(angle)
  rotate.y<-base.x*sin(angle) + base.y*cos(angle)
  meniscus.x<-center[1]+rotate.x
  meniscus.y<-center[2]+rotate.y
  polygon(meniscus.x, meniscus.y, col=liquid, border=meniscus[1],lwd=meniscus[2])
  meni.half.x<-numeric(0)
  meni.half.y<-numeric(0)
  rad.seq.half<-rev(seq(0,pi,pi/resolution))
  for(r in rad.seq.half) {
    meni.r<-(a*b)/(sqrt((a*sin(r))^2 +(b*cos(r))^2)) 
    meni.half.x<-c(meni.half.x,meni.r*cos(r))
    meni.half.y<-c(meni.half.y,meni.r*sin(r))}
  base.x<-c(meni.x,rad,rad,bottom.x,-rad,-rad,meni.half.x)
  base.y<-c(height/2+meni.y,height/2,-height/2,
            bottom.y-(0.5)*height,-height/2,
            height/2,height/2+meni.half.y)
  rotate.x<-base.x*cos(angle) - base.y*sin(angle)
  rotate.y<-base.x*sin(angle) + base.y*cos(angle)
  glass.x<-center[1]+rotate.x
  glass.y<-center[2]+rotate.y
  x.list<-numeric(0)
  y.list<-numeric(0)
  for(c in 1:N.colonies) {
    is.looking.for.colony<-TRUE
    while(is.looking.for.colony) {
      is.isolated <- TRUE
      random.angle <- runif(1,min=0,max=2*pi)
      random.radius <- runif(1,min=0,max=rad-edge_spacing)
      a <- random.radius
      b <- (1-side.view)*random.radius
      delta.x <- a*cos(random.angle)
      delta.y <- b*sin(random.angle)
      if(c>=2) {
        for(oc in 1:(c-1)) {
          inter.colony.distace <- getDistance(c(center[1]+delta.x,center[2]+delta.y),
                                              c(x.list[oc],y.list[oc]),resolution = 1)
          if(inter.colony.distace < critical.spacing) {
            is.isolated <- FALSE
          }
        }  
      }
      if(is.isolated) {
        is.looking.for.colony<-FALSE
        polygon(center[1]+delta.x+colony.radius*meni.x,
                center[2]+delta.y+colony.radius*meni.y,
                col=colony.color,border=NA)
        x.list <- c(x.list,center[1]+delta.x)
        y.list <- c(y.list,center[2]+delta.y)
      }
    }
  }
  polygon(glass.x, glass.y, col=NA, border=glass[1],lwd=glass[2])
}

###################### ------------
# Arrange Fig2 -------
matrix_layout <- as.matrix(read.csv('inputs/matrix.csv', header = F))
n.row = nrow(matrix_layout)
n.col = ncol(matrix_layout)
pdf.width = 3.25
outer_margin = 0.01
T_square_in = 0.8
B_square_in = 1.4
L_left_in = 0.2
L_right_in = 0.15
R_left_in = 0
R_right_in = 0
R_column = B_square_in + R_left_in + R_right_in
L_column = B_square_in + L_left_in + L_right_in
if (pdf.width < (R_column + L_column + outer_margin + outer_margin)){
  print('check plot widths')}
T_top_in = 0
T_bottom_in = 0.3
B_top_in = 0
B_bottom_in = 0.4
A_row = pdf.width + (pdf.width*0.8)
T_row = T_square_in + T_top_in + T_bottom_in
B_row = B_square_in + B_top_in + B_bottom_in
C_D_bracket = 9.25 - A_row - T_row - B_row - (outer_margin*2)
pdf.height = 9.25 #max 9.25 inches

# Initiate Fig2 ----
pdf(file = 'Fig2.pdf', width = pdf.width, height = pdf.height)
layout(matrix(matrix_layout, nrow = n.row, ncol = n.col), 
       widths = c(lcm(L_column*2.54),lcm(R_column*2.54)), 
       heights = c(A_row,lcm(T_row*2.54),C_D_bracket,lcm(B_row*2.54)))
par(omi=c(outer_margin,outer_margin,outer_margin,outer_margin)) # outer margins (b, l, t, r)
global.cex = 0.75
###### A panel ##### -----
## plasmid row -----
par(mai = c(0,0,0,0))
plot(c(0,1),c(-0.8,1), type = 'n', axes = F)
x_genotype_plasmid = c(0.1, 0.3, 0.5, 0.9)
y_genotype_plasmid = 0.965
plasmid_radius = 0.05
plasmid_gene_width = plasmid_radius/2.4
gene_location <- c(3,1)
A_mutation_location <- c(2.7,2.5)
B_mutation_location <- c(2.1,1.9)
C_mutation_location <- c(1.5,1.3)
gene_thickness <- 0.5
mutation_thickness <- 0.5
barcode_thickness <- 0.5
cell_width = 0.1
cell_length = 0.15
cell_spacing = 0.055
cell_thickness = 1
#  plot plasmid row -----
draw.plasmid(x=x_genotype_plasmid[1], y=y_genotype_plasmid, r=plasmid_radius, w=plasmid_gene_width,
             genes.pos=list(gene_location),
             genes.lwd=c(gene_thickness),
             genes.col=c(gene_color), backbone.col="black", backbone.width=1, res=100)
draw.plasmid(x=x_genotype_plasmid[2], y=y_genotype_plasmid, r=plasmid_radius, w=plasmid_gene_width,
             genes.pos=list(gene_location, A_mutation_location),
             genes.lwd=c(gene_thickness, mutation_thickness),
             genes.col=c(gene_color,pie_color), backbone.col="black", backbone.width=1, res=100)
draw.plasmid(x=x_genotype_plasmid[3], y=y_genotype_plasmid, r=plasmid_radius, w=plasmid_gene_width,
             genes.pos=list(gene_location, B_mutation_location),
             genes.lwd=c(gene_thickness, mutation_thickness),
             genes.col=c(gene_color,pie_color), backbone.col="black", backbone.width=1, res=100)
points(x = c(0.65,0.7,0.75), y = rep(y_genotype_plasmid, 3), pch = 16, lwd=0.5)
draw.plasmid(x=x_genotype_plasmid[4], y=y_genotype_plasmid, r=plasmid_radius, w=plasmid_gene_width,
             genes.pos=list(gene_location, A_mutation_location, B_mutation_location, C_mutation_location),
             genes.lwd=c(gene_thickness, rep(mutation_thickness,3)),
             genes.col=c(gene_color, rep(pie_color, 3)), backbone.col="black", backbone.width=1, res=100)

## barcode row -----
y_barcode_center <- 0.69
y_barcode_spacing <- 0.06
y_barcodes <- c(y_barcode_center+y_barcode_spacing, y_barcode_center, y_barcode_center-y_barcode_spacing)
x_barcodes <- x_genotype_plasmid
barcode_white_space <- c(0.65,0.15)
barcode_background_color <- 'white'
barcode_border_color <- 'white'
barcode_border_thickness <- 2
WT_barcodes <- list(c(0.06, 0.30, 0.60, 0.70, 1.00), 
                    c(0.15, 0.20, 0.50, 0.65, 0.91),
                    c(0.06, 0.46, 0.71, 0.88, 0.94))
A_barcodes <- list(c(0.18, 0.28, 0.34, 0.58, 0.88), 
                   c(0.02, 0.11, 0.26, 0.67, 0.98),
                   c(0.27, 0.44, 0.76, 0.85, 0.92))
B_barcodes <- list(c(0.00, 0.30, 0.60, 0.70, 0.98), 
                   c(0.21, 0.28, 0.48, 0.72, 0.89),
                   c(0.06, 0.38, 0.51, 0.78, 0.91))
ABC_barcodes <- list(c(0.12, 0.22, 0.42, 0.62, 0.82), 
                     c(0.02, 0.35, 0.55, 0.75, 0.95),
                     c(0.10, 0.33, 0.43, 0.53, 0.77))
WT_barcodes_locations <- list()
A_barcodes_locations <- list()
B_barcodes_locations <- list()
ABC_barcodes_locations <- list()
for (b in 1:length(WT_barcodes)){
  WT_barcodes_locations <- c(WT_barcodes_locations, list(WT_barcodes[[b]] * (barcode_white_space[1]-barcode_white_space[2]) + barcode_white_space[2]))
  A_barcodes_locations <- c(A_barcodes_locations, list(A_barcodes[[b]] * (barcode_white_space[1]-barcode_white_space[2]) + barcode_white_space[2]))
  B_barcodes_locations <- c(B_barcodes_locations, list(B_barcodes[[b]] * (barcode_white_space[1]-barcode_white_space[2]) + barcode_white_space[2]))
  ABC_barcodes_locations <- c(ABC_barcodes_locations, list(ABC_barcodes[[b]] * (barcode_white_space[1]-barcode_white_space[2]) + barcode_white_space[2]))
}

#  plot barcode row -----
for (y in c(1:length(y_barcodes))){
  locations <- list(gene_location, barcode_white_space) # all of the fragments that are not in the barcode
  for (l in WT_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_barcodes[1], y=y_barcodes[y], r=plasmid_radius, w=plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(gene_thickness, barcode_border_thickness, rep(barcode_thickness, length(WT_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep('black', length(WT_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_barcodes))){
  locations <- list(gene_location, barcode_white_space, A_mutation_location) # all of the fragments that are not in the barcode
  for (l in A_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_barcodes[2], y=y_barcodes[y], r=plasmid_radius, w=plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(gene_thickness, barcode_border_thickness, rep(mutation_thickness,1), rep(barcode_thickness, length(A_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 1), rep('black', length(A_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 1), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_barcodes))){
  locations <- list(gene_location, barcode_white_space, B_mutation_location) # all of the fragments that are not in the barcode
  for (l in B_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_barcodes[3], y=y_barcodes[y], r=plasmid_radius, w=plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(gene_thickness, barcode_border_thickness, rep(mutation_thickness,1), rep(barcode_thickness, length(B_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 1), rep('black', length(B_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 1), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_barcodes))){
  locations <- list(gene_location, barcode_white_space, A_mutation_location, B_mutation_location, C_mutation_location) # all of the fragments that are not in the barcode
  for (l in ABC_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_barcodes[4], y=y_barcodes[y], r=plasmid_radius, w=plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(gene_thickness, barcode_border_thickness, rep(mutation_thickness,3), rep(barcode_thickness, length(ABC_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 3), rep('black', length(ABC_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 3), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
points(x = c(0.65,0.7,0.75), y = rep(y_barcodes[2], 3), pch = 16, lwd=0.5)
## cell rows -----
library_cell_width = 0.1
library_cell_length = 0.15
library_cell_spacing = 0.055
library_cell_thickness = 1
in_factor = 0.75
in_plasmid_radius <- plasmid_radius * in_factor * 0.9
in_plasmid_gene_width <- plasmid_gene_width * in_factor
in_gene_thickness <- gene_thickness * in_factor
in_mutation_thickness <- mutation_thickness * in_factor
in_barcode_thickness <- barcode_thickness * in_factor
y_red_center <- 0.35
y_blue_center <- 0.07
y_red_cells <- c(y_red_center+library_cell_spacing, y_red_center, y_red_center-library_cell_spacing)
x_red_cells <- c(0.05, 0.25, 0.45, 0.85)
y_blue_cells <- c(y_blue_center+library_cell_spacing, y_blue_center, y_blue_center-library_cell_spacing)
x_blue_cells <- c(0.15, 0.35, 0.55, 0.95)
#  plot red cell row -----
for (y in c(1:length(y_red_cells))) {
  draw.cell(center = c(x_red_cells[1], y_red_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(red_Ec_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space) # all of the fragments that are not in the barcode
  for (l in WT_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_red_cells[1], y=y_red_cells[y], r=in_plasmid_radius, w=in_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(in_gene_thickness, barcode_border_thickness, rep(in_barcode_thickness, length(WT_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep('black', length(WT_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_red_cells))) {
  draw.cell(center = c(x_red_cells[2], y_red_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(red_Ec_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, A_mutation_location) # all of the fragments that are not in the barcode
  for (l in A_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_red_cells[2], y=y_red_cells[y], r=in_plasmid_radius, w=in_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(in_gene_thickness, barcode_border_thickness, rep(in_mutation_thickness,1), rep(in_barcode_thickness, length(A_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 1), rep('black', length(A_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 1), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_red_cells))) {
  draw.cell(center = c(x_red_cells[3], y_red_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(red_Ec_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, B_mutation_location) # all of the fragments that are not in the barcode
  for (l in B_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_red_cells[3], y=y_red_cells[y], r=in_plasmid_radius, w=in_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(in_gene_thickness, barcode_border_thickness, rep(in_mutation_thickness,1), rep(in_barcode_thickness, length(B_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 1), rep('black', length(B_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 1), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_red_cells))) {
  draw.cell(center = c(x_red_cells[4], y_red_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(red_Ec_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, A_mutation_location, B_mutation_location, C_mutation_location) # all of the fragments that are not in the barcode
  for (l in ABC_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_red_cells[4], y=y_red_cells[y], r=in_plasmid_radius, w=in_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(in_gene_thickness, barcode_border_thickness, rep(in_mutation_thickness,3), rep(in_barcode_thickness, length(ABC_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 3), rep('black', length(ABC_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 3), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
points(x = c(0.6,0.65,0.7), y = rep(y_red_cells[2], 3), pch = 16, lwd=0.5, col = red_Ec_color)
#  plot blue cell row -----
for (y in c(1:length(y_blue_cells))){
  draw.cell(center = c(x_blue_cells[1], y_blue_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(blue_Kp_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space) # all of the fragments that are not in the barcode
  for (l in WT_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_blue_cells[1], y=y_blue_cells[y], r=in_plasmid_radius, w=in_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(in_gene_thickness, barcode_border_thickness, rep(in_barcode_thickness, length(WT_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep('black', barcode_background_color, length(WT_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_blue_cells))){
  draw.cell(center = c(x_blue_cells[2], y_blue_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(blue_Kp_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, A_mutation_location) # all of the fragments that are not in the barcode
  for (l in A_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_blue_cells[2], y=y_blue_cells[y], r=in_plasmid_radius, w=in_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(in_gene_thickness, barcode_border_thickness, rep(in_mutation_thickness,1), rep(in_barcode_thickness, length(A_barcodes_locations[[y]]))),
                    genes.col=c(gene_color,barcode_background_color,  rep(pie_color, 1), rep('black', length(A_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 1), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_blue_cells))){
  draw.cell(center = c(x_blue_cells[3], y_blue_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(blue_Kp_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, B_mutation_location) # all of the fragments that are not in the barcode
  for (l in B_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_blue_cells[3], y=y_blue_cells[y], r=in_plasmid_radius, w=in_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(in_gene_thickness, barcode_border_thickness, rep(in_mutation_thickness,1), rep(in_barcode_thickness, length(B_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 1), rep('black', length(B_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 1), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_blue_cells))){
  draw.cell(center = c(x_blue_cells[4], y_blue_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(blue_Kp_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, A_mutation_location, B_mutation_location, C_mutation_location) # all of the fragments that are not in the barcode
  for (l in ABC_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_blue_cells[4], y=y_blue_cells[y], r=in_plasmid_radius, w=in_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(in_gene_thickness, barcode_border_thickness, rep(in_mutation_thickness,3), rep(in_barcode_thickness, length(ABC_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 3), rep('black', length(ABC_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 3), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
points(x = c(0.7,0.75, 0.8), y = rep(y_blue_cells[2], 3), pch = 16, lwd=0.5, col = blue_Kp_color)
##### B panel ##### -----
## pooled libraries ----
cell_factor = 0.65
library_factor = 0.45
library_cell_width = cell_width * cell_factor
library_cell_length = cell_length * cell_factor
library_cell_spacing = cell_spacing * cell_factor
library_cell_offset = library_cell_spacing * 0.75
library_cell_thickness = cell_thickness * cell_factor
library_plasmid_radius <- plasmid_radius * library_factor * 0.9
library_plasmid_gene_width <- plasmid_gene_width * library_factor
library_gene_thickness <- gene_thickness * library_factor
library_mutation_thickness <- mutation_thickness * library_factor
library_barcode_thickness <- barcode_thickness * library_factor
blue_offset <- 0.05
y_red_center <- -0.275
y_blue_center <- -0.275
y_red_pooled_cells <- c(y_red_center+library_cell_spacing+library_cell_offset, y_red_center+library_cell_offset, y_red_center-library_cell_spacing+library_cell_offset)
y_red_pooled_cells_offset <- c(y_red_center+library_cell_spacing, y_red_center, y_red_center-library_cell_spacing)
x_red_pooled_cells <- c(0.01, 0.03, 0.05, 0.07)
y_blue_pooled_cells <- c(y_blue_center+library_cell_spacing+library_cell_offset, y_blue_center+library_cell_offset, y_blue_center-library_cell_spacing+library_cell_offset)
y_blue_pooled_cells_offset <- c(y_blue_center+library_cell_spacing, y_blue_center, y_blue_center-library_cell_spacing)
x_blue_pooled_cells <- c(0.51, 0.53, 0.55, 0.57) + blue_offset
#  plot red pooled libraries ----
for (y in c(1:length(y_red_pooled_cells))) {
  draw.cell(center = c(x_red_pooled_cells[1], y_red_pooled_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(red_Ec_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space) # all of the fragments that are not in the barcode
  for (l in WT_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_red_pooled_cells[1], y=y_red_pooled_cells[y], r=library_plasmid_radius, w=library_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(library_gene_thickness, barcode_border_thickness, rep(library_barcode_thickness, length(WT_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep('black', length(WT_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_red_pooled_cells_offset))) {
  draw.cell(center = c(x_red_pooled_cells[2], y_red_pooled_cells_offset[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(red_Ec_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, A_mutation_location) # all of the fragments that are not in the barcode
  for (l in A_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_red_pooled_cells[2], y=y_red_pooled_cells_offset[y], r=library_plasmid_radius, w=library_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(library_gene_thickness, barcode_border_thickness, rep(library_mutation_thickness,1), rep(library_barcode_thickness, length(A_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 1), rep('black', length(A_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 1), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_red_pooled_cells))) {
  draw.cell(center = c(x_red_pooled_cells[3], y_red_pooled_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(red_Ec_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, B_mutation_location) # all of the fragments that are not in the barcode
  for (l in B_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_red_pooled_cells[3], y=y_red_pooled_cells[y], r=library_plasmid_radius, w=library_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(library_gene_thickness, barcode_border_thickness, rep(library_mutation_thickness,1), rep(library_barcode_thickness, length(B_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 1), rep('black', length(B_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 1), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_red_pooled_cells_offset))) {
  draw.cell(center = c(x_red_pooled_cells[4], y_red_pooled_cells_offset[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(red_Ec_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, A_mutation_location, B_mutation_location, C_mutation_location) # all of the fragments that are not in the barcode
  for (l in ABC_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_red_pooled_cells[4], y=y_red_pooled_cells_offset[y], r=library_plasmid_radius, w=library_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(library_gene_thickness, barcode_border_thickness, rep(library_mutation_thickness,3), rep(library_barcode_thickness, length(ABC_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 3), rep('black', length(ABC_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 3), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
#  plot blue pooled libraries -----
for (y in c(1:length(y_blue_pooled_cells))){
  draw.cell(center = c(x_blue_pooled_cells[1], y_blue_pooled_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(blue_Kp_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space) # all of the fragments that are not in the barcode
  for (l in WT_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_blue_pooled_cells[1], y=y_blue_pooled_cells[y], r=library_plasmid_radius, w=library_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(library_gene_thickness, barcode_border_thickness, rep(library_barcode_thickness, length(WT_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep('black', barcode_background_color, length(WT_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_blue_pooled_cells_offset))){
  draw.cell(center = c(x_blue_pooled_cells[2], y_blue_pooled_cells_offset[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(blue_Kp_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, A_mutation_location) # all of the fragments that are not in the barcode
  for (l in A_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_blue_pooled_cells[2], y=y_blue_pooled_cells_offset[y], r=library_plasmid_radius, w=library_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(library_gene_thickness, barcode_border_thickness, rep(library_mutation_thickness,1), rep(library_barcode_thickness, length(A_barcodes_locations[[y]]))),
                    genes.col=c(gene_color,barcode_background_color,  rep(pie_color, 1), rep('black', length(A_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 1), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_blue_pooled_cells))){
  draw.cell(center = c(x_blue_pooled_cells[3], y_blue_pooled_cells[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(blue_Kp_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, B_mutation_location) # all of the fragments that are not in the barcode
  for (l in B_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_blue_pooled_cells[3], y=y_blue_pooled_cells[y], r=library_plasmid_radius, w=library_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(library_gene_thickness, barcode_border_thickness, rep(library_mutation_thickness,1), rep(library_barcode_thickness, length(B_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 1), rep('black', length(B_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 1), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
for (y in c(1:length(y_blue_pooled_cells_offset))){
  draw.cell(center = c(x_blue_pooled_cells[4], y_blue_pooled_cells_offset[y]), length = library_cell_length, width = library_cell_width, rounded = 1, angle = 0, cytoplasm = "white", membrane = c(blue_Kp_color,library_cell_thickness))
  locations <- list(gene_location, barcode_white_space, A_mutation_location, B_mutation_location, C_mutation_location) # all of the fragments that are not in the barcode
  for (l in ABC_barcodes_locations[[y]]){ # add all of the locations for the barcodes
    locations <- c(locations, list(c(l, l)))}
  draw.plasmid.back(x=x_blue_pooled_cells[4], y=y_blue_pooled_cells_offset[y], r=library_plasmid_radius, w=library_plasmid_gene_width,
                    genes.pos=locations,
                    genes.lwd=c(library_gene_thickness, barcode_border_thickness, rep(library_mutation_thickness,3), rep(library_barcode_thickness, length(ABC_barcodes_locations[[y]]))),
                    genes.col=c(gene_color, barcode_background_color, rep(pie_color, 3), rep('black', length(ABC_barcodes_locations[[y]]))), 
                    border.col = c('black', barcode_border_color, rep('black', 3), rep('black', length(WT_barcodes_locations[[y]]))), 
                    backbone.col="black", backbone.width=1, res=100)
}
## test tubes ----
tube_height <- 0.15
tube_width <- tube_height/(0.7/0.3)
tube_linewidth <- 1
x_tubes <- list(c(0.175,0.675+blue_offset), c(0.25,0.75+blue_offset), c(0.325,0.825+blue_offset), c(0.45,0.95+blue_offset))
y_tubes <- -0.3
tubes_color <- list(rgb(245, 245, 245, maxColorValue = 255),
                    rgb(235, 235, 235, maxColorValue = 255),
                    rgb(225, 225, 225, maxColorValue = 255),
                    rgb(150, 150, 150, maxColorValue = 255))
#  plot test tubes-----
for (t in 1:length(x_tubes)){
  for (h in 1:length(x_tubes[[1]])){
    draw.tube(center = c(x_tubes[[t]][h],y_tubes), width = tube_width, height = tube_height, liquid.level = .75, 
              side.view = .75, resolution=100, angle=0, liquid=tubes_color[[t]], 
              meniscus=c("black",tube_linewidth), glass=c("black",tube_linewidth))}}
x_points_red_tubes <- (x_tubes[[3]][1]+x_tubes[[4]][1])/2
x_points_blue_tubes <- (x_tubes[[3]][2]+x_tubes[[4]][2])/2
point_offset <- 0.015
points(x = c(x_points_red_tubes - point_offset,
             x_points_red_tubes,
             x_points_red_tubes + point_offset), y = rep(y_tubes, 3), pch = 16, cex=0.5)
points(x = c(x_points_blue_tubes - point_offset,
             x_points_blue_tubes,
             x_points_blue_tubes + point_offset), y = rep(y_tubes, 3), pch = 16, cex=0.5)
## amplicons ----
y_dish <- -0.65
amplicon_length <- 0.06
amplicon_spacing <- 0.04
amplicon_gene_width <- plasmid_gene_width * 0.5
y_amplicons <-c(y_dish + amplicon_spacing, y_dish, y_dish - amplicon_spacing)
x_red_amplicon <- 0.1
x_red_amplicon_white_space <- c(x_red_amplicon-(amplicon_length/3), 
                                x_red_amplicon+(amplicon_length/3))
x_blue_amplicon <- 0.6 + blue_offset
x_blue_amplicon_white_space <- c(x_blue_amplicon-(amplicon_length/3), 
                                x_blue_amplicon+(amplicon_length/3))
amplicon_factor <- 2
amplicon_thickness <- barcode_thickness*amplicon_factor
amplicon_border_thickness <- barcode_border_thickness*amplicon_factor
#  plot amplicons -----
for (a in 1:length(y_amplicons)) {
  locations <- list(x_red_amplicon_white_space)
  for (l in WT_barcodes[[a]]){
    x_value <- l * (x_red_amplicon_white_space[2]-x_red_amplicon_white_space[1]) + x_red_amplicon_white_space[1]
    locations <- c(locations, list(c(x_value, x_value)))
    }
  draw.segment(x=x_red_amplicon, y=y_amplicons[a], r=amplicon_length, w=amplicon_gene_width,
               genes.pos=locations, 
               genes.col=c('white', rep('black', length(WT_barcodes[[a]]))),
               border.col=c('white', rep('black', length(WT_barcodes[[a]]))), 
               genes.lwd = c(amplicon_border_thickness, rep(amplicon_thickness, length(WT_barcodes[[y]]))), backbone.col="black", backbone.width=1)}
for (a in 1:length(y_amplicons)) {
  locations <- list(x_blue_amplicon_white_space)
  for (l in WT_barcodes[[a]]){
    x_value <- l * (x_blue_amplicon_white_space[2]-x_blue_amplicon_white_space[1]) + x_blue_amplicon_white_space[1]
    locations <- c(locations, list(c(x_value, x_value)))
  }
  draw.segment(x=x_blue_amplicon, y=y_amplicons[a], r=amplicon_length, w=amplicon_gene_width,
               genes.pos=locations, 
               genes.col=c('white', rep('black', length(WT_barcodes[[a]]))),
               border.col=c('white', rep('black', length(WT_barcodes[[a]]))), 
               genes.lwd = c(amplicon_border_thickness, rep(amplicon_thickness, length(WT_barcodes[[y]]))), backbone.col="black", backbone.width=1)}

## petri dishes ----
center_red_dish <- c(0.375, y_dish)
center_blue_dish <- c(0.875+blue_offset, y_dish)
dish_side_view <- 0.55
dish_liquid_level <- 1
dish_width <- 0.225
dish_height <- 0.01
dish_resolution <- 50
dish_colony_radius <- 0.075
dish_colony_spacing <- 0.02
dish_color <- 'lightgoldenrod'
dish_border_thickness <- 1
dish_number_colonies <- 15
#  plot petri dishes ----
draw.dish.colonies(center = center_red_dish,
                   width = dish_width,
                   height = dish_height,
                   liquid.level = dish_liquid_level,
                   side.view = dish_side_view,
                   resolution=dish_resolution,
                   angle=0,
                   colony.radius=dish_colony_radius,
                   edge_spacing =dish_colony_radius*.5,
                   colony.color=red_Ec_color,
                   critical.spacing=dish_colony_spacing,
                   liquid=dish_color,
                   meniscus=c("black",dish_border_thickness),
                   glass=c("black",dish_border_thickness),
                   N.colonies = dish_number_colonies)
draw.dish.colonies(center = center_blue_dish,
                   width = dish_width,
                   height = dish_height,
                   liquid.level = dish_liquid_level,
                   side.view = dish_side_view,
                   resolution=dish_resolution,
                   angle=0,
                   colony.radius=dish_colony_radius,
                   edge_spacing =dish_colony_radius*.5,
                   colony.color=blue_Kp_color,
                   critical.spacing=dish_colony_spacing,
                   liquid=dish_color,
                   meniscus=c("black",dish_border_thickness),
                   glass=c("black",dish_border_thickness),
                   N.colonies = dish_number_colonies)
## arrows ----
bracket_thickness = cell_thickness
arrow_length <- 0.075
arrow_genotype_to_barcode_buffer = 0.01
#  plot arrows ---- 
for (a in 1:length(x_genotype_plasmid)){
  arrows(x0 = x_genotype_plasmid[a], x1 = x_genotype_plasmid[a],
         y0 = y_genotype_plasmid-plasmid_radius-arrow_genotype_to_barcode_buffer, 
         y1 = y_barcode_center+arrow_genotype_to_barcode_buffer+y_barcode_spacing*2, length = arrow_length, lwd = bracket_thickness, col = 'black')
}
## brackets---------
bracket_height = 0.05
bracket_curvature = 0.1
#  plot brackets -----
brackets(x1 = x_red_cells[4]+cell_length/2, 
         x2 = x_red_cells[1]-cell_length/2, 
         y1 = y_red_cells[3]-cell_spacing, 
         y2 = y_red_cells[3]-cell_spacing, 
         ticks = c(0.94), lwd=bracket_thickness, col = red_Ec_color,
         curvature = bracket_curvature,
         h = bracket_height)
brackets(x1 = x_blue_cells[4]+cell_length/2, 
         x2 = x_blue_cells[1]-cell_length/2, 
         y1 = y_blue_cells[3]-cell_spacing, 
         y2 = y_blue_cells[3]-cell_spacing, 
         ticks = c(0.5), lwd=bracket_thickness, col = blue_Kp_color,
         curvature = bracket_curvature,
         h = bracket_height)
y_tube_offset <- tube_height*0.77
brackets(x1 = x_tubes[[4]][1], 
         x2 = (x_red_pooled_cells[2] + x_red_pooled_cells[3])/2, 
         y1 = y_tubes - y_tube_offset, 
         y2 = y_tubes - y_tube_offset, 
         ticks = c(0.175,0.875), # right position, then left position 
         lwd=bracket_thickness, col = red_Ec_color,
         h = bracket_height)
arrows(x0 = (x_red_pooled_cells[2] + x_red_pooled_cells[3])/2, x1 = (x_red_pooled_cells[2] + x_red_pooled_cells[3])/2,
       y0 = y_tubes - y_tube_offset, 
       y1 = y_tubes - library_cell_spacing*2, length = 0, lwd = bracket_thickness, col = red_Ec_color)
arrows(x0 = x_tubes[[1]][1], x1 = x_tubes[[1]][1],
       y0 = y_tubes - y_tube_offset, 
       y1 = y_tubes - y_tube_offset - bracket_height/2, length = 0, lwd = bracket_thickness, col = red_Ec_color)
arrows(x0 = x_tubes[[2]][1], x1 = x_tubes[[2]][1],
       y0 = y_tubes - y_tube_offset, 
       y1 = y_tubes - y_tube_offset - bracket_height/2, length = 0, lwd = bracket_thickness, col = red_Ec_color)
arrows(x0 = x_tubes[[3]][1], x1 = x_tubes[[3]][1],
       y0 = y_tubes - y_tube_offset, 
       y1 = y_tubes - y_tube_offset - bracket_height/2, length = 0, lwd = bracket_thickness, col = red_Ec_color)

brackets(x1 = x_tubes[[4]][2], 
         x2 = (x_blue_pooled_cells[2] + x_blue_pooled_cells[3])/2, 
         y1 = y_tubes - y_tube_offset, 
         y2 = y_tubes - y_tube_offset, 
         ticks = c(0.175,0.875), lwd=bracket_thickness, col = blue_Kp_color,
         h = bracket_height)
arrows(x0 = (x_blue_pooled_cells[2] + x_blue_pooled_cells[3])/2, x1 = (x_blue_pooled_cells[2] + x_blue_pooled_cells[3])/2,
       y0 = y_tubes - y_tube_offset, 
       y1 = y_tubes - library_cell_spacing*2, length = 0, lwd = bracket_thickness, col = blue_Kp_color)
arrows(x0 = x_tubes[[1]][2], x1 = x_tubes[[1]][2],
       y0 = y_tubes - y_tube_offset, 
       y1 = y_tubes - y_tube_offset - bracket_height/2, length = 0, lwd = bracket_thickness, col = blue_Kp_color)
arrows(x0 = x_tubes[[2]][2], x1 = x_tubes[[2]][2],
       y0 = y_tubes - y_tube_offset, 
       y1 = y_tubes - y_tube_offset - bracket_height/2, length = 0, lwd = bracket_thickness, col = blue_Kp_color)
arrows(x0 = x_tubes[[3]][2], x1 = x_tubes[[3]][2],
       y0 = y_tubes - y_tube_offset, 
       y1 = y_tubes - y_tube_offset - bracket_height/2, length = 0, lwd = bracket_thickness, col = blue_Kp_color)

## text ----
y_text_barcode_seq_offset <- 0.05
x_text_barcode_seq_offset <- 0.05
x_text_plating_offset <- 0.35
y_white_coverup <- 0.01
y_text_offset <- 0.015
#  plot text -----
mtext("a", line = -1, side=3, at = -0.02, cex = global.cex)  
mtext("b", line = -0.5, side=2, at = -0.02, las =2, padj = 1, cex = global.cex)  
y_text_add_barcodes <- y_text_offset+((y_genotype_plasmid-plasmid_radius-arrow_genotype_to_barcode_buffer)+(y_barcode_center+arrow_genotype_to_barcode_buffer+y_barcode_spacing*2))/2
arrows(x0 = x_genotype_plasmid[1], x1 = x_genotype_plasmid[1],
       y0 = y_text_add_barcodes+y_white_coverup, 
       y1 = y_text_add_barcodes-y_white_coverup*1.75, length = 0, lwd = bracket_thickness*2, col = 'white')
text(labels = 'add barcodes', x = x_genotype_plasmid[1], y = y_text_add_barcodes)  
text(labels = 'transform', x = 0.1, y = 0.54)
text(labels = 'pool strains', x = 0.06, y = -.12)
text(labels = 'sequencing for\nfrequencies', 
     x = ((x_red_pooled_cells[2] + x_red_pooled_cells[3])/2) + x_text_barcode_seq_offset, 
     y = y_tubes - y_tube_offset - bracket_height - y_text_barcode_seq_offset)
text(labels = 'plating for\ndensity', 
     x = ((x_red_pooled_cells[2] + x_red_pooled_cells[3])/2) + x_text_plating_offset, 
     y = y_tubes - y_tube_offset - bracket_height - y_text_barcode_seq_offset)

##### C panel ##### -----
# left top -----
par(mai=c(T_bottom_in,L_left_in,T_top_in,L_right_in))
Dose_response_plot(df_reformat = df_reformat,
                   Data_fit = Data_fit,
                   host = 'Ec',
                   focal_color = red_Ec_color,
                   node.radius = 0.6, 
                   node.edge.width = 0.5,
                   axis.label = T, 
                   global.cex = global.cex,
                   xyratio = B_square_in/T_square_in)
mtext("c", side=2, line = 1, at = 11, cex = global.cex, las = 1) 
# right top -----
par(mai=c(T_bottom_in,R_left_in,T_top_in,R_right_in))
Dose_response_plot(df_reformat = df_reformat,
                   Data_fit = Data_fit,
                   host = 'Kp',
                   focal_color = blue_Kp_color,
                   node.radius = 0.6, 
                   node.edge.width = 0.5, 
                   global.cex = global.cex,
                   xyratio = B_square_in/T_square_in)
##### C-D bracket #### ----
plot.new()
##### D panel ##### -----
# left bottom ------
par(mai=c(B_bottom_in,L_left_in,B_top_in,L_right_in))
fitness.landscape.graph(df_filename = 'inputs/suboptimal.csv',
                        offsets_filename = "inputs/aligned_offsets.csv",
                        host = 'Ec',
                        host.color = red_Ec_color, 
                        node.edge.width = 0.5,
                        xyratio = 1,
                        axis.label = T,
                        radius.buffer = 1,
                        node.radius = 0.15,
                        global.cex = global.cex)
mtext("d", side=2, line = 1, at = 12, cex = global.cex, las = 1) 

# right bottom ------
par(mai=c(B_bottom_in,R_left_in,B_top_in,R_right_in))
fitness.landscape.graph(df_filename = 'inputs/suboptimal.csv',
                        offsets_filename = "inputs/aligned_offsets.csv",
                        host = 'Kp',
                        host.color = blue_Kp_color,
                        node.edge.width = 0.5,
                        xyratio = 1,
                        axis.label = F,
                        radius.buffer = 1,
                        node.radius = 0.15,
                        global.cex = global.cex)
dev.off()

