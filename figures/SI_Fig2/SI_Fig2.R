# Packages-----
## Download and install the package
#install.packages("igraph")
## Load package
library(igraph)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(gridGraphics)
library(readxl)
library(scales)
library(ggpubr)
library(rstatix)  
library(ggpmisc)
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
### SI Figure  --------------
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
r5_top_margin = 0.2
r5_bottom_margin = 0.4 + tick_marks
r6_top_margin = 0
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
r5 = r5_bottom_margin + r5_top_margin + d_square_h_in  # 1.95
r6 = r5_bottom_margin + r6_top_margin + d_square_h_in
pdf.width = c1 + c2 + c3 + (outer_margin*2) # 6.5 max
pdf.height = r5 + r6 + (outer_margin*2) # 9.25 max
### Initiate  --------
pointsize = 12
globalpointsize = 9
tickpointsize = 7
pdf(file = 'SI_Fig2.pdf', width = pdf.width, height = pdf.height, pointsize = pointsize)
layout(matrix(matrix_layout, nrow = n.row, ncol = n.col),
       widths = c(lcm(c1*2.54),lcm(c2*2.54),lcm(c3*2.54)), 
       heights = c(lcm(r5*2.54),lcm(r6*2.54)))
par(omi=c(outer_margin,outer_margin,outer_margin,outer_margin)) # outer margins (b, l, t, r)
global.cex = globalpointsize / pointsize
tick.cex = tickpointsize / globalpointsize

### Panel A, B, C---------
Kp_sim_ymax = 800
par(mai=c(r5_bottom_margin,c1_left_margin,r5_top_margin,c1_right_margin), cex = global.cex)
df_Kp_only <- simulation_data_load(filename = "inputs/79.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Kp_Ec <- simulation_data_load(filename = "inputs/80.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Kp_Se <- simulation_data_load(filename = "inputs/81.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
simulation_plot(df = df_Kp_only, 
                x.main.label = 'time', 
                global.cex = global.cex, 
                y.main.label = "resistance", 
                focal_color = blue_Kp_color, 
                x_max = 60, 
                y_max = Kp_sim_ymax)
mtext("a", line = -0.5, side=3, at = -18.5)  
par(mai=c(r5_bottom_margin,c2_left_margin,r5_top_margin,c2_right_margin), cex = global.cex)
arrow_prop = 17.5
line_prop_mult = 4
buffer_prop = c(100, 40) # was 80
simulation_plot_split(df = df_Kp_Ec, 
                      df.mean = df_Kp_only,
                      x.main.label = 'time', 
                      global.cex = global.cex, 
                      y.main.label = "", 
                      combine_color = purple_Ec.Kp_color, 
                      focal_color = blue_Kp_color, 
                      transient_color = red_Ec_color, 
                      x_max = 60, 
                      y_max = Kp_sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
mtext("b", line = -0.5, side=3, at = -8)  
par(mai=c(r5_bottom_margin,c3_left_margin,r5_top_margin,c3_right_margin), cex = global.cex)
buffer_prop = c(80, 20)
simulation_plot_split(df = df_Kp_Se, 
                      df.mean = df_Kp_only,
                      x.main.label = 'time', 
                      global.cex = global.cex, 
                      y.main.label = "", 
                      combine_color = green_Se.Kp_color, 
                      focal_color = blue_Kp_color, 
                      transient_color = yellow_Se_color, 
                      x_max = 60, 
                      y_max = Kp_sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
mtext("c", line = -0.5, side=3, at = -8)  
### Panel D, E, F---------
Se_sim_ymax = 1400
par(mai=c(r5_bottom_margin,c1_left_margin,r6_top_margin,c1_right_margin), cex = global.cex)
df_Se_only <- simulation_data_load(filename = "inputs/82.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Se_Ec <- simulation_data_load(filename = "inputs/83.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Se_Kp <- simulation_data_load(filename = "inputs/84.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
simulation_plot(df = df_Se_only, 
                x.main.label = 'time', 
                global.cex = global.cex, 
                y.main.label = "resistance", 
                focal_color = yellow_Se_color, 
                x_max = 60, 
                y_max = Se_sim_ymax)
mtext("d", line = -0.5, side=3, at = -18.5)  
par(mai=c(r5_bottom_margin,c2_left_margin,r6_top_margin,c2_right_margin), cex = global.cex)
arrow_prop = 17.5
line_prop_mult = 4
buffer_prop = c(1000, 40) # was 80
simulation_plot_split(df = df_Se_Ec, 
                      df.mean = df_Se_only, 
                      x.main.label = 'time', 
                      global.cex = global.cex, 
                      y.main.label = "", 
                      combine_color = orange_Ec.Se_color, 
                      focal_color = yellow_Se_color, 
                      transient_color = red_Ec_color, 
                      x_max = 60, 
                      y_max = Se_sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
mtext("e", line = -0.5, side=3, at = -8)  
par(mai=c(r5_bottom_margin,c3_left_margin,r6_top_margin,c3_right_margin), cex = global.cex)
buffer_prop = c(200, 40)
simulation_plot_split(df = df_Se_Kp, 
                      df.mean = df_Se_only,
                      x.main.label = 'time', 
                      global.cex = global.cex, 
                      y.main.label = "", 
                      combine_color = green_Se.Kp_color, 
                      focal_color = yellow_Se_color, 
                      transient_color = blue_Kp_color, 
                      x_max = 60, 
                      y_max = Se_sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
mtext("f", line = -0.5, side=3, at = -8)  
dev.off()




