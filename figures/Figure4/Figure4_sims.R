# Packages-----
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
                            tick.cex,
                            y.main.label,
                            focal_color,
                            x_max,
                            y_max) {
  plot((0:x_max), seq(0, y_max, y_max/x_max), type="n", xlab='', ylab='', frame.plot=T, axes = F)  
  mtext(side=1, line=1.5, text =x.main.label, col='black', cex=global.cex)
  mtext(side=2, line=1.5, text =y.main.label, col='black',  cex=global.cex)
  axis(2, cex.axis=tick.cex, padj = 1.1, tck=-(0.05/square_h_in), lwd = 0, lwd.ticks = 1)
  axis(1, cex.axis=tick.cex, padj = -1.5, tck=-(0.05/square_h_in), lwd = 0, lwd.ticks = 1)
  polygon(x = c(df$Cumulative_time,rev(df$Cumulative_time)), y = c(df$Resistance.mean-df$Resistance.se, rev(df$Resistance.mean+df$Resistance.se)), col=adjustcolor(focal_color,alpha.f = 0.25), border=NA)
  lines(df$Cumulative_time, df$Resistance.mean, col = focal_color, lwd = 1.5)
}

simulation_plot_split <- function(df,
                                  df.mean,
                                  x.main.label,
                                  global.cex,
                                  tick.cex,
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
  axis(2, labels = FALSE, cex.axis=tick.cex, padj = -1.5, tck=-(0.05/square_h_in), lwd = 0, lwd.ticks = 1)
  axis(1, cex.axis=tick.cex, padj = -1.5, tck=-(0.05/square_h_in), lwd = 0, lwd.ticks = 1)
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
# Box plot functions ##############
MinMeanSEMMax <- function(x) {
  v <- c(max(x), 
         mean(x) - sd(x)/sqrt(length(x)), 
         mean(x), 
         mean(x) + sd(x)/sqrt(length(x)), 
         min(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}
boxplot_func <- function(y, positions, color) {
  bp_vals <- MinMeanSEMMax(y)
  box_width <- 0.2
  # Draw the bottom of the box
  segments(
    positions - box_width, bp_vals["lower"], positions + box_width, bp_vals["lower"],
    col = color,
    lwd = 1.5
  )
  # Draw the top of the box
  segments(
    positions - box_width, bp_vals["upper"], positions + box_width, bp_vals["upper"],
    col = color,
    lwd = 1.5
  )
  # Draw the left line of the box
  segments(
    positions - box_width, bp_vals["lower"], positions - box_width, bp_vals["upper"],
    col = color,
    lwd = 1.5
  )
  # Draw the right line of the box
  segments(
    positions + box_width, bp_vals["lower"], positions + box_width, bp_vals["upper"],
    col = color,
    lwd = 1.5
  )
  # Draw bottom whisker
  segments(
    positions - (box_width/2), bp_vals["ymin"], positions + (box_width/2), bp_vals["ymin"],
    col = color,
    lwd = 1.5
  )
  # Draw top whisker
  segments(
    positions - (box_width/2), bp_vals["ymax"], positions + (box_width/2), bp_vals["ymax"],
    col = color,
    lwd = 1.5
  )
  # Connect top whisker to top of the box
  segments(
    positions, bp_vals["upper"], positions, bp_vals["ymax"],
    col = color,
    lwd = 1.5
  )
  # Connect bottom whisker to bottom of the box
  segments(
    positions, bp_vals["lower"], positions, bp_vals["ymin"],
    col = color,
    lwd = 1.5
  )
  # Draw the middle of the box
  segments(
    positions - box_width, bp_vals["middle"], positions + box_width, bp_vals["middle"],
    col = color,
    lwd = 3
  )
}
box_plot_function <- function(df,
                              focal_host,
                              mean_value,
                              color_focal_set,
                              color_HGT_first_set,
                              color_HGT_second_set,
                              Focal.only.label,
                              HGT.first.label,
                              HGT.second.label,
                              tick.cex,
                              tick.cex.plot) {
  df_focal <- df[df$focal == focal_host, ]
  y_axis_addition <- max(df_focal$Resistance) * 0.15
  y_max <- max(df_focal$Resistance + y_axis_addition)
  plot(1, type="n", xlim=c(0.5, 3.5), ylim=c(0, y_max), xaxt="n", yaxt="n", xlab="", ylab="", axes=FALSE)
  box(bty = "L", lwd = 1)
  axis(2, labels = FALSE, padj = -1.5, tck=-(0.05/square_h_in), lwd = 0, lwd.ticks = 1)
  
  # Add x-axis labels
  ## focal only label
  mtext(side=1, line=0, at = 1, text = c(Focal.only.label[1], Focal.only.label[2]), 
        col=c("black", color_focal_set[1]), cex=tick.cex)
  
  # HGT first label
  mtext(side=1, line=0, at = 2, text = c(HGT.first.label[1], HGT.first.label[2]), 
        col=c("black", color_focal_set[1]), cex=tick.cex)
  mtext(side=1, line=0.75, at = 2, text = c(HGT.first.label[3], HGT.first.label[4]), 
        col=c("black", color_HGT_first_set[1]), cex=tick.cex)
  
  # HGT second label
  mtext(side=1, line=0, at = 3, text = c(HGT.second.label[1], HGT.second.label[2]), 
        col=c("black", color_focal_set[1]), cex=tick.cex)
  mtext(side=1, line=0.75, at = 3, text = c(HGT.second.label[3], HGT.second.label[4]), 
        col=c("black", color_HGT_second_set[1]), cex=tick.cex)
  
  # Loop over treatments and draw boxplots
  treatments <- unique(df_focal$Treatment_ID)
  positions <- 1:length(treatments)
  for (i in 1:length(treatments)) {
    treatment <- treatments[i]
    y_values <- df_focal[df_focal$Treatment_ID == treatment, "Resistance"]
    boxplot_func(y_values, positions[i], color_focal_set[i])
  }
  
  # add dashed line for the mean (replace with your desired y-value)
  abline(h=mean_value, lty=2, col=color_focal_set[1], lwd = 1.5)
  
  # add data points 
  alpha_decimal <- 0.05  # Transparency level (0.00 to 1.00)
  alpha_hex <- sprintf("%02X", round(alpha_decimal * 255))
  for (i in 1:length(treatments)) {
    treatment <- treatments[i]
    y_values <- df_focal[df_focal$Treatment_ID == treatment, "Resistance"]
    jittered_x <- rep(i, length(y_values)) + runif(length(y_values), -0.2, 0.2)
    points(jittered_x, y_values, col=paste0(color_focal_set[i], alpha_hex), pch=16)}
  
  # add brackets
  ## left brackets
  segments(1, y_max-(y_axis_addition*0.75), 1, y_max-(y_axis_addition*0.6), col = 'black', lwd = 1)
  segments(1, y_max-(y_axis_addition*0.25), 1, y_max-(y_axis_addition*0.1), col = 'black', lwd = 1)
  ## right brackets
  segments(2, y_max-(y_axis_addition*0.25), 2, y_max-(y_axis_addition*0.1), col = 'black', lwd = 1)
  segments(3, y_max-(y_axis_addition*0.75), 3, y_max-(y_axis_addition*0.6), col = 'black', lwd = 1)
  ## flat lines
  segments(1, y_max-(y_axis_addition*0.6), 3, y_max-(y_axis_addition*0.6), col = 'black', lwd = 1)
  segments(1, y_max-(y_axis_addition*0.1), 2, y_max-(y_axis_addition*0.1), col = 'black', lwd = 1)
  ## labels
  text(1.5, y_max - (y_axis_addition*0.25), "ns", pos = 3, col = 'black', cex=tick.cex.plot)
  text(2, y_max - (y_axis_addition*0.7), "ns", pos = 3, col = 'black', cex=tick.cex.plot)
}
# Arrange plot --------------
### Figure  --------------
matrix_layout <- as.matrix(read.csv('inputs/matrix.csv', header = F))
n.row = nrow(matrix_layout)
n.col = ncol(matrix_layout)
outer_margin = 0.01
tick_marks = 0.05
plot_spacer = 0.1
pointsize = 12
globalpointsize = 9
tickpointsize = 7
global.cex = globalpointsize / pointsize
tick.cex = tickpointsize / pointsize
tick.cex.sim = tickpointsize / globalpointsize
# columns
## generic 
c_right_margins = plot_spacer 
c_left_margins = 0 + tick_marks
## specific
c1_left_margin = 0.3 + tick_marks
c4_right_margin = 0.1
c4_left_margin = 0.1
# rows
## generic 
r_bottom_margin = 0.4 + tick_marks
r_top_margin = 0
## specific
r1_top_margin = 0.2
# plot square
square_w_in = 1.3
square_h_in = 1.3
c1 = c1_left_margin + c_right_margins + square_w_in
c2 = c_left_margins + c_right_margins + square_w_in  
c3 = c_left_margins + c_right_margins + square_w_in 
c4 = c4_left_margin + c4_right_margin + square_w_in
r1 = r_bottom_margin + r1_top_margin + square_h_in
r = r_bottom_margin + r_top_margin + square_h_in
pdf.width = c1 + c2 + c3 + c4 + (outer_margin*2) # 6.5 max
pdf.height = r1 + r + r +(outer_margin*2) # 9.25 max
### Initiate  --------
pdf(file = 'Figure4.pdf', width = pdf.width, height = pdf.height, pointsize = pointsize)
layout(matrix(matrix_layout, nrow = n.row, ncol = n.col),
       widths = c(lcm(c1*2.54),lcm(c2*2.54),lcm(c3*2.54),lcm(c4*2.54)), 
       heights = c(lcm(r1*2.54),lcm(r*2.54),lcm(r*2.54)))
par(omi=c(outer_margin,outer_margin,outer_margin,outer_margin)) # outer margins (b, l, t, r)
subplot_vertical_location = 0.2
subplot_simfocal_labl = -14
subplot_barplot_label = 0.15
### Panel A, B, C---------
par(mai=c(r_bottom_margin,c1_left_margin,r1_top_margin,c_right_margins), cex = global.cex)
df_Ec_only <- simulation_data_load(filename = "inputs/76.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Ec_Kp <- simulation_data_load(filename = "inputs/77.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Ec_Se <- simulation_data_load(filename = "inputs/78.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
y_max <- max(read.csv("inputs/76.csv")$Resistance)
y_axis_addition <- y_max * 0.15
sim_ymax <- y_max + y_axis_addition
sim_ymax = 875
simulation_plot(df = df_Ec_only, 
                x.main.label = 'time', 
                global.cex = global.cex,
                tick.cex = tick.cex.sim,
                y.main.label = "resistance", 
                focal_color = red_Ec_color, 
                x_max = 60, 
                y_max = sim_ymax)
mtext("a", line = subplot_vertical_location, side=3, at = subplot_simfocal_labl)  
par(mai=c(r_bottom_margin,c_left_margins,r1_top_margin,c_right_margins), cex = global.cex)
arrow_prop = 17.5
line_prop_mult = 4
buffer_prop = c(1000, 40) # was 80
simulation_plot_split(df = df_Ec_Kp,
                      df.mean = df_Ec_only, 
                      x.main.label = 'time', 
                      global.cex = global.cex,
                      tick.cex = tick.cex.sim, 
                      y.main.label = "", 
                      combine_color = purple_Ec.Kp_color, 
                      focal_color = red_Ec_color, 
                      transient_color = blue_Kp_color, 
                      x_max = 60, 
                      y_max = sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
#mtext("b", line = -0.5, side=3, at = -5)  
par(mai=c(r_bottom_margin,c_left_margins,r1_top_margin,c_right_margins), cex = global.cex)
buffer_prop = c(80, 20)
simulation_plot_split(df = df_Ec_Se, 
                      df.mean = df_Ec_only,
                      x.main.label = 'time', 
                      global.cex = global.cex,
                      tick.cex = tick.cex.sim, 
                      y.main.label = "", 
                      combine_color = orange_Ec.Se_color, 
                      focal_color = red_Ec_color, 
                      transient_color = yellow_Se_color, 
                      x_max = 60, 
                      y_max = sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
#mtext("c", line = -0.5, side=3, at = -5)  
### Panel D ########
par(mai=c(r_bottom_margin,c4_left_margin,r1_top_margin,c4_right_margin), cex = global.cex)
# Load data (replace with your data loading code)
df_Ec_only_box <- read.csv("inputs/76.csv")
df_Ec_Kp_box <- read.csv("inputs/77.csv")
df_Ec_Se_box <- read.csv("inputs/78.csv")
df_Ec_only_box$focal <- 'Ec'
df_Ec_Kp_box$focal <- 'Ec'
df_Ec_Se_box$focal <- 'Ec'
df_Ec_rbind <- rbind(df_Ec_only_box, df_Ec_Kp_box, df_Ec_Se_box)
df_Ec <- df_Ec_rbind %>%
  filter(Cumulative_time == 60)

Ec_mean <- read.csv("inputs/76_averages.csv")
Ec.mean <- Ec_mean %>% filter(Cumulative_time == 60) %>% select(Resistance.mean)

color_Ec = c(red_Ec_color, purple_Ec.Kp_color, orange_Ec.Se_color)
color_Kp = c(blue_Kp_color, purple_Ec.Kp_color, green_Se.Kp_color)
color_Se = c(yellow_Se_color, orange_Ec.Se_color, green_Se.Kp_color)

Ec.only.main.label = expression(phantom('Ec') * ' only')
Ec.only.sp.label = expression('Ec' * phantom(' only'))
Ec.only.label = c(Ec.only.main.label, Ec.only.sp.label)

Ec.Kp.main.label.line1 = expression(phantom('Ec') * ' with')
Ec.Kp.sp.Ec.label.line1 = expression('Ec' * phantom(' with'))
Ec.Kp.label.line1 = c(Ec.Kp.main.label.line1, Ec.Kp.sp.Ec.label.line1)
Ec.Kp.main.label.line2 = expression('HGT to ' * phantom('Kp'))
Ec.Kp.sp.Kp.label.line2 = expression(phantom('HGT to ') * 'Kp')
Ec.Kp.label.line2 = c(Ec.Kp.main.label.line2, Ec.Kp.sp.Kp.label.line2)
Ec.Kp.label <- c(Ec.Kp.label.line1, Ec.Kp.label.line2)

Ec.Se.main.label.line1 = expression(phantom('Ec') * ' with HGT')
Ec.Se.sp.Ec.label.line1 = expression('Ec' * phantom(' with HGT'))
Ec.Se.label.line1 = c(Ec.Se.main.label.line1, Ec.Se.sp.Ec.label.line1)
Ec.Se.main.label.line2 = expression('to ' * phantom('Se'))
Ec.Se.sp.Se.label.line2 = expression(phantom('to ') * 'Se')
Ec.Se.label.line2 = c(Ec.Se.main.label.line2, Ec.Se.sp.Se.label.line2)
Ec.Se.label <- c(Ec.Se.label.line1, Ec.Se.label.line2)

box_plot_function(df = df_Ec, 
                  focal_host = 'Ec',
                  mean_value = Ec.mean,
                  color_focal_set = color_Ec,
                  color_HGT_first_set = color_Kp,
                  color_HGT_second_set = color_Se, 
                  Focal.only.label = Ec.only.label,
                  HGT.first.label = Ec.Kp.label,
                  HGT.second.label = Ec.Se.label,
                  tick.cex = tick.cex, 
                  tick.cex.plot = tick.cex.sim)
mtext("b", line = subplot_vertical_location, side=3, at = subplot_barplot_label)  

### Panel E, F, G---------
par(mai=c(r_bottom_margin,c1_left_margin,r_top_margin,c_right_margins), cex = global.cex)
df_Kp_only <- simulation_data_load(filename = "inputs/79.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Kp_Ec <- simulation_data_load(filename = "inputs/80.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Kp_Se <- simulation_data_load(filename = "inputs/81.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
Kp_sim_ymax = 800
simulation_plot(df = df_Kp_only, 
                x.main.label = 'time', 
                global.cex = global.cex,
                tick.cex = tick.cex.sim, 
                y.main.label = "resistance", 
                focal_color = blue_Kp_color, 
                x_max = 60, 
                y_max = Kp_sim_ymax)
mtext("c", line = subplot_vertical_location, side=3, at = subplot_simfocal_labl)  
par(mai=c(r_bottom_margin,c_left_margins,r_top_margin,c_right_margins), cex = global.cex)
arrow_prop = 17.5
line_prop_mult = 4
buffer_prop = c(100, 40) # was 80
simulation_plot_split(df = df_Kp_Ec, 
                      df.mean = df_Kp_only,
                      x.main.label = 'time', 
                      global.cex = global.cex,
                      tick.cex = tick.cex.sim, 
                      y.main.label = "", 
                      combine_color = purple_Ec.Kp_color, 
                      focal_color = blue_Kp_color, 
                      transient_color = red_Ec_color, 
                      x_max = 60, 
                      y_max = Kp_sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
#mtext("f", line = -0.5, side=3, at = -5)  
par(mai=c(r_bottom_margin,c_left_margins,r_top_margin,c_right_margins), cex = global.cex)
buffer_prop = c(80, 20)
simulation_plot_split(df = df_Kp_Se, 
                      df.mean = df_Kp_only,
                      x.main.label = 'time', 
                      global.cex = global.cex,
                      tick.cex = tick.cex.sim, 
                      y.main.label = "", 
                      combine_color = green_Se.Kp_color, 
                      focal_color = blue_Kp_color, 
                      transient_color = yellow_Se_color, 
                      x_max = 60, 
                      y_max = Kp_sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
#mtext("g", line = -0.5, side=3, at = -5)  
### Panel H ########
par(mai=c(r_bottom_margin,c4_left_margin,r_top_margin,c4_right_margin), cex = global.cex)
df_Kp_only_box <- read.csv("inputs/79.csv")
df_Kp_Ec_box <- read.csv("inputs/80.csv")
df_Kp_Se_box <- read.csv("inputs/81.csv")
df_Kp_only_box$focal <- 'Kp'
df_Kp_Ec_box$focal <- 'Kp'
df_Kp_Se_box$focal <- 'Kp'
df_Kp_rbind <- rbind(df_Kp_only_box, df_Kp_Ec_box, df_Kp_Se_box)
df_Kp <- df_Kp_rbind %>%
  filter(Cumulative_time == 60)

Kp_mean <- read.csv("inputs/79_averages.csv")
Kp.mean <- Kp_mean %>% filter(Cumulative_time == 60) %>% select(Resistance.mean)

Kp.only.main.label = expression(phantom('Kp') * ' only')
Kp.only.sp.label = expression('Kp' * phantom(' only'))
Kp.only.label = c(Kp.only.main.label, Kp.only.sp.label)

Kp.Ec.main.label.line1 = expression(phantom('Kp') * ' with')
Kp.Ec.sp.Kp.label.line1 = expression('Kp' * phantom(' with'))
Kp.Ec.label.line1 = c(Kp.Ec.main.label.line1, Kp.Ec.sp.Kp.label.line1)
Kp.Ec.main.label.line2 = expression('HGT to ' * phantom('Ec'))
Kp.Ec.sp.Ec.label.line2 = expression(phantom('HGT to ') * 'Ec')
Kp.Ec.label.line2 = c(Kp.Ec.main.label.line2, Kp.Ec.sp.Ec.label.line2)
Kp.Ec.label <- c(Kp.Ec.label.line1, Kp.Ec.label.line2)

Kp.Se.main.label.line1 = expression(phantom('Kp') * ' with HGT')
Kp.Se.sp.Kp.label.line1 = expression('Kp' * phantom(' with HGT'))
Kp.Se.label.line1 = c(Kp.Se.main.label.line1, Kp.Se.sp.Kp.label.line1)
Kp.Se.main.label.line2 = expression('to ' * phantom('Se'))
Kp.Se.sp.Se.label.line2 = expression(phantom('to ') * 'Se')
Kp.Se.label.line2 = c(Kp.Se.main.label.line2, Kp.Se.sp.Se.label.line2)
Kp.Se.label <- c(Kp.Se.label.line1, Kp.Se.label.line2)

box_plot_function(df = df_Kp, 
                  focal_host = 'Kp',
                  mean_value = Kp.mean,
                  color_focal_set = color_Kp,
                  color_HGT_first_set = color_Ec,
                  color_HGT_second_set = color_Se, 
                  Focal.only.label = Kp.only.label,
                  HGT.first.label = Kp.Ec.label,
                  HGT.second.label = Kp.Se.label,
                  tick.cex = tick.cex, 
                  tick.cex.plot = tick.cex.sim)
mtext("d", line = subplot_vertical_location, side=3, at = subplot_barplot_label)  

### Panel I, J, K---------
par(mai=c(r_bottom_margin,c1_left_margin,r_top_margin,c_right_margins), cex = global.cex)
df_Se_only <- simulation_data_load(filename = "inputs/82.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Se_Ec <- simulation_data_load(filename = "inputs/83.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
df_Se_Kp <- simulation_data_load(filename = "inputs/84.csv", reference_filename = "inputs/MIC_simulation_file_3barcodes.csv")
Se_sim_ymax = 1400
simulation_plot(df = df_Se_only, 
                x.main.label = 'time', 
                global.cex = global.cex,
                tick.cex = tick.cex.sim, 
                y.main.label = "resistance", 
                focal_color = yellow_Se_color, 
                x_max = 60, 
                y_max = Se_sim_ymax)
mtext("e", line = subplot_vertical_location, side=3, at = subplot_simfocal_labl)  
par(mai=c(r_bottom_margin,c_left_margins,r_top_margin,c_right_margins), cex = global.cex)
arrow_prop = 17.5
line_prop_mult = 4
buffer_prop = c(1000, 40) # was 80
simulation_plot_split(df = df_Se_Ec, 
                      df.mean = df_Se_only, 
                      x.main.label = 'time', 
                      global.cex = global.cex,
                      tick.cex = tick.cex.sim, 
                      y.main.label = "", 
                      combine_color = orange_Ec.Se_color, 
                      focal_color = yellow_Se_color, 
                      transient_color = red_Ec_color, 
                      x_max = 60, 
                      y_max = Se_sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
# mtext("j", line = -0.5, side=3, at = -5)  
par(mai=c(r_bottom_margin,c_left_margins,r_top_margin,c_right_margins), cex = global.cex)
buffer_prop = c(200, 40)
simulation_plot_split(df = df_Se_Kp, 
                      df.mean = df_Se_only,
                      x.main.label = 'time', 
                      global.cex = global.cex,
                      tick.cex = tick.cex.sim, 
                      y.main.label = "", 
                      combine_color = green_Se.Kp_color, 
                      focal_color = yellow_Se_color, 
                      transient_color = blue_Kp_color, 
                      x_max = 60, 
                      y_max = Se_sim_ymax,
                      arrow_proportion = arrow_prop,
                      line_proportion_multiplier = line_prop_mult,
                      buffer_proportion = buffer_prop)
# mtext("k", line = -0.5, side=3, at = -5)  
### Panel L ########
par(mai=c(r_bottom_margin,c4_left_margin,r_top_margin,c4_right_margin), cex = global.cex)
df_Se_only_box <- read.csv("inputs/82.csv")
df_Se_Ec_box <- read.csv("inputs/83.csv")
df_Se_Kp_box <- read.csv("inputs/84.csv")
df_Se_only_box$focal <- 'Se'
df_Se_Ec_box$focal <- 'Se'
df_Se_Kp_box$focal <- 'Se'
df_Se_rbind <- rbind(df_Se_only_box, df_Se_Ec_box, df_Se_Kp_box)
df_Se <- df_Se_rbind %>%
  filter(Cumulative_time == 60)

Se_mean <- read.csv("inputs/82_averages.csv")
Se.mean <- Se_mean %>% filter(Cumulative_time == 60) %>% select(Resistance.mean)

Se.only.main.label = expression(phantom('Se') * ' only')
Se.only.sp.label = expression('Se' * phantom(' only'))
Se.only.label = c(Se.only.main.label, Se.only.sp.label)

Se.Ec.main.label.line1 = expression(phantom('Se') * ' with')
Se.Ec.sp.Kp.label.line1 = expression('Se' * phantom(' with'))
Se.Ec.label.line1 = c(Se.Ec.main.label.line1, Se.Ec.sp.Kp.label.line1)
Se.Ec.main.label.line2 = expression('HGT to ' * phantom('Ec'))
Se.Ec.sp.Ec.label.line2 = expression(phantom('HGT to ') * 'Ec')
Se.Ec.label.line2 = c(Se.Ec.main.label.line2, Se.Ec.sp.Ec.label.line2)
Se.Ec.label <- c(Se.Ec.label.line1, Se.Ec.label.line2)

Se.Kp.main.label.line1 = expression(phantom('Se') * ' with HGT')
Se.Kp.sp.Se.label.line1 = expression('Se' * phantom(' with HGT'))
Se.Kp.label.line1 = c(Se.Kp.main.label.line1, Se.Kp.sp.Se.label.line1)
Se.Kp.main.label.line2 = expression('to ' * phantom('Kp'))
Se.Kp.sp.Kp.label.line2 = expression(phantom('to ') * 'Kp')
Se.Kp.label.line2 = c(Se.Kp.main.label.line2, Se.Kp.sp.Kp.label.line2)
Se.Kp.label <- c(Se.Kp.label.line1, Se.Kp.label.line2)

box_plot_function(df = df_Se, 
                  focal_host = 'Se',
                  mean_value = Se.mean,
                  color_focal_set = color_Se,
                  color_HGT_first_set = color_Ec,
                  color_HGT_second_set = color_Kp, 
                  Focal.only.label = Se.only.label,
                  HGT.first.label = Se.Ec.label,
                  HGT.second.label = Se.Kp.label,
                  tick.cex = tick.cex, 
                  tick.cex.plot = tick.cex.sim)
mtext("f", line = subplot_vertical_location, side=3, at = subplot_barplot_label)  
dev.off()






# trash ###########
# Set up the plotting area
par(mfrow=c(1, 1))  # 1 row and 1 column
dev.off()
plot.new()

#function arguments
focal_host <- 'Ec'
mean_value <- Ec.mean
color_focal_set <- color_Ec
color_HGT_first_set <- color_Kp
color_HGT_second_set <- color_Se
Focal.only.label <- Ec.only.label
HGT.first.label <- Ec.Kp.label
HGT.second.label <- Ec.Se.label
box_plot_function(df = df, 
                  focal_host = 'Ec',
                  mean_value = Ec.mean,
                  color_focal_set = color_Ec,
                  color_HGT_first_set = color_Kp,
                  color_HGT_second_set = color_Se, 
                  Focal.only.label = Ec.only.label,
                  HGT.first.label = Ec.Kp.label,
                  HGT.second.label = Ec.Se.label)




library(ggplot2)
library(cowplot)

# Create a ggplot2 plot
gg_plot <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()

# Create a base R plot using the plot function
plot(mtcars$wt, mtcars$mpg, xlab = "Weight", ylab = "Miles Per Gallon")
p1 <- recordPlot() 

# Combine plots using cowplot
combined_plot <- plot_grid(p1, E, K, S, labels = c("A", "B"), ncol = 2)

# Display the combined plot
print(combined_plot)

p1 <- recordPlot() 
