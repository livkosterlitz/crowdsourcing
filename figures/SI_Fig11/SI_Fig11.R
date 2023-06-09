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
library(patchwork)
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
color_Ec = c(red_Ec_color, purple_Ec.Kp_color, orange_Ec.Se_color)
color_Kp = c(blue_Kp_color, purple_Ec.Kp_color, green_Se.Kp_color)
color_Se = c(yellow_Se_color, orange_Ec.Se_color, green_Se.Kp_color)

# General plotting customization -------------------------------------------------------
Figure_width = 3.25 # in width for one column figure in PNAS
Figure_height = 4
Figure_label_size = 12
text_size_axis_title = 9
text_size_axis_tick_labels = 7
axis_line_size = 0.75/(ggplot2::.pt*72.27/96) # pt size for the axis lines
axis_tick_size = 0.75/(ggplot2::.pt*72.27/96)
axis_tick_lengths = 0.05
y_axis_min = 0
y_axis_max = 900
y_axis_number_of_ticks = 8
y_axis_label = expression(paste("resistance"))
point_jitter_width = 0.15
point_size = 1
margins<- c(0.1,0.1,0.1,0.1)
box_line_size = 0.75/(ggplot2::.pt*72.27/96) 

# Functions -------------------------------------------------------
MinMeanSEMMax <- function(x) {
  v <- c(min(x), 
         mean(x) - sd(x)/sqrt(length(x)), 
         mean(x), 
         mean(x) + sd(x)/sqrt(length(x)), 
         max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

Parameter_sweep_load_data <- function(treatments_filename, treatment_information, sweep_variable) {
  treatment_information <- read.csv(treatment_information)
  csvfiles <- list.files(treatments_filename, full.names = TRUE)
  csvfileslist <- lapply(csvfiles, read.csv, header = TRUE, stringsAsFactors = FALSE)
  for (l in 1:length(csvfileslist)) {
    csvfileslist[[l]][sweep_variable] <- treatment_information[sweep_variable][treatment_information$Treatment_ID==csvfileslist[[l]]$Treatment_ID[1],]
    csvfileslist[[l]]$focal <- treatment_information$Focal.Species[treatment_information$Treatment_ID==csvfileslist[[l]]$Treatment_ID[1]]
    csvfileslist[[l]]$treatment_name <- treatment_information$Species[treatment_information$Treatment_ID==csvfileslist[[l]]$Treatment_ID[1]]
    csvfileslist[[l]] <- csvfileslist[[l]] %>% filter(Cumulative_time == treatment_information$Cumulative_time[treatment_information$Treatment_ID==csvfileslist[[l]]$Treatment_ID[1]])
  }
  df <- bind_rows(csvfileslist)
  return(df)
}

# time step data -------------------------------------------------------
df <- Parameter_sweep_load_data(treatments_filename = "timesteps", 
                                treatment_information = "Treatment_master_DHFR.csv", 
                                sweep_variable = "time_step")

df$treatment_name <- factor(df$treatment_name, 
                            levels = c("Ec1", "Ec1,Ec2,Ec1", "Ec1,Ec3,Ec1", "Ec2", "Ec2,Ec1,Ec2", "Ec2,Ec3,Ec2", "Ec3", "Ec3,Ec1,Ec3", "Ec3,Ec2,Ec3"),
                            labels = c("Ec1", "Ec1 \nwith \nEc2", "Ec1 \nwith \nEc3", "Ec2", "Ec2 \nwith \nEc1", "Ec2 \nwith \nEc3", "Ec3", "Ec3 \nwith \nEc1", "Ec3 \nwith \nEc2"))

my_comparisons_Ec <- list(c("Ec1", "Ec1 \nwith \nEc3"),
                          c("Ec1", "Ec1 \nwith \nEc2"))

my_comparisons_Kp <- list(c("Ec2", "Ec2 \nwith \nEc3"),
                          c("Ec2", "Ec2 \nwith \nEc1"))

my_comparisons_Se <- list(c("Ec3", "Ec3 \nwith \nEc2"),
                          c("Ec3", "Ec3 \nwith \nEc1"))

baseline <- 36

# plot -------------------------------------------------------
compare_means(Resistance ~ treatment_name, data = df %>% filter(focal == 'Ec1'), group.by = 'time_step', ref.group = "Ec1", method = "wilcox.test", p.adjust.method = "bonferroni")

Ec.max <- as.numeric(df %>% filter(focal == 'Ec1') %>% summarise(maxR = max(Resistance)))
E <- ggplot(data = df %>% filter(focal == 'Ec1'), 
            aes(x=treatment_name, y=Resistance, color = treatment_name))+
  stat_summary(data = df %>% filter(focal == 'Ec1', time_step !=baseline), mapping = aes(treatment_name, Resistance, color = treatment_name), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(data = df %>% filter(focal == 'Ec1', time_step !=baseline), width = point_jitter_width, size = point_size, alpha = 0.1) +
  stat_summary(data = df %>% filter(focal == 'Ec1', time_step ==baseline), color = 'dark grey', mapping = aes(treatment_name, Resistance, color = treatment_name), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(data = df %>% filter(focal == 'Ec1', time_step ==baseline), color = 'dark grey', width = point_jitter_width, size = point_size, alpha = 0.1) +
  facet_grid(~time_step)+
  stat_compare_means(comparisons = my_comparisons_Ec, method = "wilcox.test", aes(label = ..p.signif..), method.args = list(p.adjust.method = "bonferroni"), size = 3) +
  ylab(y_axis_label)+
  scale_y_continuous(limits = c(0, Ec.max+(Ec.max*.2)))+
  theme(axis.title.x = element_blank()) +
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text.y = element_text(size = text_size_axis_tick_labels))+
  theme(axis.text.x = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))+
  scale_color_manual(values = color_Ec, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank())

compare_means(Resistance ~ treatment_name, data = df %>% filter(focal == 'Ec2'), group.by = 'time_step', ref.group = "Ec2", method = "wilcox.test", p.adjust.method = "bonferroni")

Kp.max <- as.numeric(df %>% filter(focal == 'Ec2') %>% summarise(maxR = max(Resistance)))
K <- ggplot(data = df %>% filter(focal == 'Ec2'), 
            aes(x=treatment_name, y=Resistance, color = treatment_name))+
  stat_summary(data = df %>% filter(focal == 'Ec2', time_step !=baseline), mapping = aes(treatment_name, Resistance, color = treatment_name), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(data = df %>% filter(focal == 'Ec2', time_step !=baseline), width = point_jitter_width, size = point_size, alpha = 0.1) +
  stat_summary(data = df %>% filter(focal == 'Ec2', time_step ==baseline), color = 'dark grey', mapping = aes(treatment_name, Resistance, color = treatment_name), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(data = df %>% filter(focal == 'Ec2', time_step ==baseline), color = 'dark grey', width = point_jitter_width, size = point_size, alpha = 0.1) +
  facet_grid(~time_step)+
  stat_compare_means(comparisons = my_comparisons_Kp, method = "wilcox.test", aes(label = ..p.signif..), method.args = list(p.adjust.method = "bonferroni"), size = 3) +
  ylab(y_axis_label)+
  scale_y_continuous(limits = c(0, Kp.max+(Kp.max*.2)))+
  theme(axis.title.x = element_blank()) +
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text.y = element_text(size = text_size_axis_tick_labels))+
  theme(axis.text.x = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))+
  scale_color_manual(values = color_Kp, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank())

compare_means(Resistance ~ treatment_name, data = df %>% filter(focal == 'Ec3'), group.by = 'time_step', ref.group = "Ec3", method = "wilcox.test", p.adjust.method = "bonferroni")

Se.max <- as.numeric(df %>% filter(focal == 'Ec3') %>% summarise(maxR = max(Resistance)))
S <- ggplot(data = df %>% filter(focal == 'Ec3'), 
            aes(x=treatment_name, y=Resistance, color = treatment_name))+
  stat_summary(data = df %>% filter(focal == 'Ec3', time_step !=baseline), mapping = aes(treatment_name, Resistance, color = treatment_name), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(data = df %>% filter(focal == 'Ec3', time_step !=baseline), width = point_jitter_width, size = point_size, alpha = 0.1) +
  stat_summary(data = df %>% filter(focal == 'Ec3', time_step ==baseline), color = 'dark grey', mapping = aes(treatment_name, Resistance, color = treatment_name), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(data = df %>% filter(focal == 'Ec3', time_step ==baseline), color = 'dark grey', width = point_jitter_width, size = point_size, alpha = 0.1) +
  facet_grid(~time_step)+
  stat_compare_means(comparisons = my_comparisons_Se, method = "wilcox.test", aes(label = ..p.signif..), method.args = list(p.adjust.method = "bonferroni"), size = 3) +
  ylab(y_axis_label)+
  scale_y_continuous(limits = c(0, Se.max+(Se.max*.2)))+
  theme(axis.title.x = element_blank()) +
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text.y = element_text(size = text_size_axis_tick_labels))+
  theme(axis.text.x = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))+
  scale_color_manual(values = color_Se, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank())

Fig <- (E / K / S) +
  plot_annotation(tag_levels = list(c('a', 'b', 'c'))) & 
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = Figure_label_size)) 

save_plot("SI_Fig11.pdf", plot = Fig, base_width = 6.5, base_height = 7.5)
file.remove("Rplots.pdf")
