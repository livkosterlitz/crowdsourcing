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

# Arrange plot --------------
### Panel G
# General plotting customization -------------------------------------------------------
Figure_width = 3.25 
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


# Box plot customization -------------------------------------------------
box_line_size = 0.75/(ggplot2::.pt*72.27/96) 

# Load conjugation data -------------------------------------------------------
df_Ec1_only <- read.csv("inputs/74.csv")
df_Ec1_Ec2 <- read.csv("inputs/75.csv")
df_Ec1_Ec3 <- read.csv("inputs/76.csv")
df_Ec1_only$focal <- 'Ec1'
df_Ec1_Ec2$focal <- 'Ec1'
df_Ec1_Ec3$focal <- 'Ec1'

df_Ec2_only <- read.csv("inputs/92.csv")
df_Ec2_Ec1 <- read.csv("inputs/93.csv")
df_Ec2_Ec3 <- read.csv("inputs/94.csv")
df_Ec2_only$focal <- 'Ec2'
df_Ec2_Ec1$focal <- 'Ec2'
df_Ec2_Ec3$focal <- 'Ec2'

df_Ec3_only <- read.csv("inputs/110.csv")
df_Ec3_Ec1 <- read.csv("inputs/111.csv")
df_Ec3_Ec2 <- read.csv("inputs/112.csv")
df_Ec3_only$focal <- 'Ec3'
df_Ec3_Ec1$focal <- 'Ec3'
df_Ec3_Ec2$focal <- 'Ec3'


df_rbind <- rbind(df_Ec1_only, 
                  df_Ec1_Ec2, 
                  df_Ec1_Ec3, 
                  df_Ec2_only, 
                  df_Ec2_Ec1, 
                  df_Ec2_Ec3,
                  df_Ec3_only,
                  df_Ec3_Ec1,
                  df_Ec3_Ec2)

df <- df_rbind %>%
  filter(Cumulative_time == 36)

# Statistics -------------------------------------------------------
df$Treatment_ID <- factor(df$Treatment_ID, levels= c("1", "2", "3", "9", "10", "11", "17", "18", "19"), 
                          labels = c("Ec1 only", "Ec1 with \nHGT to \nEc2", "Ec1 with \nHGT to \nEc3",
                                     "Ec2 only", "Ec2 with \nHGT to \nEc1", "Ec2 with \nHGT to \nEc3",
                                     "Ec3 only", "Ec3 with \nHGT to \nEc1", "Ec3 with \nHGT to \nEc2"))

my_comparisons_Ec1 <- list(c("Ec1 only", "Ec1 with \nHGT to \nEc3"),
                           c("Ec1 only", "Ec1 with \nHGT to \nEc2"))

my_comparisons_Ec2 <- list(c("Ec2 only", "Ec2 with \nHGT to \nEc3"),
                           c("Ec2 only", "Ec2 with \nHGT to \nEc1"))

my_comparisons_Ec3 <- list(c("Ec3 only", "Ec3 with \nHGT to \nEc2"),
                           c("Ec3 only", "Ec3 with \nHGT to \nEc1"))

color_Ec1 = c(red_Ec_color, purple_Ec.Kp_color, orange_Ec.Se_color) 
color_Ec2 = c(blue_Kp_color, purple_Ec.Kp_color, green_Se.Kp_color) 
color_Ec3 = c(yellow_Se_color, orange_Ec.Se_color, green_Se.Kp_color)

MinMeanSEMMax <- function(x) {
  v <- c(min(x), 
         mean(x) - sd(x)/sqrt(length(x)), 
         mean(x), 
         mean(x) + sd(x)/sqrt(length(x)), 
         max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

# Plot -------------------------------------------------------
compare_means(Resistance ~ Treatment_ID, data = df %>% filter(focal == 'Ec1'), ref.group = "Ec1 only", method = "wilcox.test", p.adjust.method = "bonferroni")
Ec1_mean <- read.csv("inputs/74_averages.csv")
Ec1.mean <- Ec1_mean %>% filter(Cumulative_time == 36) %>% select(Resistance.mean)
Ec1 <- ggplot(data = df %>% filter(focal == 'Ec1'), 
            aes(x=Treatment_ID, y=Resistance, color = Treatment_ID))+
  stat_summary(data = df %>% filter(focal == 'Ec1'), mapping = aes(Treatment_ID, Resistance, color = Treatment_ID), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(width = point_jitter_width, size = point_size, alpha = 0.1) +
  #stat_compare_means(aes(label=..p.adj..), method = "wilcox.test", ref.group = "Ec1 only", method.args = list(p.adjust.method = "bonferroni")) +
  stat_compare_means(comparisons = my_comparisons_Ec1, method = "wilcox.test", aes(label = ..p.signif..)) +
  ylab(y_axis_label)+
  geom_hline(yintercept = Ec1.mean$Resistance.mean, linetype = 'dashed', color = red_Ec_color)+
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
  #scale_fill_manual(values = color)+
  scale_color_manual(values = color_Ec1, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank())
Ec1

compare_means(Resistance ~ Treatment_ID, data = df %>% filter(focal == 'Ec2'), ref.group = "Ec2 only", method = "wilcox.test", p.adjust.method = "bonferroni")
Ec2_mean <- read.csv("inputs/92_averages.csv")
Ec2.mean <- Ec2_mean %>% filter(Cumulative_time == 36) %>% select(Resistance.mean)
Ec2 <- ggplot(data = df %>% filter(focal == 'Ec2'), 
            aes(x=Treatment_ID, y=Resistance, color = Treatment_ID))+
  stat_summary(data = df %>% filter(focal == 'Ec2'), mapping = aes(Treatment_ID, Resistance, color = Treatment_ID), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(width = point_jitter_width, size = point_size, alpha = 0.1) +
  stat_compare_means(comparisons = my_comparisons_Ec2, method = "wilcox.test", aes(label = ..p.signif..)) +
  #stat_compare_means(aes(label=..p.adj..), method = "wilcox.test", ref.group = "Ec2 only", method.args = list(p.adjust.method = "bonferroni")) +
  ylab(y_axis_label)+
  geom_hline(yintercept = Ec2.mean$Resistance.mean, linetype = 'dashed', color = blue_Kp_color)+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text.y = element_text(size = text_size_axis_tick_labels))+
  theme(axis.text.x = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))+
  #scale_fill_manual(values = color)+
  scale_color_manual(values = color_Ec2, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank())
Ec2

compare_means(Resistance ~ Treatment_ID, data = df %>% filter(focal == 'Ec3'), ref.group = "Ec3 only", method = "wilcox.test", p.adjust.method = "bonferroni")
Ec3_mean <- read.csv("inputs/110_averages.csv")
Ec3.mean <- Ec3_mean %>% filter(Cumulative_time == 36) %>% select(Resistance.mean)
Ec3 <- ggplot(data = df %>% filter(focal == 'Ec3'), 
              aes(x=Treatment_ID, y=Resistance, color = Treatment_ID))+
  stat_summary(data = df %>% filter(focal == 'Ec3'), mapping = aes(Treatment_ID, Resistance, color = Treatment_ID), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(width = point_jitter_width, size = point_size, alpha = 0.1) +
  stat_compare_means(comparisons = my_comparisons_Ec3, method = "wilcox.test", aes(label = ..p.signif..)) +
  #stat_compare_means(aes(label=..p.adj..), method = "wilcox.test", ref.group = "Ec2 only", method.args = list(p.adjust.method = "bonferroni")) +
  ylab(y_axis_label)+
  geom_hline(yintercept = Ec3.mean$Resistance.mean, linetype = 'dashed', color = yellow_Se_color)+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text.y = element_text(size = text_size_axis_tick_labels))+
  theme(axis.text.x = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))+
  scale_color_manual(values = color_Ec3, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank())
Ec3


Fig <- (Ec1 + Ec2 + Ec3) +
  plot_annotation(tag_levels = list(c('a', 'b', 'c'))) & 
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = Figure_label_size)) 

save_plot("SI_Fig6.pdf", plot = Fig, base_width = Figure_width*2, base_height = Figure_height)
file.remove("Rplots.pdf")



