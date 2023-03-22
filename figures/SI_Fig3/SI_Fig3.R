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
# General plotting customization -------------------------------------------------------
Figure_width = 3.25 # in width for one column figure in PNAS
#plot_width = 3
#plot_height = plot_width
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
df_Ec_only <- read.csv("inputs/76.csv")
df_Ec_Kp <- read.csv("inputs/77.csv")
df_Ec_Se <- read.csv("inputs/78.csv")
df_Ec_only$focal <- 'Ec'
df_Ec_Kp$focal <- 'Ec'
df_Ec_Se$focal <- 'Ec'

df_Kp_only <- read.csv("inputs/79.csv")
df_Kp_Ec <- read.csv("inputs/80.csv")
df_Kp_Se <- read.csv("inputs/81.csv")
df_Kp_only$focal <- 'Kp'
df_Kp_Ec$focal <- 'Kp'
df_Kp_Se$focal <- 'Kp'

df_Se_only <- read.csv("inputs/82.csv")
df_Se_Ec <- read.csv("inputs/83.csv")
df_Se_Kp <- read.csv("inputs/84.csv")
df_Se_only$focal <- 'Se'
df_Se_Ec$focal <- 'Se'
df_Se_Kp$focal <- 'Se'

df_rbind <- rbind(df_Ec_only, df_Ec_Kp, df_Ec_Se, df_Kp_only, df_Kp_Ec, df_Kp_Se, df_Se_only, df_Se_Ec, df_Se_Kp)

df <- df_rbind %>%
  filter(Cumulative_time == 60)

# Statistics -------------------------------------------------------
df$Treatment_ID <- factor(df$Treatment_ID, levels= c("76", "77", "78", "79", "80", "81", "82", "83", "84"), 
                          labels = c("Ec only", "Ec with \nHGT to \nKp", "Ec with \nHGT to \nSe",
                                     "Kp only", "Kp with \nHGT to \nEc", "Kp with \nHGT to \nSe",
                                     "Se only", "Se with \nHGT to \nEc", "Se with \nHGT to \nKp"))


color_Ec = c(red_Ec_color, purple_Ec.Kp_color, orange_Ec.Se_color)
color_Kp = c(blue_Kp_color, purple_Ec.Kp_color, green_Se.Kp_color)
color_Se = c(yellow_Se_color, orange_Ec.Se_color, green_Se.Kp_color)

my_comparisons_Ec <- list(c("Ec only", "Ec with \nHGT to \nSe"),
                          c("Ec only", "Ec with \nHGT to \nKp"))

my_comparisons_Kp <- list(c("Kp only", "Kp with \nHGT to \nSe"),
                          c("Kp only", "Kp with \nHGT to \nEc"))

my_comparisons_Se <- list(c("Se only", "Se with \nHGT to \nKp"),
                          c("Se only", "Se with \nHGT to \nEc"))

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
compare_means(Resistance ~ Treatment_ID, data = df %>% filter(focal == 'Ec'), ref.group = "Ec only", method = "wilcox.test", p.adjust.method = "bonferroni")
Ec_mean <- read.csv("inputs/76_averages.csv")
Ec.mean <- Ec_mean %>% filter(Cumulative_time == 60) %>% select(Resistance.mean)
E <- ggplot(data = df %>% filter(focal == 'Ec'), 
            aes(x=Treatment_ID, y=Resistance, color = Treatment_ID))+
  stat_summary(data = df %>% filter(focal == 'Ec'), mapping = aes(Treatment_ID, Resistance, color = Treatment_ID), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(width = point_jitter_width, size = point_size, alpha = 0.1) +
  stat_compare_means(comparisons = my_comparisons_Ec, method = "t.test", aes(label = ..p.signif..), method.args = list(p.adjust.method = "bonferroni")) +
  geom_hline(yintercept = Ec.mean$Resistance.mean, linetype = 'dashed', color = red_Ec_color)+
  ylab(y_axis_label)+
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

compare_means(Resistance ~ Treatment_ID, data = df %>% filter(focal == 'Kp'), ref.group = "Kp only", method = "wilcox.test", p.adjust.method = "bonferroni")
Kp_mean <- read.csv("inputs/79_averages.csv")
Kp.mean <- Kp_mean %>% filter(Cumulative_time == 60) %>% select(Resistance.mean)
K <- ggplot(data = df %>% filter(focal == 'Kp'), 
            aes(x=Treatment_ID, y=Resistance, color = Treatment_ID))+
  stat_summary(data = df %>% filter(focal == 'Kp'), mapping = aes(Treatment_ID, Resistance, color = Treatment_ID), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(width = point_jitter_width, size = point_size, alpha = 0.1) +
  stat_compare_means(comparisons = my_comparisons_Kp, method = "t.test", aes(label = ..p.signif..), method.args = list(p.adjust.method = "bonferroni")) +
  geom_hline(yintercept = Kp.mean$Resistance.mean, linetype = 'dashed', color = blue_Kp_color)+
  ylab(y_axis_label)+
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
  scale_color_manual(values = color_Kp, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank())

compare_means(Resistance ~ Treatment_ID, data = df %>% filter(focal == 'Se'), ref.group = "Se only", method = "wilcox.test", p.adjust.method = "bonferroni")
Se_mean <- read.csv("inputs/82_averages.csv")
Se.mean <- Se_mean %>% filter(Cumulative_time == 60) %>% select(Resistance.mean)
S <- ggplot(data = df %>% filter(focal == 'Se'), 
            aes(x=Treatment_ID, y=Resistance, color = Treatment_ID))+
  stat_summary(data = df %>% filter(focal == 'Se'), mapping = aes(Treatment_ID, Resistance, color = Treatment_ID), fun.data=MinMeanSEMMax, geom="boxplot", size = box_line_size)+ 
  geom_jitter(width = point_jitter_width, size = point_size, alpha = 0.1) +
  stat_compare_means(comparisons = my_comparisons_Se, method = "t.test", aes(label = ..p.signif..), method.args = list(p.adjust.method = "bonferroni")) +
  geom_hline(yintercept = Se.mean$Resistance.mean, linetype = 'dashed', color = yellow_Se_color)+
  ylab(y_axis_label)+
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
  scale_color_manual(values = color_Se, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank())

Fig <- (E + K + S) +
  plot_annotation(tag_levels = list(c('a', 'b', 'c'))) & 
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = Figure_label_size)) 

save_plot("SI_Fig3.pdf", plot = Fig, base_width = Figure_width*2, base_height = Figure_height)
file.remove("Rplots.pdf")





