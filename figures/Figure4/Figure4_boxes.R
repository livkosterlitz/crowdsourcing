# packages #####
library(tidyverse)
library(ggpubr)
library(cowplot)
library(egg)
library(RColorBrewer)
source("src/Random_Landscape_Functions.R")
theme_set(theme_cowplot() +
            theme(
              text = element_text(size = 9),
              axis.title = element_text(size = 9),
              axis.text = element_text(size = 7),
              legend.title = element_text(size = 9),
              legend.text = element_text(size = 9),
              plot.title = element_text(size = 9),
              plot.subtitle = element_text(size = 9),
              plot.caption = element_text(size = 9),
              strip.text = element_text(size = 9)
            )
)
# colors #####
red_Ec_color <- rgb(255, 48, 48, maxColorValue = 255)
blue_Kp_color <- rgb(0, 0, 238, maxColorValue = 255)
yellow_Se_color <- rgb(242, 184, 0, maxColorValue = 255)
purple_Ec.Kp_color <- rgb(178, 58, 238, maxColorValue = 255)
orange_Ec.Se_color <- rgb(238, 118, 0, maxColorValue = 255)
green_Se.Kp_color <- rgb(70, 178, 0, maxColorValue = 255)
color_Ec = c(red_Ec_color, purple_Ec.Kp_color, orange_Ec.Se_color)
color_Kp = c(blue_Kp_color, purple_Ec.Kp_color, green_Se.Kp_color)
color_Se = c(yellow_Se_color, orange_Ec.Se_color, green_Se.Kp_color)


df_Ec_only_box_lines <- read.csv("inputs/76_lines.csv")
lines <- ggplot(df_Ec_only_box_lines %>% filter(Replicate < 11), aes(x = Cumulative_time, y = Resistance, group = Replicate)) +
  geom_line(color = 'grey')+
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, by = 200))        
lines         
p_width = 1.3
p_height = 1.3
glines <- set_panel_size(lines, width  = unit(p_width, "in"), height = unit(p_height, "in"))

save_plot("Figure4_lines.pdf", plot = glines, base_width = 3, base_height = 3)



# load data
df_Ec_only_box <- read.csv("inputs/76.csv")
df_Ec_Kp_box <- read.csv("inputs/77.csv")
df_Ec_Se_box <- read.csv("inputs/78.csv")
df_Ec_only_box$focal <- 'Ec'
df_Ec_Kp_box$focal <- 'Ec'
df_Ec_Se_box$focal <- 'Ec'
df_Ec_rbind <- rbind(df_Ec_only_box, df_Ec_Kp_box, df_Ec_Se_box)
mapping <- c('76' = 'Ec', '77' = 'Ec-Kp', '78' ='Ec-Se')
df_Ec <- df_Ec_rbind %>%
  filter(Cumulative_time == 60) %>% 
  mutate(Treatment_ID = recode(Treatment_ID, !!!mapping)) %>%
  mutate(Treatment_ID = factor(Treatment_ID, levels = c('Ec', 'Ec-Kp', 'Ec-Se')))

Ec_mean <- read.csv("inputs/76_averages.csv")
Ec.mean <- Ec_mean %>% filter(Cumulative_time == 60) %>% select(Resistance.mean)
Ec <- ggplot(data = df_Ec, aes(x = Treatment_ID, y = Resistance, color = Treatment_ID)) +
  geom_violin(width = 0.9)+
  geom_jitter(width = 0.15, size = 1, alpha = 0.01)+
  stat_summary(fun = "mean", geom = "point", size = 2)+
  geom_hline(yintercept = Ec.mean[[1]], linetype = 'dashed', color = color_Ec[1])+
  theme(legend.position = "none") +
  ylab(" ")+
  xlab(" ")+
  scale_color_manual(values = color_Ec)+
  scale_y_continuous(limits = c(0, 2200), breaks = seq(0, 2200, by = 200))+
  theme(axis.title.y = element_blank())
Ec

df_Kp_only_box <- read.csv("inputs/79.csv")
df_Kp_Ec_box <- read.csv("inputs/80.csv")
df_Kp_Se_box <- read.csv("inputs/81.csv")
df_Kp_only_box$focal <- 'Kp'
df_Kp_Ec_box$focal <- 'Kp'
df_Kp_Se_box$focal <- 'Kp'
df_Kp_rbind <- rbind(df_Kp_only_box, df_Kp_Ec_box, df_Kp_Se_box)
mapping <- c('79' = 'Kp', '80' = 'Kp-Ec', '81' ='Kp-Se')
df_Kp <- df_Kp_rbind %>%
  filter(Cumulative_time == 60) %>% 
  mutate(Treatment_ID = recode(Treatment_ID, !!!mapping)) %>%
  mutate(Treatment_ID = factor(Treatment_ID, levels = c('Kp', 'Kp-Ec', 'Kp-Se')))

Kp_mean <- read.csv("inputs/79_averages.csv")
Kp.mean <- Kp_mean %>% filter(Cumulative_time == 60) %>% select(Resistance.mean)
Kp <- ggplot(data = df_Kp, aes(x = Treatment_ID, y = Resistance, color = Treatment_ID)) +
  geom_violin(width = 0.9)+
  geom_jitter(width = 0.15, size = 1, alpha = 0.01)+
  stat_summary(fun = "mean", geom = "point", size = 2)+
  geom_hline(yintercept = Kp.mean[[1]], linetype = 'dashed', color = color_Kp[1])+
  theme(legend.position = "none") +
  ylab(" ")+
  xlab(" ")+
  scale_color_manual(values = color_Kp)+
  scale_y_continuous(limits = c(0, 240), breaks = seq(0, 2000, by = 200))+
  theme(axis.title.y = element_blank())
Kp

df_Se_only_box <- read.csv("inputs/82.csv")
df_Se_Ec_box <- read.csv("inputs/83.csv")
df_Se_Kp_box <- read.csv("inputs/84.csv")
df_Se_only_box$focal <- 'Se'
df_Se_Ec_box$focal <- 'Se'
df_Se_Kp_box$focal <- 'Se'
mapping <- c('82' = 'Se', '83' = 'Se-Ec', '84' ='Se-Kp')
df_Se_rbind <- rbind(df_Se_only_box, df_Se_Ec_box, df_Se_Kp_box)
df_Se <- df_Se_rbind %>%
  filter(Cumulative_time == 60) %>% 
  mutate(Treatment_ID = recode(Treatment_ID, !!!mapping)) %>%
  mutate(Treatment_ID = factor(Treatment_ID, levels = c('Se', 'Se-Ec', 'Se-Kp')))

Se_mean <- read.csv("inputs/82_averages.csv")
Se.mean <- Se_mean %>% filter(Cumulative_time == 60) %>% select(Resistance.mean)
Se <- ggplot(data = df_Se, aes(x = Treatment_ID, y = Resistance, color = Treatment_ID)) +
  geom_violin(width = 0.9)+
  geom_jitter(width = 0.15, size = 1, alpha = 0.01)+
  stat_summary(fun = "mean", geom = "point", size = 2)+
  geom_hline(yintercept = Se.mean[[1]], linetype = 'dashed', color = color_Se[1])+
  theme(legend.position = "none") +
  ylab(" ")+
  xlab(" ")+
  scale_color_manual(values = color_Se)+
  scale_y_continuous(limits = c(0, 3350), breaks = seq(0, 4000, by = 200))+
  theme(axis.title.y = element_blank())
Se


# assembly ####
p_width = 0.9.3
p_height = 1.3
gEc <- set_panel_size(Ec, width  = unit(p_width, "in"), height = unit(p_height, "in"))
gKp <- set_panel_size(Kp, width  = unit(p_width, "in"), height = unit(p_height, "in"))
gSe <- set_panel_size(Se, width  = unit(p_width, "in"), height = unit(p_height, "in"))
gEKS <- plot_grid(gEc,gKp,gSe, nrow = 3)

save_plot("Figure4_boxes.pdf", plot = gEKS, base_width = 6.5, base_height = 6.5)
