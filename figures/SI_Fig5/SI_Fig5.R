# packages ####
library(tidyverse)
library(ggpubr)
library(cowplot)
library(egg)
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
# colors ####
red_Ec_color <- rgb(255, 48, 48, maxColorValue = 255)
blue_Kp_color <- rgb(0, 0, 238, maxColorValue = 255)
yellow_Se_color <- rgb(242, 184, 0, maxColorValue = 255)

# load data ####
MIC_df <- read.csv("data/MIC_data.csv")
barcode_df <- read.csv("data/MIC_simulation_file_3barcodes.csv")

# MIC analysis
MIC <- MIC_df %>%
  group_by(Species, Shorthand) %>%
  filter(!(is.na(MIC))) %>%
  summarise(N = n(),
            mean_MIC = mean(MIC),
            se_MIC = sd(MIC)/sqrt(N))%>%
  select(-N)

# Barcode analysis
barcode <- barcode_df %>%
  select(Species, Shorthand, fit.1, fit.2) %>%
  pivot_longer(cols = c(fit.1, fit.2), names_to = 'replicate') %>%
  group_by(Species, Shorthand) %>%
  summarise(N = n(),
            mean_barcode = mean(value),
            se_barcode = sd(value)/sqrt(N)) %>%
  select(-N)

total <- left_join(barcode, MIC)

Fig <- ggplot(data = total, aes(x = mean_barcode, y = mean_MIC, color = Species)) +
  geom_errorbar(aes(ymin=mean_MIC-se_MIC, ymax=mean_MIC+se_MIC), width=0) +
  geom_errorbar(aes(xmin=mean_barcode-se_barcode, xmax=mean_barcode+se_barcode), width=0) +
  geom_point(aes(color = Species), size = 1)+
  scale_x_continuous(trans = 'log2', limits = c(-Inf, 2048),
                     breaks = c(0.25, 4, 64, 1024), 
                     labels = c("0.25", "4", "64", "1024"))+
  scale_y_continuous(trans = 'log2', limits = c(-Inf, 2048),
                     breaks = c(0.25, 4, 64, 1024), 
                     labels = c("0.25", "4", "64", "1024"))+
  ylab('MIC resistance level') +
  xlab('pooled assay resistance level')+
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x)+
  stat_cor(method = "pearson", size = 2)+
  scale_color_manual(values = c(red_Ec_color, blue_Kp_color, yellow_Se_color))+ 
  guides(color = FALSE)
  
# assembly ####
p_width = 3
p_height = 3
gFig <- set_panel_size(Fig, width  = unit(p_width, "in"), height = unit(p_height, "in"))

save_plot("SI_Fig2.pdf", plot = gFig, base_width = 4, base_height = 4)

