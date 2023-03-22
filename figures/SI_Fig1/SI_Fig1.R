## Packages ----------------------------------------------------------------
library(drc)
library(tidyverse)
library(ggnewscale)
library(broom)
library(modelr)
library(patchwork)
library(cowplot)

## Load data ---------------------------------------------------------------
AppGrowthRate_data_all <- read.csv("inputs/Outliers.csv")
AppGrowthRate_data <- AppGrowthRate_data_all %>%
  filter(Shorthand %in% c('G', 'pGM', 'E'))

AppGrowthRate_AllData <- AppGrowthRate_data %>%
  select(Species, Genotype, Shorthand, Barcode, Concentration, Cmax, AppGrowthRate, Outlier)
AppGrowthRate_truncated <- AppGrowthRate_AllData %>%
  filter(Concentration <= Cmax) %>%
  select(Species, Genotype, Shorthand, Barcode, Outlier, Concentration, AppGrowthRate) %>%
  mutate(Concentration = ifelse(Concentration == 0, 0.01575, Concentration))

## Functions for 4 parameter model fitting and plotting  -------------------------------------------------------
drm.func.Se <- function(x) {
  drm(AppGrowthRate ~ Concentration,
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                 fixed = c(NA, 0.0194509, NA, NA)), data = x)}
predict.fun <- function(x) {
  add_predictions(data.frame(Concentration = c(seq(0.0156,1, 0.001),seq(1.001,515, 0.001))), x)}
coefs.fun <- function(x) {coef(x) %>% tidy}

Se_outlierRemoved <- AppGrowthRate_truncated %>%
  filter(Outlier == F) 

Se_trunc_4param <- Se_outlierRemoved %>%
  group_by(Shorthand, Barcode) %>%
  nest() %>%
  mutate(drmod = map(data, drm.func.Se),
         pred = map(drmod, predict.fun),
         coefs = map(drmod, coefs.fun))

## General plotting customization -------------------------------------------------------
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
point_jitter_width = 0.15
point_size = 1
margins<- c(0.1,0.1,0.1,0.1)
box_line_size = 0.75/(ggplot2::.pt*72.27/96) 

## Plots ------------------
Top <- ggplot() +
  geom_line(data = AppGrowthRate_truncated, aes(Concentration, AppGrowthRate, color = Genotype, group = Barcode, linetype = Outlier), linewidth = 0.75) +
  geom_point(data = AppGrowthRate_truncated, aes(Concentration, AppGrowthRate, color = Genotype), size = 1.5) +
  scale_x_continuous(trans = 'log2', limits = c(0.0156, AppGrowthRate_AllData$Cmax)) +
  scale_y_continuous(limits = c(-0.1, 0.6)) +
  theme_bw() +
  ylab(label = "Approximate growth rate")+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text.y = element_text(size = text_size_axis_tick_labels))+
  theme(axis.text.x = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))

Bottom <- ggplot() +
  geom_line(aes(Concentration, pred, color = Shorthand, group = Barcode), data = Se_trunc_4param %>% unnest(pred), size = 1) +
  geom_vline(aes(xintercept = x, color = Shorthand, group = Barcode), linetype = 5, data = Se_trunc_4param %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)")) +
  scale_x_continuous(trans = 'log2', limits = c(0.0156, AppGrowthRate_AllData$Cmax))+
  theme_bw()+
  scale_y_continuous(limits = c(-0.1, 0.6)) +
  theme(legend.position = "none")+
  ylab(label = "Approximate growth rate")+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text.y = element_text(size = text_size_axis_tick_labels))+
  theme(axis.text.x = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))

Fig <- (Top / Bottom) +
  plot_annotation(tag_levels = list(c('a', 'b'))) & 
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = Figure_label_size)) 

save_plot("SI_Fig1.pdf", plot = Fig, base_width = 6.5*0.75, base_height = 7.5*0.75)
file.remove("Rplots.pdf")