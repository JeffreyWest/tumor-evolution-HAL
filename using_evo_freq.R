

# install.packages('devtools')
# library(devtools)
# install_github('MathOnco/EvoFreq')

library(EvoFreq)
library(ggplot2)
library(colormap)
library(magick)
library(ggraph)
library(igraph)
library(reshape2)
library(plyr)
library(colorspace)
library(gridExtra)

# custom color maps
# custom_colors <- c("#3874b1", "#c6382c","#4f9f39", "#bdbe3a","#8e66ba","#f08627","#53bbce", "#d67bbf","#85584c", "#b2c5e6","#f39c97", "#a6de90","#dcdc93","#c2aed3","#f6bf7e","#a9d8e4","#eeb8d1","#be9d92","#c7c7c7","#7f7f7f")
# c_pallete <- custom_colors

setwd('~/Documents/GitHub/tumor-evolution-HAL/data-output/')
clone_history_file <- paste("./evofreq_dataframe.csv", sep = "")

clone_df <- read.csv(clone_history_file, check.names = F, header = T)
parent_list <- clone_df$Parent
attr_list <- clone_df$Drivers
clone_list <- row.names(clone_df)
time_pts <- as.numeric(colnames(clone_df))
time_pts <- which(!is.na(time_pts))

pos_df1 <- get_evofreq(clone_df[time_pts],clone_list,parent_list, clone_cmap = "rainbow_soft")
plot_evofreq(pos_df1)
