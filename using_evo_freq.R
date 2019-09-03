## to install EvoFreq:
# install.packages('devtools')
library(devtools)
install_github('MathOnco/EvoFreq')

library(EvoFreq)

# default coloring
mypath <- "~/Documents/GitHub/tumor-evolution-HAL/data-output/phylogeny_tracker.csv"
hal_info <- read.HAL(mypath)
print(hal_info$evofreq_plot)

# color by the number of passenger mutations
hal_info <- read.HAL(mypath, fill_name = "Passengers")
print(hal_info$evofreq_plot)

# color by the same color scheme as HAL
hal_info <- read.HAL(mypath,fill_name = "Color")
print(hal_info$evofreq_plot)

# alternatively, you can update the color scheme of an existing plot as follows:
updated_plot <- update_colors(hal_info$freq_frame, hal_info$clones, hal_info$attributes$Color)
plot_evofreq(updated_plot)
