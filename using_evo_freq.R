## remove Evofreq
detach("package:EvoFreq", unload = TRUE) # May throw error, it is ok
remove.packages("EvoFreq") # May throw error, it is ok

## to install EvoFreq:
# install.packages('devtools')
library(devtools)
install_github('MathOnco/EvoFreq')
               # , ref = "HAL_integration", force=T)

library(colormap)
library(ggplot2)
# read.HAL <- function(path_to_file){
#   ### FOR TESTING ###
#   # path_to_file <- f
#   ###
#   clone_df <- read.csv(f, check.names = F, stringsAsFactors = F)  
#   df_cols <- colnames(clone_df)
#   time_col_idx <- suppressWarnings(which(!is.na(as.numeric(df_cols))))
#   time_cols <- df_cols[time_col_idx]
#   size_df <- clone_df[, time_cols]
#   
#   attribute_col_idx <- which(!df_cols %in% c(time_cols, "CloneID", "ParentID"))
#   if(length(attribute_col_idx) > 0){
#     attribute_df <- cbind(data.frame("CloneID"=clone_df$CloneID),clone_df[, attribute_col_idx]) 
#   }else{
#     attribute_df <- NULL
#   }
#   
#   return(list("clones"=clone_df$CloneID, "parents"=clone_df$ParentID, "size_df"=size_df, "attributes"=attribute_df)) 
# }



library(EvoFreq)
# library(ggplot2)
# library(colormap)

# custom color maps
driver_colors <- c("#3874b1", "#c6382c","#4f9f39", "#bdbe3a","#8e66ba","#f08627","#53bbce", "#d67bbf","#85584c", "#b2c5e6","#f39c97", "#a6de90","#dcdc93","#c2aed3","#f6bf7e","#a9d8e4","#eeb8d1","#be9d92","#c7c7c7","#7f7f7f")
passenger_colors <- colormap("density")

# set working directory, read in clone information output from HAL
setwd('~/Documents/GitHub/tumor-evolution-HAL/')
# clone_history_file <- paste("./evofreq_dataframe.csv", sep = "")
clone_history_file <- paste("./test_new.csv", sep = "")


clone_df <- read.csv(clone_history_file, check.names = F, header = T)
parent_list <- clone_df$Parent


# map the number of passengers/drivers to the appropriate color
color_by_passengers <- passenger_colors[clone_df$Passengers+1]
color_by_drivers <- driver_colors[clone_df$Drivers]

# get list of clones & list of time points
clone_list <- row.names(clone_df)
time_pts <- as.numeric(colnames(clone_df))
time_pts <- which(!is.na(time_pts))


evo_df <- get_evofreq(clone_df[time_pts],clone_list,parent_list, fill_value = color_by_drivers, threshold=0.0)
plot_evofreq(evo_df)

tree_info <- get_evogram(clone_df[time_pts],as.numeric(clone_list),as.numeric(parent_list), fill_value = color_by_drivers)
plot_evogram(tree_info$dendro_pos, tree_info$links)







# f <- "~/Documents/GitHub/tumor-evolution-HAL/data-output/evofreq_dataframe.csv"




# passengers <- hal_info$attributes$Passengers
# drivers <- hal_info$attributes$Drivers
# driver_colors <- c("#3874b1", "#c6382c","#4f9f39", "#bdbe3a","#8e66ba","#f08627","#53bbce", "#d67bbf","#85584c", "#b2c5e6","#f39c97", "#a6de90","#dcdc93","#c2aed3","#f6bf7e","#a9d8e4","#eeb8d1","#be9d92","#c7c7c7","#7f7f7f")
# passenger_colors <- colormap("density")


# hal_plot_df <- get_evofreq(size_df = hal_info$size_df, clones=hal_info$clones, parents=hal_info$parents, fill_value = passengers, threshold = 0, clone_cmap = "jet")

# color_by_passengers <- passenger_colors[passengers+1]
# color_by_drivers <- driver_colors[drivers]


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









