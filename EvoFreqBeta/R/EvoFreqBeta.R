# Created By Chandler Gatenbee And Ryan Schenck
# 19 July 2018
# Distributed under the
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# library(ggplot2)
# library(colormap)
# library(magick)
# library(ggraph)
# library(igraph)
# library(reshape2)

get_children <- function(c, clone_list, parent_list){
  children_idx <- which(parent_list==c)
  children <- clone_list[children_idx]
  return(children)
}

traverse <- function(c, clones, parents , fnx){
  f <- match.fun(fnx)
  f(c, clones)
  for(n in get_children(c, clones, parents)){
    traverse(n, clones, parents, fnx)
  }
}

get_all_idx <- function(c, c_list, p_list){
  all_children_idx <- c()
  add_child_idx <- function(i, c_list){
    child_idx <- which(c_list==i)
    all_children_idx <<- append(all_children_idx, child_idx)
  }
  traverse(c, c_list, p_list, add_child_idx)

  return(all_children_idx)
}

check_edges <- function(clones, parents){
  ##### FOR TESTING ###
  # clones <- clone_list
  # parents <- parent_list
  ######

  for(cid in clones){
    clone_idx <- which(clones==cid)
    parent <- parents[clone_idx]
    unique_parents <- unique(parent)
    n_unique_parents <- length(unique_parents)
    if(n_unique_parents != 1){
      stop(paste("clone", cid, "has", n_unique_parents, "parent(s), but it should only have 1"))
    }
  }
}

check_freq_mat <- function(freq_mat, clones, parents){
  ### FOR TESTING ###
  # freq_mat <- mutation_df
  # clones <- clone_list
  # parents <- parent_list
  #########
  for(i in seq(1, length(clones))){
    cln <- clones[i]
    cln_freq <- freq_mat[i, ]
    all_children_idx <- get_all_idx(cln, clones, parents)
    children <- clones[all_children_idx]
    children_idx <- all_children_idx[clones[all_children_idx] != cln]
    if(length(children == 0)){
      next
    }
    children_df <- freq_mat[children_idx, ]
    parent_greater_than_children <- apply(children_df, 1, function(x){all(cln_freq>=x)})
    if(all(parent_greater_than_children)==F){
      too_big_idx <- which(parent_greater_than_children==F)
      clones_too_big <- clones[too_big_idx]
      stop(paste(cln, "has descendent mutations greater than that are greater than it's size. Descendents are:", clones_too_big))
    }
  }

}

get_mutation_freq_df <- function(clone_size_df, clones, parents){
  ### FOR TESTING ####
  # clone_size_df <- clone_df
  # clones <- clone_list
  # parents <- parent_list
  ####

  check_edges(clones, parents)

  greater_than_one_time_pt <- ncol(clone_size_df)>1
  mutation_df <- matrix(NA, nrow = nrow(clone_size_df), ncol = ncol(clone_size_df))
  for(i in seq(1, length(clones))){
    cln <- clones[i]
    all_children_idx <- get_all_idx(cln, clones, parents)

    child_sizes_df <- clone_size_df[all_children_idx, ]
    if(greater_than_one_time_pt){
      total_cells_carrying_mutation <- colSums(child_sizes_df) ### recursion already includes adding the current clone
    }else{
      total_cells_carrying_mutation <- sum(child_sizes_df)
    }

    mutation_df[i, ] <- total_cells_carrying_mutation

  }

  # check_freq_mat(mutation_df, clones, parents)
  return(mutation_df)
}

order_clones_at_time <- function(clones_at_time, parents, clones, root_id){
  ### possible that parent and child arose within the same period. Need to ensure that parents are drawn first

  #### For testing ###
  # parents <- filtered_parents
  # clones <- filtered_clones
  # clones_at_time <- c(0.2707922, 1.0000000)
  # root_id <-1
  # root_id <-0.2707922
  ####

  n_clones <- length(clones_at_time)
  n_ancestors_in_t <- rep(-1, n_clones)
  for(ct_i in seq(n_clones)){
    ancestors_at_time <- 0
    ct <- clones_at_time[ct_i]
    clone_idx <- which(clones==ct)
    if(length(clone_idx)==0){
      next
    }
    ct_parent <- parents[clone_idx]
    while(ct_parent != root_id){
      if(ct_parent %in% clones_at_time){
        ancestors_at_time <- ancestors_at_time + 1
      }
      clone_idx <- which(clones==ct_parent)
      ct_parent <- parents[clone_idx]
      if(length(clone_idx)==0){
        break
      }
    }
    n_ancestors_in_t[ct_i] <- ancestors_at_time
  }
  new_order <- order(n_ancestors_in_t, decreasing = F)
  clones_at_time_ordered <- clones_at_time[new_order]

  return(clones_at_time_ordered)


}


predict_sizes <- function(x1, y1, x2, y2, new_x){
  ### FOR Testing ###
  # y1 <- 1.00
  # y2 <- 0.02
  # x1 <- 0
  # x2 <- 30
  # new_x <- 15
  ####


  fit <- lm(c(y1, y2) ~ c(x1, x2))
  new_y <- fit$coefficients[[1]] + fit$coefficients[[2]]*new_x
  return(new_y)
}



interp_mut_mat <- function(mut_mat, time_pts= NULL){
  #For each clone, estimate intermediate sizes

  ### TODO FOR TESTING ###
  # mut_mat <- clone_df
  # time_pts <- timepoints
  ####

  if(is.null(time_pts)){
    time_pts <- seq(1, ncol(mut_mat))
  }


  intermediate_vals <- mapply(i=seq(2, length(time_pts)), FUN=function(i){(time_pts[i]-time_pts[i-1])/2 + time_pts[i-1]})
  new_times <- rep(NA, length(time_pts)+length(intermediate_vals))
  og_time_idx <- seq(1, length(new_times),2)
  new_time_idx <- seq(2, length(new_times)-1,2)

  new_times[og_time_idx] <- time_pts
  new_times[new_time_idx] <- intermediate_vals
  n_new_times <- length(new_times)
  ## Increase width of matrix
  new_mat <- matrix(NA, nrow=nrow(mut_mat), ncol=n_new_times)
  colnames(new_mat) <- new_times
  row.names(new_mat) <- row.names(mut_mat)

  for(cidx in seq(nrow(new_mat))){
    # cidx <- 2
    clone_sizes <- mut_mat[cidx, ]
    intermediate_sizes <- mapply(i=seq(2, length(time_pts)), function(i){
      if(clone_sizes[i]== 0 | clone_sizes[i-1]== 0){
        return(0)
      }

      fit <- lm(c(clone_sizes[i], clone_sizes[i-1]) ~ c(time_pts[i], time_pts[i-1]))
      new_y <- fit$coefficients[[1]] + fit$coefficients[[2]]* intermediate_vals[i-1]
      return(new_y)
    })

    new_mat[cidx, og_time_idx] <- clone_sizes
    new_mat[cidx, new_time_idx] <- intermediate_sizes

    ### Remove points that were predicted immediately before creation and after extinction
    # zero_idx <- which(new_mat[cidx,]==0)
    #
    # if(length(zero_idx) > 0){
    #   before_creation_idx <- zero_idx[1] + 1
    #   after_extinction_idx <- zero_idx[length(zero_idx)] - 1
    #   if(before_creation_idx > 1){
    #     new_mat[cidx, before_creation_idx] <- 0
    #   }
    #   if(after_extinction_idx < n_new_times){
    #     new_mat[cidx, after_extinction_idx] <- 0
    #   }
    # }

  }

  return(new_mat)

}


get_pos <- function(clones, parents, mut_mat, cmap_vals, og_time_pts=NULL){
  ### FOR TESTING ###
  # parents <- parent_list
  # clones <- clone_list
  # # mut_mat <- clone_df
  # mut_mat <- freq_mat
  # # mut_mat <- clone_df
  # cmap_vals <- rep(1, length(clone_list))
  # color_vals <- filtered_info$color_vals
  # og_time_pts <- timepoints
  ###############
  roots_root_idx <- which(!parents %in% clones)
  roots_root <- parents[roots_root_idx]
  root_id <- clones[roots_root_idx]


  if(is.null(og_time_pts)){
    og_time_pts <- seq(1, ncol(mut_mat))
  }

  og_mut_mat <- mut_mat
  mut_mat <- interp_mut_mat(og_mut_mat, og_time_pts)
  time_pts <- as.numeric(colnames(mut_mat))


  origin_times <- apply(mut_mat, 1, function(x){which(x>0)[1]})
  unique_ordered_origin_times <- sort(unique(origin_times))
  y_btm_list <- list()
  y_top_list <- list()

  clone_id <- roots_root
  clone_str_id <- as.character(clone_id)
  n_time_pts <- ncol(mut_mat)

  max_size <- max(mut_mat)
  y_btm <- rep(0, n_time_pts)
  y_top <- rep(1, n_time_pts)
  y_btm_list[[clone_str_id]] <- y_btm
  y_top_list[[clone_str_id]] <- y_top

  # if(is.null(time_pts)){
  # x_vals <- seq(1, n_time_pts)
  # }else{
  x_vals <- time_pts

  clone_x <- c(x_vals, rev(x_vals))

  clone_str_id <- as.character(clone_id)
  clone_pos_list <- list()
  draw_order <- 0

  # clone_tree_df_list <- list()
  for(ot in unique_ordered_origin_times){
    parents_at_time <- unique(parents[origin_times==ot])
    n_parents <- length(parents_at_time)
    if(n_parents > 1){
      parents_at_time <- order_clones_at_time(parents_at_time, parents = parents, clones=clones, root_id=root_id)
    }
    for(parent in parents_at_time){
      parent_str_id <- as.character(parent)
      children_idx <- which(parents==parent)
      if(length(children_idx)==1){
        total_child_area <- mut_mat[children_idx,]
      }else{
        total_child_area <- as.numeric(colSums(mut_mat[children_idx,]))
      }

      if(parent==roots_root){
        parent_area <- rep(1, n_time_pts)
        n_children <- 1
      }else{
        parent_idx <- which(clones==parent)
        parent_area <- as.numeric(mut_mat[parent_idx,])
        # children_sizes <- mut_mat[children_idx,]
        # if(length(children_idx) > 1){
        #   # print(children_sizes)
        #   n_children <- apply(children_sizes, MARGIN = 2, function(x){length(x[x>0])})
        # }else{
        #   n_children <- rep(0, length(children_sizes))
        #   n_children[children_sizes > 0] <- 1
        #   print(n_children)
        # }
        # # print(n_children)

      }

      n_children <- length(children_idx)
      spacing <- (parent_area - total_child_area)/(n_children+1)


      current_btm <- y_btm_list[[parent_str_id]]

      loc <- spacing + current_btm
      for(cidx in children_idx){
        clone_id <-clones[cidx]
        clone_id_str <- as.character(clone_id)

        if(clone_id_str %in% names(clone_pos_list)){
          next
        }
        clone_color <- cmap_vals[cidx]

        child_size <- mut_mat[cidx,]

        non_zero_idx <- which(child_size > 0)
        # print(length)

        ### Allows one to see clones that existed for only 1 time point
        if(length(non_zero_idx) == 1){
          end_time_idx <- max(non_zero_idx) + 1
          start_time_idx <- min(non_zero_idx) - 1
          if(end_time_idx < n_time_pts){
            non_zero_idx <- c(non_zero_idx, end_time_idx)
          }
          ### Allows one to see that the population was not present in the previous time step
          if(start_time_idx > 0){
            non_zero_idx <- c(start_time_idx, non_zero_idx)
          }

        }

        bottom <- loc
        loc <- loc + child_size
        top <- loc
        loc <- loc + spacing

        y_btm_list[[clone_id_str]] <- bottom
        y_top_list[[clone_id_str]] <- top

        clone_x <- c(x_vals[non_zero_idx], rev(x_vals[non_zero_idx]))
        clone_y <- c(bottom[non_zero_idx], rev(top[non_zero_idx]))
        clone_sizes <- c(child_size[non_zero_idx], rev(child_size[non_zero_idx]))
        clone_shape_df <- data.frame(x=clone_x, y=clone_y, color=clone_color, clone_id = clone_id, parent= parent, origin_time=ot, draw_order = draw_order, size=clone_sizes)
        draw_order <- draw_order + 1
        clone_pos_list[[clone_id_str]] <- clone_shape_df

        # clone_tree_df <- data.frame(x= c(x_vals[non_zero_idx[1]], x_vals[tail(non_zero_idx, n=1)]), y=c(bottom[non_zero_idx[1]], bottom[tail(non_zero_idx, n=1)]), clone_id = clone_id, color=clone_color)
        # clone_tree_df_list[[clone_id_str]] <- clone_tree_df
      }
    }
  }

  clone_pos_df <- do.call(rbind, clone_pos_list)
  clone_pos_df$draw_order <- factor(clone_pos_df$draw_order, ordered=T, levels = sort(unique(clone_pos_df$draw_order), decreasing = F))
  clone_pos_df$color <- as.character(clone_pos_df$color)
  clone_pos_df <- clone_pos_df[order(clone_pos_df$draw_order),]
  clone_pos_df$x <- as.numeric(clone_pos_df$x)
  clone_pos_df$y <- as.numeric(clone_pos_df$y)

  # all_clone_tree_df <- do.call(rbind, clone_tree_df_list)

  # p <- ggplot(clone_pos_df, aes(x=x, y=y, fill=color, group=draw_order)) +
  # geom_polygon(size=0.1, color="grey75") +
  # scale_fill_identity() +
  # facet_wrap(~clone_id) +
  # theme_classic()

  # tree_p <- ggplot(all_clone_tree_df, aes(x=x, y=y, group=clone_id)) +
  # geom_line()

  return(clone_pos_df)

}

# get_pos <- function(clones, parents, mut_mat, cmap_vals, time_pts=NULL){
#   ### FOR TESTING ###
#   # parents <- parent_list
#   # clones <- clone_list
#   # mut_mat <- clone_df
#   # mut_mat <- freq_mat
#   # color_vals <- filtered_info$color_vals
#   ###############
#   roots_root_idx <- which(!parents %in% clones)
#   roots_root <- parents[roots_root_idx]
#   root_id <- clones[roots_root_idx]
#
#   origin_times <- apply(mut_mat, 1, function(x){which(x>0)[1]})
#   unique_ordered_origin_times <- sort(unique(origin_times))
#   y_btm_list <- list()
#   y_top_list <- list()
#
#   clone_id <- roots_root
#   clone_str_id <- as.character(clone_id)
#   n_time_pts <- ncol(mut_mat)
#
#   max_size <- max(mut_mat)
#   y_btm <- rep(0, n_time_pts)
#   y_top <- rep(1, n_time_pts)
#   y_btm_list[[clone_str_id]] <- y_btm
#   y_top_list[[clone_str_id]] <- y_top
#
#   if(is.null(time_pts)){
#     x_vals <- seq(1, n_time_pts)
#   }else{
#     x_vals <- time_pts
#   }
#
#   clone_x <- c(x_vals, rev(x_vals))
#
#   clone_str_id <- as.character(clone_id)
#   clone_pos_list <- list()
#   draw_order <- 0
#   for(ot in unique_ordered_origin_times){
#     parents_at_time <- unique(parents[origin_times==ot])
#     n_parents <- length(parents_at_time)
#     if(n_parents > 1){
#       parents_at_time <- order_clones_at_time(parents_at_time, parents = parents, clones=clones, root_id=root_id)
#     }
#     for(parent in parents_at_time){
#       parent_str_id <- as.character(parent)
#       children_idx <- which(parents==parent)
#       if(length(children_idx)==1){
#         total_child_area <- mut_mat[children_idx,]
#       }else{
#         total_child_area <- as.numeric(colSums(mut_mat[children_idx,]))
#       }
#
#       if(parent==roots_root){
#         parent_area <- rep(1, n_time_pts)
#         n_children <- 1
#       }else{
#         parent_idx <- which(clones==parent)
#         parent_area <- as.numeric(mut_mat[parent_idx,])
#         # children_sizes <- mut_mat[children_idx,]
#         # if(length(children_idx) > 1){
#         #   # print(children_sizes)
#         #   n_children <- apply(children_sizes, MARGIN = 2, function(x){length(x[x>0])})
#         # }else{
#         #   n_children <- rep(0, length(children_sizes))
#         #   n_children[children_sizes > 0] <- 1
#         #   print(n_children)
#         # }
#         # # print(n_children)
#
#       }
#
#       n_children <- length(children_idx)
#       spacing <- (parent_area - total_child_area)/(n_children+1)
#
#
#       current_btm <- y_btm_list[[parent_str_id]]
#
#       loc <- spacing + current_btm
#       for(cidx in children_idx){
#         clone_id <-clones[cidx]
#         clone_id_str <- as.character(clone_id)
#
#         if(clone_id_str %in% names(clone_pos_list)){
#           next
#         }
#         clone_color <- cmap_vals[cidx]
#
#         child_size <- mut_mat[cidx,]
#
#         non_zero_idx <- which(child_size > 0)
#
#         ### Allows one to see clones that existed for only 1 time point
#         if(non_zero_idx == 1){
#           end_time_idx <- max(non_zero_idx) + 1
#           start_time_idx <- min(non_zero_idx) - 1
#           if(end_time_idx < n_time_pts){
#             non_zero_idx <- c(non_zero_idx, end_time_idx)
#           }
#           ### Allows one to see that the population was not present in the previous time step
#           if(start_time_idx > 0){
#             non_zero_idx <- c(start_time_idx, non_zero_idx)
#           }
#
#         }
#
#         bottom <- loc
#         loc <- loc + child_size
#         top <- loc
#         loc <- loc + spacing
#
#         y_btm_list[[clone_id_str]] <- bottom
#         y_top_list[[clone_id_str]] <- top
#
#         clone_x <- c(x_vals[non_zero_idx], rev(x_vals[non_zero_idx]))
#         clone_y <- c(bottom[non_zero_idx], rev(top[non_zero_idx]))
#         clone_sizes <- c(child_size[non_zero_idx], rev(child_size[non_zero_idx]))
#
#         clone_shape_df <- data.frame(x=clone_x, y=clone_y, color=clone_color, clone_id = clone_id, parent= parent, origin_time=ot, draw_order = draw_order, size=clone_sizes)
#         draw_order <- draw_order + 1
#         clone_pos_list[[clone_id_str]] <- clone_shape_df
#       }
#     }
#   }
#
#   clone_pos_df <- do.call(rbind, clone_pos_list)
#
#   clone_pos_df$draw_order <- factor(clone_pos_df$draw_order, ordered=T, levels = sort(unique(clone_pos_df$draw_order), decreasing = F))
#   clone_pos_df$color <- as.character(clone_pos_df$color)
#   clone_pos_df <- clone_pos_df[order(clone_pos_df$draw_order),]
#
#   # p <- ggplot(clone_pos_df, aes(x=x, y=y, fill=color, group=draw_order)) +
#   #   geom_polygon(size=0.1, color="grey75") +
#   #   # scale_fill_identity() +
#   #   # facet_wrap(~clone_id) +
#   #   theme_classic()
#
#   return(clone_pos_df)
#
# }


filter_freq_mat <- function(clones, parents, freq_mat, color_vals, threshold){
  max_freqs <- apply(freq_mat, 1, max)
  above_thresh_idx <- which(max_freqs > threshold)

  filtered_freq_mat <- freq_mat[above_thresh_idx, ]
  filtered_clones <- clones[above_thresh_idx]
  filtered_parents <- parents[above_thresh_idx]
  if(!is.null(color_vals)){
    filtered_color_vals <- color_vals[above_thresh_idx]
  }else{
    filtered_color_vals <- NULL
  }

  return_vals <- list('freq_mat'=filtered_freq_mat, 'clones'= filtered_clones, 'parents'=filtered_parents, 'color_vals'=filtered_color_vals, 'thresh_idx'=above_thresh_idx)
}

get_freq_dynamics <- function(size_df, clones, parents, time_pts=NULL, attribute_vals=NULL, attribute_cmap='viridis', clone_cmap='rainbow_soft', threshold=0.01, scale_by_sizes_at_time = FALSE, data_type="size", interpolation_steps = 0){
  ### FOR TESTING ###
  # clone_cmap <- "rainbow_soft"
  # size_df <- clone_df
  # threshold <- 0.0
  # clones <- clone_list
  # parents <- parent_list
  ####
  if(data_type=="size"){
    ### data are sizes of each clone. Here, the frequency of each mutation will be determined
    freq_df <- get_mutation_freq_df(size_df, clones = clones, parents = parents)
    freq_mat <- as.matrix(freq_df)

  }else{
    check_edges(clones, parents)
    check_freq_mat(size_df, clones, parents)
    freq_mat <- size_df

  }

  if(scale_by_sizes_at_time){
    max_sizes_at_each_time <- apply(freq_mat, 2, max)
    freq_mat <- sweep(freq_mat, MARGIN = 2, max_sizes_at_each_time, FUN = "/")
    # freq_mat <- apply(freq_mat, MARGIN = 1, function(x){x/max_sizes_at_each_time})
    # freq_mat <- t(freq_mat)
  }else{
    freq_mat <- freq_mat/max(freq_mat)
  }

  ### Get color values
  if(is.null(attribute_vals)){
    if(clone_cmap[1] %in% names(colormap::colormaps)){
      cmap_vals <- get_clone_color(length(clones), cmap = clone_cmap)
    }else{
      cmap_vals <- clone_cmap
    }
  }else{
    cmap_vals <- get_attribute_colors(attribute_vals, cmap=attribute_cmap)
  }

  ### Filter values based on threshold
  filtered_info <- filter_freq_mat(clones, parents, freq_mat, cmap_vals, threshold)
  filtered_freq_mat <- filtered_info$freq_mat
  filtered_clones <- filtered_info$clones
  filtered_parents <- filtered_info$parents
  filtered_colors <- filtered_info$color_vals
  row.names(filtered_freq_mat) <- filtered_clones

  plot_pos_df <- get_pos(filtered_clones, filtered_parents, filtered_freq_mat, filtered_colors, time_pts)
  if(interpolation_steps > 0){
    plot_pos_df <- smooth_pos(plot_pos_df, n_intermediate_steps = interpolation_steps)
  }

  return(plot_pos_df)

}

get_clone_color <- function(n_clones, cmap="rainbow_soft"){
  ### FOR TESTING ###
  # n_clones <- length(clones)
  ####
  # n_shades <- n_clones
  pop_colors <-  colormap::colormap(colormap=colormap::colormaps[cmap][[1]], nshades=n_clones+1)[1:n_clones]
  # print(pop_colors)
  # ramp_colors <-  colormap(colormap=colormaps[cmap][[1]], nshades=n_clones)
  # n_ramp_reps <- ceiling(n_clones/n_shades)
  # pop_colors <- rep(ramp_colors, n_ramp_reps)[1:n_clones]
  # pop_colors <-sample(ramp_colors, n_clones, replace = T)
  return(pop_colors)
}

get_attribute_colors <- function(x, min_x =NULL, max_x = NULL, n_color_bins = 100, cmap="viridis"){
  ### FOR TESTING ###
  # x <- clone_antigenicity
  # root_idx <- which(clone_antigenicity==0)
  ####
  if(is.null(min_x)){
    min_x <- min(x)
  }
  if(is.null(max_x)){
    max_x <- max(x)
  }

  bin_breaks <- seq(min_x, max_x , length.out = n_color_bins)
  bin_number <- findInterval(x, bin_breaks, rightmost.closed = F)
  cmap_colors <- colormap::colormap(colormap::colormaps[cmap][[1]], nshades=n_color_bins)
  colors <- cmap_colors[bin_number]
  return(colors)
}

smooth_pos <- function(sparse_pos_df, n_intermediate_steps=5){
  interp_df_list <- list()
  unique_clones <- unique(sparse_pos_df$clone_id)
  for(cid in unique_clones){
    clone_pos_df <- subset(sparse_pos_df, clone_id==cid)
    nx <- length(unique(clone_pos_df$x))

    if(nx >= 4){
      forward_idx <- which(duplicated(clone_pos_df$x)==F)
      forward_df <- clone_pos_df[forward_idx,]
      ### need at least 4 points for spline interpolation
      new_x <- seq(min(forward_df$x), max(forward_df$x), length.out = length(forward_df$x)*n_intermediate_steps)

      # forward_spline <- smooth.spline(forward_df$x, forward_df$y)
      # new_forward <- predict(forward_spline, new_x)
      # new_forward_df <- data.frame('x'=new_forward$x, 'y'=new_forward$y, 'color'=clone_pos_df$color[1], 'clone_id'=clone_pos_df$clone_id[1], 'parent'=clone_pos_df$parent[1], 'origin_time'=clone_pos_df$origin_time[1], 'draw_order'=clone_pos_df$draw_order[1])

      func = splinefun(x=forward_df$x, y=forward_df$y, method="monoH.FC",  ties = mean)
      new_forward <- func(new_x)
      new_forward_df <- data.frame('x'=new_x, 'y'=new_forward, 'color'=clone_pos_df$color[1], 'clone_id'=clone_pos_df$clone_id[1], 'parent'=clone_pos_df$parent[1], 'origin_time'=clone_pos_df$origin_time[1], 'draw_order'=clone_pos_df$draw_order[1])



      reverse_idx <- which(duplicated(clone_pos_df$x)==T)
      reverse_df <- clone_pos_df[reverse_idx,]

      # reverse_spline <- smooth.spline(rev(reverse_df$x), rev(reverse_df$y))
      # new_reverse <- predict(reverse_spline, new_x)
      # new_reverse_df <- data.frame('x'=rev(new_reverse$x), 'y'=rev(new_reverse$y), 'color'=clone_pos_df$color[1], 'clone_id'=clone_pos_df$clone_id[1], 'parent'=clone_pos_df$parent[1], 'origin_time'=clone_pos_df$origin_time[1], 'draw_order'=clone_pos_df$draw_order[1])

      rev_func = splinefun(x=rev(reverse_df$x), y=rev(reverse_df$y), method="monoH.FC",  ties = mean)
      new_reverse <- rev_func(new_x)

      new_reverse_df <- data.frame('x'=rev(new_x), 'y'=rev(new_reverse), 'color'=clone_pos_df$color[1], 'clone_id'=clone_pos_df$clone_id[1], 'parent'=clone_pos_df$parent[1], 'origin_time'=clone_pos_df$origin_time[1], 'draw_order'=clone_pos_df$draw_order[1])
      new_clone_pos_df <- rbind(new_forward_df, new_reverse_df)
    }else{
      new_clone_pos_df <- clone_pos_df[colnames(clone_pos_df) != "size"]
    }
    interp_df_list[[cid]] <- new_clone_pos_df

  }
  interp_df <- do.call(rbind, interp_df_list)

  ### interpolation may cause points to go out of bounds
  interp_df$y[interp_df$y < 0] <- 0
  interp_df$y[interp_df$y > 1] <- 1

  return(interp_df)

}

plot_freq_dynamics <- function(pos_df, n_time_pts=NULL, start_time=NULL, end_time=NULL, bw=0.05, bc="grey75"){
  ### FOR TESTING ###
  # n_time_pts <- 20
  # start_time=NULL
  # end_time=NULL
  # bw=0.05
  #####

  unique_time_pts <- unique(pos_df$x)
  if(is.null(n_time_pts)){
    n_time_pts <- length(unique_time_pts)
  }
  if(is.null(start_time)){
    start_time <- 0
  }
  if(is.null(end_time)){
    end_time <- max(pos_df$x)
  }

  view_time_pts <- unique_time_pts[unique_time_pts <= end_time]
  time_pt_idx <- unique(floor(seq(1, length(view_time_pts), length.out = n_time_pts)))

  time_pts <- view_time_pts[time_pt_idx]
  if(!end_time %in% time_pts){
    time_pts <- c(time_pts, end_time)
  }

  view_df <- pos_df[pos_df$x %in% time_pts, ]

  ggevodyn <- ggplot2::ggplot(view_df, ggplot2::aes(x=x, y=y, group=draw_order, fill=color)) +
    ggplot2::geom_polygon(size=bw, color=bc) +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_void()
  return(ggevodyn)
}

animate_freq_dynamics <- function(freq_df, bw=0.05, bc="grey75", fps=10, ani_f_out = "evo_dynamics.gif", step_size=1, start_time=NULL, end_time=NULL, interpolation_steps = 0){
  ### FOR TESTING ###
  # freq_df <- pos_df
  # step_size <- 10
  # start_time=NULL
  # end_time=NULL
  # ani_f_out = "fish.mp4"
  # out_w <- 6
  # out_h <- 3
  ############
  final_time <- max(pos_df$x)
  if(is.null(start_time)){
    start_time <- 1
  }
  if(is.null(end_time)){
    end_time <- final_time
  }

  time_pts <- seq(start_time, end_time, by = step_size)
  if(!final_time %in% time_pts){
    time_pts <- c(time_pts, final_time)
  }

  if (file.exists(ani_f_out)) file.remove(ani_f_out)
  movie_format <- strsplit(ani_f_out, ".", fixed = T)[[1]]
  movie_format <- movie_format[length(movie_format)]
  if(movie_format == "gif"){
    img <- image_graph(600, 340, res = 1000)
    out <- lapply(time_pts, function(x){
      fish_plot_at_time <- plot_freq_dynamics(freq_df, end_time = x, bw=bw, bc=bc)
      print(fish_plot_at_time)
    })
    dev.off()
    animation <- image_animate(img, fps = fps, dispose = "previous")
    image_write(animation, ani_f_out)
    return(animation)
  }else{
    tmp_img_dir <- ".temp_fish_imgs/"
    dir.create(tmp_img_dir, showWarnings = F)
    n_time_pts <- length(time_pts)

    img_prefix <- paste0("%0",nchar(as.character(n_time_pts)), 'd.png')
    for(i in seq(n_time_pts)){
      end_time <- time_pts[i]

      fish_plot_at_time <- plot_freq_dynamics(freq_df, end_time = end_time, bw=bw, bc=bc)
      f_out <- paste0(tmp_img_dir, sprintf(img_prefix, i))
      ggsave(f_out, fish_plot_at_time, width = 6, height = 3)
    }
    img_dir_prefix <- paste0(tmp_img_dir, img_prefix)
    animate_string <- paste("ffmpeg -i", img_dir_prefix,
                            "-pix_fmt", "yuv420p",
                            # "-r", "2",
                            # "-crf", "15",
                            paste0("-vf fps=", fps),
                            "-b", "5000k",
                            "-vcodec", "libx264",
                            "-profile:v", "high",
                            "-preset", "veryslow",
                            ani_f_out)
    system(animate_string)
    unlink(tmp_img_dir, recursive = T)
    return(NULL)

  }

}

get_last_idx_above_0 <- function(x){
  idx <- tail(which(x > 0), 1)
  if(length(idx)==0){
    ###clone was never above 0
    idx <- 0
  }
  return(idx)
}

get_last_idx_at_0 <- function(x){
  idx <- tail(which(x == 0), 1)
  if(length(idx)==0){
    ### clone was never 0
    idx <- 0
  }
  return(idx)
}

scale_values <- function(x, out_range=c(0, 1)){
  a <- min(out_range)
  b <- max(out_range)
  in_min <- min(x, na.rm = T)
  in_max <- max(x, na.rm = T)

  scaled_x <- (b-a)*(x-in_min)/(in_max - in_min) + a

  return(scaled_x)
}

plot_freq_dynamics_alpha <- function(pos_df, n_time_pts=NULL, start_time=NULL, end_time=NULL, bw=0.05, bc="grey75"){
  ### FOR TESTING ###
  # n_time_pts <- 20
  # start_time=NULL
  # end_time=NULL
  # bw=0.05
  #####

  unique_time_pts <- unique(pos_df$x)
  if(is.null(n_time_pts)){
    n_time_pts <- length(unique_time_pts)
  }
  if(is.null(start_time)){
    start_time <- 0
  }
  if(is.null(end_time)){
    end_time <- max(pos_df$x)
  }

  view_time_pts <- unique_time_pts[unique_time_pts <= end_time]
  time_pt_idx <- unique(floor(seq(1, length(view_time_pts), length.out = n_time_pts)))

  time_pts <- view_time_pts[time_pt_idx]
  if(!end_time %in% time_pts){
    time_pts <- c(time_pts, end_time)
  }

  view_df <- pos_df[pos_df$x %in% time_pts, ]

  ggevodyn <- ggplot2::ggplot(view_df, ggplot2::aes(x=x, y=y, group=draw_order, fill=color, alpha=alpha)) +
    ggplot2::geom_polygon(size=bw, color=bc) +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_void()
  return(ggevodyn)
}

view_as_hierarchy <- function(size_df, clones, parents, data_type="size", threshold=0.01, show_extinct=F, min_node_size=1, max_node_size=5, attribute_vals=NULL, attribute_cmap='viridis', clone_cmap='rainbow_soft', defualt_color_val="black", scale_node_by_size=T){
  ### FOR TESTING ###

  # clone_history_dir <- "../AdinCaProject/chandler/AdinCa/adincar_model/adincar_model_java_code/get_smart_some_comp_time_lag/sweeps/0.5_0.7000000000000001_1.0/6/"
  # clone_history_file <- paste0(clone_history_dir, "clone_history.csv")
  #
  # clone_history_df <- read.csv(clone_history_file, check.names = F, header = T, stringsAsFactors = F)
  # time_idx <- which(!is.na(as.numeric(colnames(clone_history_df))))
  # time_pts <- as.numeric(colnames(clone_history_df)[time_idx])
  # size_df <- clone_history_df[time_idx]
  # clones <- clone_history_df$head
  # parents <- clone_history_df$tail


  # size_df <- clone_df
  # clones <- clone_list
  # parents <- parent_list
  # attribute_vals <- clone_antigenicity
  # threshold <- 0.01
  # node_sizes <- 5
  # clone_cmap='rainbow_soft'
  #####

  roots_root_idx <- which(!parents %in% clones)
  roots_root <- parents[roots_root_idx]
  root_id <- clones[roots_root_idx]



  if(data_type=="size"){
    ### data are sizes of each clone. Here, the frequency of each mutation will be determined
    freq_df <- get_mutation_freq_df(size_df, clones = clones, parents = parents)
    freq_mat <- as.matrix(freq_df)/max(freq_df)

  }else{
    check_freq_mat(size_df, clones, parents)
    freq_mat <- size_df
  }

  if(is.null(attribute_vals)){
    if(clone_cmap[1] %in% names(colormap::colormaps)){

      cmap_vals <- get_clone_color(length(clones), cmap = clone_cmap)
    }else{
      cmap_vals <- clone_cmap
    }
  }else{
    cmap_vals <- get_attribute_colors(attribute_vals, cmap=attribute_cmap)
  }

  filtered_info <- filter_freq_mat(clones, parents, freq_mat, cmap_vals, threshold)
  filtered_parents <- filtered_info$parents
  filtered_colors <- filtered_info$color_vals
  filtered_clones <- filtered_info$clones
  filtered_idx <- filtered_info$thresh_idx
  final_sizes <- freq_mat[filtered_idx, ncol(freq_mat)]

  if(show_extinct==F){
    last_idx_above_0 <- as.numeric(unlist(apply(filtered_info$freq_mat,1,get_last_idx_above_0)))
    numeric_times <- as.numeric(colnames(size_df))
    last_times_at_0 <- numeric_times[last_idx_above_0]
    not_extinct_idx <- which(as.numeric(last_times_at_0) == max(last_times_at_0))


    extant_parents <- filtered_parents[not_extinct_idx]
    extant_clones <- filtered_clones[not_extinct_idx]
    extant_colors <- filtered_colors[not_extinct_idx]
    extant_sizes <- final_sizes[not_extinct_idx]

    if(scale_node_by_size & length(unique(extant_sizes))>1){
      node_sizes <- scale_values(extant_sizes, out_range=c(min_node_size, max_node_size))
    }else{
      node_sizes <- max_node_size
    }

    vert_df <- data.frame("name"=extant_clones, "size"=node_sizes)
    if(!is.null(extant_colors)){
      vert_df$color <- extant_colors
    }else{
      vert_df$color <- "black"
    }

    edge_df <- data.frame("from"=extant_parents, "to"= extant_clones)
    edge_df <- subset(edge_df, to != root_id)

  }else{
    if(scale_node_by_size & length(unique(final_sizes))>1){
      node_sizes <- scale_values(final_sizes, out_range=c(min_node_size, max_node_size))
      node_sizes[final_sizes==0] <- 0 ###Makes clear that some populations are actually extinct
    }else{
      node_sizes <- max_node_size
    }

    vert_df <- data.frame("name"=filtered_clones, "size"=node_sizes)
    if(!is.null(filtered_colors)){
      vert_df$color <- filtered_colors
    }else{
      vert_df$color <- "black"
    }
    edge_df <- data.frame("from"=filtered_parents, "to"= filtered_clones)
    edge_df <- subset(edge_df, to != root_id)
  }

  clone_graph <- graph_from_data_frame(edge_df, vertices = vert_df)
  dendro <- ggraph(clone_graph, 'igraph', algorithm = 'tree', circular = F) +
    geom_edge_diagonal() +
    geom_node_point(aes(size=size, color=color)) +
    scale_color_identity() +
    ggforce::theme_no_axes()

  return(dendro)

}

#' Convert XYZ
#'
#' Builds...
#' @param edges_df Description here.
#' @param long_pop_sizes_df Description here.
#' @param time_col_name File to plot the animation to.
#' @param clone_col_name Description here.
#' @param parent_col_name Description here.
#' @param size_col_name Description here.
#' @keywords long_to_wide_size_df
#' @export
#' @examples
#' ## Using Provided Example Data
#' ## Data is in a long format consisting of two files:
#' ## 1. Population Data Frame
#' head(pop_df)
#' ## 2. Edges Data Frame
#' head(edges_df)
#'
#' ## For EvoFreq to use this data it must convert it.
#' wide_df <- long_to_wide_size_df(edges_df, pop_df, "Generation", "Identity", "Parent", "Population")
#'
#' ## Pull apart the necessary information for creating the frequency dynamics object
#' clone_list <- wide_df$clones
#' parent_list <- wide_df$parents
#' clone_df <- wide_df$wide_size_df
#' timepoints = sort(unique(pop_df$Generation))
#'
#' ## Get the frequency dynamics object using get_freq_dynamics()
#'
#' ## Plot the frequency dynamics using plot_freq_dynamics()
long_to_wide_size_df <- function(edges_df, long_pop_sizes_df, time_col_name, clone_col_name, parent_col_name, size_col_name){
  ### For testing ###
  # library(ggmuller)
  # long_pop_sizes_df <- example_pop_df
  # edges_df <- example_edges
  # time_col_name <- "Generation"
  # parent_col_name <- "Parent"
  # clone_col_name <- "Identity"
  # size_col_name <- "Population"
  #####
  long_pop_sizes_df <- as.data.frame(long_pop_sizes_df)
  wide_size_df <- reshape(as.data.frame(long_pop_sizes_df), timevar = time_col_name, idvar = clone_col_name, direction = "wide")

  parents <- edges_df[, parent_col_name]
  clones <- edges_df[, clone_col_name]
  clones_in_size_df <- unique(long_pop_sizes_df[,clone_col_name])
  clones_not_in_clone_list <- clones_in_size_df[!clones_in_size_df %in% clones]
  if(length(clones_not_in_clone_list)>0){
    if(length(clones_not_in_clone_list)>1){
      stop("More thant 1 clone in size data frame but not in clone list")
    }
    warning("Clone in size data frame but not in clone list. Assuming this is the root")
    rand_parent_id <- runif(1)
    clones <- c(clones, clones_not_in_clone_list)
    parents <- c(parents, rand_parent_id)
    edges_df <- data.frame(clones, parents)
    names(edges_df) <- c(clone_col_name, parent_col_name)
  }

  wide_size_df <- merge(wide_size_df, edges_df)

  ordered_clones <- wide_size_df[,clone_col_name]
  ordered_parents <- wide_size_df[,parent_col_name]
  ordered_size_df <- wide_size_df[,!colnames(wide_size_df) %in% c(clone_col_name, parent_col_name)]

  return(list("wide_size_df"=ordered_size_df, "parents"=ordered_parents, "clones"=ordered_clones))
}
