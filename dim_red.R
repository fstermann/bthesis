source("auxillary.R")

# Dimension Reduction
library(uwot)
library(Rtsne)

#===============================================================================
# INCLUDED FUNCTIONS
#
# -get_embedding
# -get_umap_embedding
# -get_tsne_embedding
# -get_pca_embedding
# 
# -show_all_dim_red
# -show_all_scenarios
# -show_all_correlations
# -show_umap_correlations
# -show_tsne_correlations
# -show_pca_correlations
# -show_umap_parameters
# 
# -get_param_grid
#
#===============================================================================


#===============================================================================
# .----------------------------.
# |     EMBEDDING FUNCTIONS    |
# '----------------------------'
#===============================================================================

get_embedding <- function(data, method, params, seed = 1234) {
  if (method == "UMAP") {
    emb <- get_umap_embedding(data = data, n_neighbors = params$n_neighbors, min_dist = params$min_dist, seed = seed)
  } else if (method == "t-SNE") {
    emb <- get_tsne_embedding(data = data, tsne_iter = params$tsne_iter, perplexity = params$perplexity, seed = seed)
  } else if (method == "PCA") {
    emb <- get_pca_embedding(data = data, center = params$pca_center, scale. = params$pca_scale, seed = seed)
  }
  return(emb)
}

get_umap_embedding <- function(data, n_neighbors = 15, min_dist = 0.01, seed = 1234) {
  set.seed(seed)
  umap_raw <- umap(data, n_neighbors = n_neighbors, min_dist = min_dist)
  embedding <- data.frame(x = umap_raw[,1], 
                          y = umap_raw[,2])
  return(embedding)
}

get_tsne_embedding <- function(data, tsne_iter = 1000, perplexity = 30, seed = 1234) {
  set.seed(seed)

  tsne_raw <- Rtsne(data, perplexity = 30, max_iter = tsne_iter, pca = FALSE, check_duplicates = FALSE)
  embedding <- data.frame(x = tsne_raw$Y[,1], 
                          y = tsne_raw$Y[,2])
  return(embedding)
}

get_pca_embedding <- function(data, center = TRUE, scale. = TRUE, seed = 1234) {
  set.seed(seed)
  pca_raw <- prcomp(data, center = center, scale. = scale.)
  embedding <- data.frame(x = pca_raw$x[,1], 
                          y = pca_raw$x[,2])
  return(embedding)
}


#===============================================================================
# .-----------------------.
# |     PLOT FUNCTIONS    |
# '-----------------------'
#===============================================================================
  
# Compare UMAP, t-SNE, PCA
show_all_dim_red <- function(data, data_obj, with_table = TRUE, group = NULL) {
  gp <- data_obj$GenerationParameter
  vp <- data_obj$VariableParameters
  
  umap_raw <- get_umap_embedding(data)
  tsne_raw <- get_tsne_embedding(data)
  pca_raw <- get_pca_embedding(data)
  if (is.null(group)) {group <- factor(c(rep(1, sum(gp[[3]][["ns"]]))))}
  
  projection <- data.frame(x_umap = umap_raw$x, 
                           y_umap = umap_raw$y,
                           x_tsne = tsne_raw$x, 
                           y_tsne = tsne_raw$y,
                           PC1 = pca_raw$x,
                           PC2 = pca_raw$y,
                           Group = group)
  
  p_umap <- ggplot(projection, aes(x = x_umap, y = y_umap, col = Group, alpha = 0.5)) + geom_point() + theme(legend.position = "none") + scale_color_manual(values = TRUE_COLORS)
  p_tsne <- ggplot(projection, aes(x = x_tsne, y = y_tsne, col = Group, alpha = 0.5)) + geom_point() + theme(legend.position = "none") + scale_color_manual(values = TRUE_COLORS)
  p_pca <- ggplot(projection, aes(x = PC1, y = PC2, col = Group, alpha = 0.5)) + geom_point() + scale_color_manual(values = TRUE_COLORS)

  if (with_table) {
    grid_params <- get_param_grid(data_obj)
    
    p_all <- grid.arrange(grid_params, p_umap, p_tsne, p_pca, layout_matrix = cbind(1, 2, 3, 4))
  } else {
    p_all <- grid.arrange(p_umap, p_tsne, p_pca, ncol = 3)
  }

  return(p_all)
}

# --> Compare all given scenarios
show_all_scenarios <- function(read_params, scenarios = c("01", "02", "03"), group = NULL) {
  plots <- list()
  
  for (sc in scenarios) {
    read_params$Scenario <- sc
    data_files <- get_data(read_params)
    data <- data_files$data
    data_obj <- data_files$data_obj
    
    p_sc <- show_all_dim_red(data, data_obj, with_table = FALSE, group = group)
    plots[[sc]] <- grid.arrange(p_sc, left = paste0("Scenario ", sc))
  }
  
  p_all <- grid.arrange(grobs = plots, nrow = length(scenarios))
  
  return(p_all)
}

# --> Compare all correlations for single scenario
show_all_correlations <- function(read_params, correlations, tsne_iter = 50, group = NULL) {
  plots <- list()
  
  for (i in 1:length(correlations)) {
    read_params$MinRho <- correlations[[i]]$MinRho
    read_params$MaxRho <- correlations[[i]]$MaxRho
    read_params$Method <- correlations[[i]]$Method
    
    data_files <- get_data(read_params)
    data <- data_files$data
    data_obj <- data_files$data_obj
    
    gp <- data_obj$GenerationParameter
    
    if (is.null(group)) group <- factor(c(rep(1, sum(gp[[3]][["ns"]]))))
    if (i == 1) grid_params <- get_param_grid(data_obj)
    
    p_r <- show_all_dim_red(data, data_obj, with_table = FALSE, group = group)
    plots[[toString(i)]] <- grid.arrange(p_r, left = paste0("Correlation ", gp[[2]]$MinRho, "-", gp[[2]]$MaxRho))
  }
  
  grid_plots <- grid.arrange(grobs = plots, ncol = 1)
  
  p_all <- grid.arrange(grid_params, grid_plots, layout_matrix = rbind(c(1, 2, 2, 2)))
  return(p_all)
}

# Compare all correlations for single scenario, only UMAP
show_umap_correlations <- function(read_params, correlations, n_neighbors = 15, min_dist = 0.01, group = NULL) {
  plots <- list()
  
  for (i in 1:length(correlations)) {
    read_params$MinRho <- correlations[[i]]$MinRho
    read_params$MaxRho <- correlations[[i]]$MaxRho
    read_params$Method <- correlations[[i]]$Method
    
    data_files <- get_data(read_params)
    data <- data_files$data
    data_obj <- data_files$data_obj
    
    gp <- data_obj$GenerationParameter
    
    if (is.null(group)) group <- factor(c(rep(1, sum(gp[[3]][["ns"]]))))
    if (i == 1) grid_params <- get_param_grid(data_obj)
    
    umap_raw <- get_umap_embedding(data, n_neighbors = n_neighbors, min_dist = min_dist, seed = read_params$Seed)
    projection <- data.frame(x_umap = umap_raw[,1], 
                             y_umap = umap_raw[,2],
                             Group = group)
    p_umap <- ggplot(projection, aes(x = x_umap, y = y_umap, col = Group, alpha = 0.5)) + 
      geom_point() + theme(legend.position = "none") + scale_color_manual(values = TRUE_COLORS) +
      annotate("text", x = Inf, y = Inf, label = paste0("Rho: [",gp[[2]]$MinRho,", ",gp[[2]]$MaxRho,"]"), hjust = "inward", vjust = "inward")
    
    plots[[toString(i)]] <- p_umap
  }
  grid_plots <- grid.arrange(grobs = plots)
  
  p_all <- grid.arrange(grid_params, grid_plots, layout_matrix = rbind(c(1, 2, 2, 2)))
  
  return(p_all)
}

# Compare all correlations for single scenario, only TSNE
show_tsne_correlations <- function(read_params, correlations, tsne_iter = 50, perplexity = 30, group = NULL) {
  plots <- list()
  
  for (i in 1:length(correlations)) {
    read_params$MinRho <- correlations[[i]]$MinRho
    read_params$MaxRho <- correlations[[i]]$MaxRho
    read_params$Method <- correlations[[i]]$Method
    
    data_files <- get_data(read_params)
    data <- data_files$data
    data_obj <- data_files$data_obj
    
    gp <- data_obj$GenerationParameter
    
    if (is.null(group)) group <- factor(c(rep(1, sum(gp[[3]][["ns"]]))))
    if (i == 1) grid_params <- get_param_grid(data_obj)
    
    tsne_raw <- get_tsne_embedding(data, tsne_iter = tsne_iter, perplexity = perplexity, seed = read_params$Seed)
    projection <- data.frame(x_tsne = tsne_raw$x, 
                             y_tsne = tsne_raw$y,
                             Group = group)
    p_tsne <- ggplot(projection, aes(x = x_tsne, y = y_tsne, col = Group, alpha = 0.5)) + 
      geom_point() + theme(legend.position = "none") + scale_color_manual(values = TRUE_COLORS) +
      annotate("text", x = Inf, y = Inf, label = paste0("Rho: [",gp[[2]]$MinRho,", ",gp[[2]]$MaxRho,"]"), hjust = "inward", vjust = "inward")
    
    plots[[toString(i)]] <- p_tsne
  }
  grid_plots <- grid.arrange(grobs = plots)
  
  p_all <- grid.arrange(grid_params, grid_plots, layout_matrix = rbind(c(1, 2, 2, 2)))
  
  return(p_all)
}

# Compare all correlations for single scenario, only PCA
show_pca_correlations <- function(read_params, correlations, center = TRUE, scale. = TRUE,  group = NULL) {
  plots <- list()
  
  for (i in 1:length(correlations)) {
    read_params$MinRho <- correlations[[i]]$MinRho
    read_params$MaxRho <- correlations[[i]]$MaxRho
    read_params$Method <- correlations[[i]]$Method
    
    data_files <- get_data(read_params)
    data <- data_files$data
    data_obj <- data_files$data_obj
    
    gp <- data_obj$GenerationParameter

    if (is.null(group)) {group <- factor(c(rep(1, sum(gp[[3]][["ns"]]))))}
    if (i == 1) grid_params <- get_param_grid(data_obj)
    
    pca_raw <- get_pca_embedding(data, center = center, scale. = scale., seed = read_params$Seed)
    projection <- data.frame(PC1 = pca_raw$x, 
                             PC2 = pca_raw$y,
                             Group = group)
    p_umap <- ggplot(projection, aes(x = PC1, y = PC2, col = Group, alpha = 0.5)) + 
      geom_point() + theme(legend.position = "none") + scale_color_manual(values = TRUE_COLORS) +
      annotate("text", x = Inf, y = Inf, label = paste0("Rho: [",gp[[2]]$MinRho,", ",gp[[2]]$MaxRho,"]"), hjust = "inward", vjust = "inward")
    
    plots[[toString(i)]] <- p_umap
  }
  grid_plots <- grid.arrange(grobs = plots)
  
  p_all <- grid.arrange(grid_params, grid_plots, layout_matrix = rbind(c(1, 2, 2, 2)))
  
  return(p_all)
}

# Compare all UMAP parameters for single dataset
show_umap_parameters <- function(read_params, n_neighbors_range = c(2, 4, 8, 16, 32, 64), min_dist_range = c(0, 0.25, 0.5, 0.75, 1), group = NULL, seed = 1234) {
  plots <- list()
  plots_mindist <- list()
  
  data_files <- get_data(read_params)
  data <- data_files$data
  data_obj <- data_files$data_obj
  gp <- data_obj$GenerationParameter
  if (is.null(group)) {group <- factor(c(rep(1, sum(gp[[3]][["ns"]]))))}
  
  for (i in 1:length(n_neighbors_range)) {
    n_neighbors <- n_neighbors_range[i]
    
    for (j in 1:length(min_dist_range)) {
      min_dist <- min_dist_range[j]
      
      if (i == 1 & j == 1) {grid_params <- get_param_grid(data_obj)}
      
      set.seed(seed)
      umap_raw <- umap(data, n_neighbors = n_neighbors, min_dist = min_dist)
      projection <- data.frame(x_umap = umap_raw[,1], 
                               y_umap = umap_raw[,2],
                               Group = group)
      p_umap <- ggplot(projection, aes(x = x_umap, y = y_umap, col = Group, alpha = 0.5)) + 
        geom_point() + theme(legend.position = "none") + scale_color_manual(values = TRUE_COLORS) + 
        theme(axis.title.y = element_blank(), axis.title.x = element_blank())
      
      if (i == 1) {p_umap <- p_umap + ggtitle(paste0("min_dist: ",min_dist)) + theme(plot.title = element_text(hjust = 0.5))} 
      
      plots[[paste0(i,"_",j)]] <- p_umap 
    }
    
    plots_mindist[[toString(i)]] <- grid.arrange(grobs = plots, ncol = length(min_dist_range), left = paste0("n_neighbors: ",n_neighbors))
    plots <- list()
  }
  grid_plots <- grid.arrange(grobs = plots_mindist, nrow = length(n_neighbors_range))
  
  p_all <- grid.arrange(grid_params, grid_plots, layout_matrix = rbind(c(1, 2, 2, 2)))
  
  return(p_all)
}


#===============================================================================
# .------------------------.
# |     OTHER FUNCTIONS    |
# '------------------------'
#===============================================================================

get_param_grid <- function(data_obj, max_vars = 10, font_size = 8) {
  gp <- data_obj$GenerationParameter
  vp <- data_obj$VariableParameters
  
  # Scenario Parameter
  df_scenario <- data.frame(1:length(gp[[3]][["ns"]]), gp[[3]][["ns"]], 
                            paste0("[", gp[[1]][["MinMeans"]], ", ", gp[[1]][["MaxMeans"]], "]"),
                            paste0("[", gp[[1]][["MinDispersions"]], ", ", gp[[1]][["MaxDispersions"]], "]"),
                            paste0("[", gp[[1]][["MinZProbs"]], ", ", gp[[1]][["MaxZProbs"]], "]"))
  colnames(df_scenario) <-  c("Group", "n", "Means", "Dispersions", "ZProbs")
  
  p_scenario_params <- tableGrob(df_scenario, rows = NULL, theme = ttheme_default(base_size = font_size))
  
  # Variable Parameter
  n_equal <- gp[[3]][["EqualVars"]]
  n_diff <- gp[[3]][["DiffVars"]]
  
  if (n_diff > max_vars) {
    n_split <- round(max_vars/2, 0)
    df_variable <- rbind(rep("...",times = ncol(vp)))
    colnames(df_variable) <- colnames(vp)
    df_variable <- rbind(vp[(n_equal+1):(n_equal+n_split),], df_variable, vp[(n_equal+n_diff-n_split+1):(n_equal+n_diff),])
    
    p_variable_params <- tableGrob(df_variable, rows = NULL, theme = ttheme_default(base_size = font_size))
  } else {
    p_variable_params <- tableGrob(vp[(n_equal+1):(n_equal+n_diff),], rows = NULL, theme = ttheme_default(base_size = font_size))
  }
  
  # Combine Grids
  grid_params <- grid.arrange(p_scenario_params, p_variable_params, 
                              top = textGrob(paste0("Scenario ", gp[[1]][["Scenario"]], 
                                                    "\np = [",n_equal,"/",n_diff, "]"), 
                                             gp=gpar(fontsize=font_size+4)), 
                              layout_matrix = cbind(c(1, 2, 2, 2)))
  
  return(grid_params)
}
