source("auxillary.R")
source("dim_red.R")

library(mclust)
library(data.table)
library(dplyr)
library(tidyr)

#===============================================================================
# INCLUDED FUNCTIONS
#
# -get_kmeans
# -get_hclust
# -get_mclust
# -get_clustering_results
# -is_exception
#
# -get_ARI
# -get_ARI_ForSeeds
# --get_ARI_ForParams
# -get_ARI_boxplots_byP
# -get_ARI_boxplots_byDR
#
# -show_CL_correlations
# --show_CL_embedding
#
# -get_bottom_grid
# -get_x_titles_grid
# -get_y_titles_grid
#
# -get_ARI_plot_DF
# -transform_df
#===============================================================================


#===============================================================================
# .----------------------------.
# |     EMBEDDING FUNCTIONS    |
# '----------------------------'
#===============================================================================

get_kmeans <- function(data, centers = 2, nstart = 50, seed = 1234) {
  set.seed(seed)
  ret <- kmeans(x = data, centers = centers, nstart = nstart)
  return(ret$cluster)
}

get_hclust <- function(data, k = 2, method = "ward.D2", dmeth = "euclidean", seed = 1234) {
  set.seed(seed)
  ret <- hclust(dist(data, method = dmeth), method = method)
  ret <- cutree(ret, k)
  return(ret)
}

get_mclust <- function(data, G = 2, seed = 1234, restrict = FALSE) {
  set.seed(seed)
  # Exlcude VVE model if restrict = TRUE
  m <- mclust.options("emModelNames")
  if (restrict) m <- m[!m=="VVE"]
  invisible(capture.output(ret <- Mclust(data, G = 2, modelNames = m)))

  ret <- list("Cluster" = ret$classification, "Model" = ret$modelName)

  return(ret)
}

# Apply Clustering method to embedding or data
get_clustering_results <- function(embedding, data, group, seed = 1234, restrict_mclust = FALSE) {
  kmeans_dim_red <- get_kmeans(embedding, seed = seed)
  hclust_dim_red <- get_hclust(embedding, seed = seed)
  mclust_dim_red <- get_mclust(embedding, seed = seed, restrict = restrict_mclust)
  
  df <- data.frame(EmbX = embedding$x, EmbY = embedding$y, 
                   TrueGroup = group,
                   kmeansGroup = factor(kmeans_dim_red),
                   hclustGroup = factor(hclust_dim_red),
                   mclustGroup = factor(mclust_dim_red$Cluster),
                   mclustModelN = mclust_dim_red$Model)
  return(df)
}

is_exception <- function(sc, min_rho, p, dimred, seed) {
  # These datasets/embedding caused the mclust algorithm to crash.
  # I contacted the maintainer Luca Scrucca, this apparently is a known 
  # problem. The VVE model with 2 Clusters can take a very long time
  # to compute in some special cases, which is why i chose to exclude the
  # model in these cases.
  
  is_exception <- all(sc == "01", min_rho == 0.7, p == 2000, dimred == "UMAP", seed == 30092018) |
                  all(sc == "03", min_rho == 0.7, p == 250, dimred == "UMAP", seed == 30092019) |
                  all(sc == "03", min_rho == 0.7, p == 1000, dimred == "UMAP", seed == 30092064)
}

#===============================================================================
# .---------------------.
# |         ARI         |
# '---------------------'
#===============================================================================

# Get the Adjusted Rand Index for each Clustering Method
get_ARI <- function(clustering_results) {
  cr <- clustering_results
  results <- c(adjustedRandIndex(cr$TrueGroup, cr$kmeansGroup),
               adjustedRandIndex(cr$TrueGroup, cr$hclustGroup),
               adjustedRandIndex(cr$TrueGroup, cr$mclustGroup))
  
  results <- round(results, 2)
  names(results) <- c("kmeans", 
                      "hclust", 
                      "mclust")
  
  df <- data.frame(rbind(results))
  
  return(df)
}

# Get ARI for different Seeds
get_ARI_ForSeeds <- function(seeds, read_params, dimred_params, correlations, ps, 
                             dimreds, clusterings, iterations = 10) {
  dfs <- list()
  for (seed in seeds) {
    read_params$Seed <- seed
    for (i in 0:(iterations-1)) {
      dfs[[paste0(toString(seed), "+", i)]] <- get_ARI_ForParams(read_params = read_params, dimred_params = dimred_params, correlations = correlations, ps = ps, 
                                                                 dimreds = dimreds, clusterings = clusterings, seed = seed+i)
    }
  }
  
  dfs[["Average"]] <- data.frame(rbindlist(dfs)[,lapply(.SD, function(x) round(mean(x), 2)), list(Rho, p)])
  dfs[["Min"]] <- data.frame(rbindlist(dfs[1:length(seeds)])[,lapply(.SD, min), list(Rho, p)])
  dfs[["Max"]] <- data.frame(rbindlist(dfs[1:length(seeds)])[,lapply(.SD, max), list(Rho, p)])
  
  return(dfs)
}

# Get ARI for all DimRed, Clustering, P's, Correlations
get_ARI_ForParams <- function(read_params, dimred_params, correlations, ps, dimreds, clusterings, seed = 1234) {
  methods <- paste(rep(dimreds, each = length(clusterings)), clusterings, sep = "+\n")
  df <- data.frame(matrix(vector(), nrow = 0, ncol = 2 + length(methods),
                   dimnames = list(c(), c("Rho", "p", methods))))
  
  for (correlation in correlations) {
    read_params$MinRho <- correlation$MinRho
    read_params$MaxRho <- correlation$MaxRho
    read_params$Method <- correlation$Method
    
    for (p in all_ps) {
      read_params$DiffVars <- p
      data_files <- get_data(read_params)
      data <- data_files$data
      
      vals <- c(paste0("[", sprintf("%.1f", correlation$MinRho), "-", sprintf("%.1f", correlation$MaxRho), "]"), p)

      for (dimred in dimreds) {
        cat("Calculating ARI for Scenario: ", read_params$Scenario, 
            ", Rho: ", paste0("[", sprintf("%.1f", correlation$MinRho), "-", sprintf("%.1f", correlation$MaxRho), "]"), 
            ", p: ", p, ", DimRed: ", dimred, ", Seed: ", seed, "          \r", sep = "")
        
        embedding <- get_embedding(data =  data, params = dimred_params, method = dimred, seed = seed)
        if (dimred_params$scale_embedding & dimred != "UMAP") {embedding <- data.frame(scale(embedding))}
        
        group <- factor(c(rep(1, times = read_params$ns[1]), rep(2, times = read_params$ns[2])))

        cr <- get_clustering_results(embedding, data, group, seed, restrict_mclust = is_exception(read_params$Scenario, 
                                                                                                  read_params$MinRho,
                                                                                                  p, dimred, seed))
        ari <- get_ARI(cr)
        
        vals <- c(vals, ari[,clusterings])
      }
      df[nrow(df)+1,] <- vals
    }
  }
  
  return(df)
}

# Get ARI Boxplots for each P
get_ARI_boxplots_byP <- function(ari_df_sc1, ari_df_sc2, ari_df_sc3) {
  df <- get_ARI_plot_DF(ari_df_sc1, ari_df_sc2, ari_df_sc3)
  plots <- list()
  
  for (p in unique(df1$p)) {
    plt <- ggplot(df[df$p == p,], aes(x = Clustering, y = ARI, fill = Rho)) + 
      geom_boxplot() +
      scale_fill_manual(values = RHO_COLORS) + facet_grid(rows = vars(Scenario), cols = vars(DimRed)) +
      scale_y_continuous(breaks = c(0.25, 0.50, 0.75, 1.00)) +
      xlab("") +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size = FONT_SIZE_PLOT),
            strip.text = element_text(face = "bold")) 

    grid_text <- textGrob(paste0("p = ", p), gp = gpar(fontsize = FONT_SIZE_TEXT))
    grid_legend <- get_shared_boxplot_legend(unique(df$Rho), RHO_COLORS, FONT_SIZE_TEXT)
    
    plots[[toString(p)]] <- grid.arrange(grobs = list(plt, grid_text, grid_legend), 
                                         layout_matrix = rbind(c(1, 1, 1), c(1, 1, 1), 
                                                               c(1, 1, 1), c(1, 1, 1),
                                                               c(1, 1, 1), c(1, 1, 1),
                                                               c(2, 3, 3)))
  }
  
  return(plots)
}
# Get ARI Boxplots for each DimRed
get_ARI_boxplots_byDR <- function(ari_df_sc1, ari_df_sc2, ari_df_sc3) {
  df <- get_ARI_plot_DF(ari_df_sc1, ari_df_sc2, ari_df_sc3)
  plots <- list()
  
  for (dimred in c("UMAP", "t-SNE", "PCA")) {
    plt <- ggplot(df[df$DimRed == dimred,], aes(x = Rho, y = ARI, fill = p)) + 
      geom_boxplot() +
      scale_fill_manual(values = P_COLORS) + 
      facet_grid(rows = vars(Scenario), cols = vars(Clustering)) +
      scale_y_continuous(breaks = c(0.25, 0.50, 0.75, 1.00)) +
      xlab("") +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size = FONT_SIZE_PLOT),
            strip.text = element_text(face = "bold")) 
    
    grid_text <- textGrob(dimred, gp = gpar(fontsize = FONT_SIZE_TEXT))
    grid_legend <- get_shared_boxplot_legend_p(names(P_COLORS), P_COLORS, FONT_SIZE_TEXT)
    
    plots[[toString(dimred)]] <- grid.arrange(grobs = list(plt, grid_text, grid_legend), 
                                         layout_matrix = rbind(c(1, 1, 1), c(1, 1, 1), 
                                                               c(1, 1, 1), c(1, 1, 1),
                                                               c(1, 1, 1), c(1, 1, 1),
                                                               c(2, 3, 3)))
  }
  
  return(plots)
}

#===============================================================================
# .----------------------------.
# |         CLUSTERING         |
# '----------------------------'
#===============================================================================

# Show all embedding for each correlation with clustering results
show_CL_correlations <- function(dim_red_method, dimred_params, read_params, correlations, group = NULL, seed = 1234) {
  plots <- list()
  df_ari <- list()
  df_ari_raw <- list()
  
  txt_methods <- paste0(dim_red_method, " + ", c("Wahre Cluster", "kmeans", "hclust", "mclust"))
  txt_correlations <- paste0("[", sprintf("%.1f", sapply(all_correlations, function(x) x$MinRho)), "-", 
                             sprintf("%.1f", sapply(all_correlations, function(x) x$MaxRho)), "]")
  
  for (i in 1:length(correlations)) {
    read_params$MinRho <- correlations[[i]]$MinRho
    read_params$MaxRho <- correlations[[i]]$MaxRho
    read_params$Method <- correlations[[i]]$Method
    
    data_files <- get_data(read_params)
    data <- data_files$data
    data_obj <- data_files$data_obj
    
    gp <- data_obj$GenerationParameter
    vp <- data_obj$VariableParameter
    
    # Get Embedding
    embedding <- get_embedding(data, method = dim_red_method, params = dimred_params, seed = seed)
    if (dimred_params$scale_embedding) {embedding <- data.frame(scale(embedding))}
    
    # Get Clustering Results
    cl_results <- get_clustering_results(embedding, data, group, seed = seed)
    ari <- get_ARI(cl_results)
    df_ari[[i]] <- cbind(data.frame("Rho" = txt_correlations[i]), ari[,1:3])
    
    plots[[i]] <- show_CL_embedding(results = cl_results, dim_red_method = dim_red_method)
  }
  
  df_ari <- do.call("rbind", df_ari)
  
  grid_plots <- grid.arrange(grobs = plots, ncol = 1)
  grid_bottom <- get_bottom_grid(data_obj, df_ari, group, font_size = FONT_SIZE_TEXT)
  grid_x_titles <- get_x_titles_grid(txt_methods, font_size = FONT_SIZE_PLOT)
  grid_y_titles <- get_y_titles_grid(txt_correlations, font_size = FONT_SIZE_PLOT)
  
  p_all <- grid.arrange(grid_x_titles, grid_plots, grid_y_titles, grid_bottom, 
                        layout_matrix = rbind(c(rep(1, 14), NA),
                                              c(rep(2, 14), 3), c(rep(2, 14), 3), c(rep(2, 14), 3),
                                              c(rep(2, 14), 3), c(rep(2, 14), 3), c(rep(2, 14), 3),
                                              c(rep(2, 14), 3), c(rep(2, 14), 3), c(rep(2, 14), 3),
                                              rep(4, 15), rep(4, 15)))
  return(p_all) 
}

# Show a single embedding with clustering results
show_CL_embedding <- function(results = NULL, embedding = NULL, data = NULL, group = NULL, dim_red_method = "") {
  if (is.null(results)) {
    df <- get_clustering_results(embedding, data, group, seed = seed)
  } else {
    df <- results
  }
  methods <- paste0(dim_red_method, " + ", c("Wahre Cluster", "kmeans", "hclust", "mclust"))
  
  df_comb <- data.frame("EmbX" = rep(df$EmbX, 4),
                        "EmbY" = rep(df$EmbY, 4),
                        "Group" = factor(c(df$TrueGroup, df$kmeansGroup, df$hclustGroup, df$mclustGroup)), 
                        "Method" = factor(rep(methods, each = 500), levels = methods))
  df_text <- data.frame(EmbX = Inf, EmbY = Inf,
                        Method = factor(paste0(dim_red_method, " + mclust")),
                        Group = factor(1))
  
  p_all <- ggplot(df_comb, aes(x = EmbX, y = EmbY, col = Group)) + 
    geom_point(size = 2, aes(alpha = 0.8)) + THEME_BLANK + scale_color_manual(values = TRUE_COLORS) + facet_wrap(~ Method, nrow = 1) +
    geom_text(data = df_text, label = df$mclustModelN[1], col = "black", vjust = 1.1, hjust = 1.1, size = 10) +
    theme(panel.spacing.x = unit(c(2, 0.5, 0.5), "cm")) + THEME_STRIP_BLANK
  
  return(p_all)
}

#===============================================================================
# .----------------------------------.
# |         HELPER FUNCTIONS         |
# '----------------------------------'
#===============================================================================

get_bottom_grid <- function(data_obj, df_ari, group, font_size = 8) {
  gp <- data_obj$GenerationParameter
  vp <- data_obj$VariableParameters
  p <- gp[[3]][["DiffVars"]]
  
  # Label
  grid_label <- textGrob(paste0("Szenario ", gp[[1]][["Scenario"]], 
                                "\np = ", p,
                                "\nseed = ", data_obj$Seed), 
                         gp = gpar(fontsize = FONT_SIZE_TEXT))
  
  # ARI
  title_ari <- textGrob("ARI", gp=gpar(fontsize = FONT_SIZE_TEXT), rot = 90)
  table_ari <- tableGrob(df_ari, rows = NULL, theme = ttheme_default(base_size = FONT_SIZE_TABLE, padding=unit(c(5, 7), "mm")))
  
  padding <- unit(1.5, "cm")
  table_ari <- gtable_add_cols(table_ari, 
                           widths = grobWidth(title_ari) + padding,
                           pos = 0)
  table_ari <- gtable_add_grob(table_ari, title_ari,
                           t=1, l=1, 
                           b = nrow(table_ari))
  
  # Legend
  grid_legend <- get_shared_point_legend(group, TRUE_COLORS, FONT_SIZE_TEXT)
  
  # All
  grid_bottom <- grid.arrange(grid_label, table_ari, grid_legend, 
                              layout_matrix = rbind(c(1, 2, 3)))
  
  return(grid_bottom)
}

get_x_titles_grid <- function(titles, font_size = 8) {
  grobs <- lapply(titles, textGrob, gp=gpar(fontsize = font_size, fontface = "bold")) #, vjust = 2)
  grid <- grid.arrange(grobs = grobs, nrow = 1)
  return(grid)
}

get_y_titles_grid <- function(titles, font_size = 8) {
  grobs <- lapply(titles, textGrob, gp=gpar(fontsize = font_size, fontface = "bold"), rot = -90) #, vjust = 2)
  grid <- grid.arrange(grobs = grobs, ncol = 1)
  return(grid)
}

# Helper to get ARI Data for plots
get_ARI_plot_DF <- function(ari_df_sc1, ari_df_sc2, ari_df_sc3) {
  df1 <- do.call(rbind, lapply(ari_df_sc1, transform_df)); df1$Scenario <- "01"
  df2 <- do.call(rbind, lapply(ari_df_sc2, transform_df)); df2$Scenario <- "02"
  df3 <- do.call(rbind, lapply(ari_df_sc3, transform_df)); df3$Scenario <- "03"
  df <- do.call(rbind, list(df1, df2, df3))
  
  return(df)
}

transform_df <- function(ari_df) {
  df <- data.frame(Rho = rep(ari_df$Rho, ncol(ari_df)-2), p = rep(ari_df$p, ncol(ari_df)-2),
                   ARI = unname(unlist(ari_df[,3:11])), Method = rep(colnames(ari_df[,3:11]), each = nrow(ari_df)))
  
  df <- separate(data = df, col = Method, into = c("DimRed", "Clustering"), sep = "\\.\\.")
  df$p <- factor(df$p, levels = unique(df$p))
  df$DimRed[df$DimRed == "t.SNE"] <- "t-SNE"
  df$DimRed <- factor(df$DimRed, levels = c("UMAP", "t-SNE", "PCA"))
  df$Clustering <- factor(df$Clustering, levels = c("kmeans", "hclust", "mclust"))
  
  df <- df %>% group_by(p, Rho) %>% mutate(Max = ifelse(ARI == max(ARI), 1, 0)) %>% ungroup()
  
  return(df)
}
