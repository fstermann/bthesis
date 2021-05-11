source("auxillary.R")

library(hexbin)

#===============================================================================
# INCLUDED FUNCTIONS
#
# -get_rho_stats
# --get_RMSR
# --get_RIR
# -get_rho_scatter_plot
# -get_rho_scatter_plots
# -get_rho_scatter_group_plots
# -get_rho_density_boxplots
#
# -get_rho_bounds_plot
# --read_rho_bounds
# 
# -get_variance_boxplots
#
#===============================================================================


#===============================================================================
# .---------------------------.
# |     RHO STAT FUNCTIONS    |
# '---------------------------'
#===============================================================================


get_rho_stats <- function(data, data_obj, group = 0) {
  corAll <- cor(data)
  corAll <- corAll[upper.tri(corAll)]
  if (group > 0) {
    tgtRho <- data_obj$Rhos[[group]]
  } else {
    tgtRho <- data_obj$Rho 
  }
  corTGT <- tgtRho[upper.tri(tgtRho)]
  
  min_rho <- data_obj$GenerationParameter[[2]]$MinRho
  max_rho <- data_obj$GenerationParameter[[2]]$MaxRho
  
  stats <- rbind(round(summary(corAll), 3), 
                 round(summary(corTGT), 3))
  stats <- data.frame("Type" = c("Simulated", "Target"),
                      stats,
                      "RMSR" = c(get_RMSR(tgtRho, cor(data), 3), 0),
                      "RhosInRange" = c(get_RIR(cor(data), min_rho, max_rho),
                                        get_RIR(tgtRho, min_rho, max_rho)))
  return(stats)
}

# Root Mean Squared Residuals
get_RMSR <- function(rho_target, rho, digits = 5) {
  rho_residuals <- rho - rho_target
  RMSR <- sqrt(sum(rho_residuals[lower.tri(rho_residuals)] * rho_residuals[lower.tri(rho_residuals)]) / sum(1:(ncol(rho)-1)))
  RMSR <- round(RMSR, digits)
  return(RMSR)
}

# Percentage of pairwise correlations between min_rho and max_rho
get_RIR <- function(rho, min_rho, max_rho) {
  rho <- rho[lower.tri(rho)]
  rir <- 1-(sum(rho < min_rho)+sum(rho > max_rho))/length(rho)
  return(paste0(round(rir*100, 2), "%"))
}

get_rho_scatter_plot <- function(data, data_obj, only_scatter = FALSE, group = NULL) {
  corAll <- cor(data)
  corAll <- corAll[upper.tri(corAll)]
  
  if (is.null(group)) {
    corTGT <- data_obj$Rho[upper.tri(data_obj$Rho)]
  } else {
    corTGT <- data_obj$Rhos[[group]][upper.tri(data_obj$Rhos[[group]])]
  }
  
  df <- data.frame("TargetCorrelation" = corTGT,
                   "SimulatedCorrelation" = corAll)
  xymin <- 0 #min(c(df$TargetCorrelation, df$SimulatedCorrelation))
  xymax <- 1 #max(c(df$TargetCorrelation, df$SimulatedCorrelation))
  
  p_rho <- ggplot(df) + stat_binhex(aes(x=TargetCorrelation, y=SimulatedCorrelation, alpha=..count..), bins = 100, fill="blue") +
    xlim(xymin, xymax) + ylim(xymin, xymax) + geom_abline(intercept = 0, slope = 1, col = "red") + theme(aspect.ratio=1)
  
  if (only_scatter) return(p_rho)

  df_dens <- data.frame("Rho" = c(df$TargetCorrelation, df$SimulatedCorrelation),
                        "Group" = rep(c("Target", "Simulated"), each = length(df$TargetCorrelation)))
  p_rho_dens <- ggplot(df_dens, aes(x = Rho, col = Group, fill = Group)) + geom_density(size = 1, alpha = 0.3)
  
  p_all <- grid.arrange(p_rho, p_rho_dens, ncol = 2)
  return(p_all)
}

get_rho_scatter_plots <- function(read_params, correlations) {
  plots <- list()
  
  first = TRUE
  for (r in correlations) {
    read_params$MinRho <- r$MinRho
    read_params$MaxRho <- r$MaxRho
    read_params$Method <- r$Method
    data_files <- get_data(read_params)
    data <- data_files$data
    data_obj <- data_files$data_obj
    
    p_rho <- get_rho_scatter_plot(data, data_obj, TRUE)
    p_rho <- p_rho + ggtitle(paste0("Correlation ", r[1], "-", r[2]))
    
    stats_df <- get_rho_stats(data, data_obj)
    p_stats1 <- tableGrob(stats_df[,1:7], rows = NULL, theme = ttheme_default(base_size = 8))
    p_stats2 <- tableGrob(stats_df[,c(1,8:9)], rows = NULL, theme = ttheme_default(base_size = 8))
    
    plots[[toString(r)]] <- grid.arrange(p_rho, p_stats1, p_stats2, layout_matrix = cbind(c(1, 1, 1, 2, 3)))
  }
  
  p_all <- grid.arrange(grobs = plots, nrow = 1)
  
  return(p_all)
}

get_rho_scatter_group_plots <- function(read_params, correlations, group = NULL) {
  plots <- list()
  
  first = TRUE
  for (r in correlations) {
    read_params$MinRho <- r$MinRho
    read_params$MaxRho <- r$MaxRho
    read_params$Method <- r$Method
    data_files <- get_data(read_params)
    data <- data_files$data
    data_obj <- data_files$data_obj
    
    p_rho <- get_rho_scatter_plot(data, data_obj, TRUE)
    p_rho <- p_rho + ggtitle(paste0("Correlation ", r$MinRho, "-", r$MaxRho))
    
    p_rhoG1 <- get_rho_scatter_plot(data[group == 1, ], data_obj, TRUE, group = 1)
    p_rhoG2 <- get_rho_scatter_plot(data[group == 2, ], data_obj, TRUE, group = 2)
    RIRG1 <- c(get_RIR(cor(data[group == 1,]), r$MinRho, r$MaxRho), get_RIR(data_obj$Rhos[[1]], r$MinRho, r$MaxRho))
    RIRG2 <- c(get_RIR(cor(data[group == 2,]), r$MinRho, r$MaxRho), get_RIR(data_obj$Rhos[[2]], r$MinRho, r$MaxRho))
    
    stats_df <- get_rho_stats(data, data_obj)
    p_stats1 <- tableGrob(stats_df[,1:7], rows = NULL, theme = ttheme_default(base_size = 12))
    p_stats2 <- tableGrob(cbind(stats_df[,c(1,8:9)], data.frame("RIRG1" = RIRG1, "RIRG2" = RIRG2)), rows = NULL, theme = ttheme_default(base_size = 12))
    
    plots[[toString(r)]] <- grid.arrange(p_rho, p_rhoG1, p_rhoG2, p_stats1, p_stats2, layout_matrix = cbind(c(1, 1, 2, 2, 3, 3, 4, 5)))
  }
  
  p_all <- grid.arrange(grobs = plots, nrow = 1)
  
  return(p_all)
}

# Get Correlation density and boxplots
get_rho_density_boxplots <- function(read_params, correlations) {
  RhoStrings <- unlist(lapply(correlations, function(x) paste0("[", sprintf("%.1f", x$MinRho), "-", sprintf("%.1f", x$MaxRho), "]")))
  
  n_rho <- length(correlations)
  stats <- list("All" = vector("list", n_rho), "G1" = vector("list", n_rho), "G2" = vector("list", n_rho))
  df <- data.frame(matrix(vector(), nrow = 0, ncol = 3,
                          dimnames = list(c(), c("RhoValues", "Rho", "Group"))))

  for (i in 1:n_rho) {
    read_params$MinRho <- correlations[[i]]$MinRho
    read_params$MaxRho <- correlations[[i]]$MaxRho
    read_params$Method <- correlations[[i]]$Method
    
    data_files <- get_data(read_params)
    data <- data_files$data
    data_obj <- data_files$data_obj
    
    ns <- data_obj$GenerationParameter[[3]][["ns"]]
    
    # Get Stats
    df_rho <- data.frame("Rho" = rep(RhoStrings[i], 2))
    stats[["All"]][[i]] <- rbind(cbind(df_rho, get_rho_stats(data, data_obj)), "-")
    stats[["G1"]][[i]] <- rbind(cbind(df_rho, get_rho_stats(data[1:ns[1],], data_obj, group = 1)), "-")
    stats[["G2"]][[i]] <- rbind(cbind(df_rho, get_rho_stats(data[(ns[1]+1):(ns[1]+ns[2]),], data_obj, group = 2)), "-")
    
    # Get Correlations
    corData <- cor(data); corData <- corData[upper.tri(corData)]
    corData1 <- cor(data[1:ns[1],]); corData1 <- corData1[upper.tri(corData1)]
    corData2 <- cor(data[(ns[1]+1):(ns[1]+ns[2]),]); corData2 <- corData2[upper.tri(corData2)]
    
    df <- rbind(df, data.frame("RhoValues" = c(corData, corData1, corData2), 
                               "Rho" = RhoStrings[i], 
                               "Group" = rep(c("Gesamt", "Gruppe 1", "Gruppe 2"), each = length(corData))))
  }
  
  df$Group <- factor(df$Group, levels = c("Gesamt", "Gruppe 1", "Gruppe 2"))
  df$Rho <- factor(df$Rho, levels = rev(RhoStrings))
  suppressWarnings(rho_breaks <- unique(na.omit(as.numeric(unlist(correlations, use.names = FALSE)))))
  
  # Draw plot
  plt <- ggplot(df, aes(x = RhoValues, fill = Rho)) + 
    annotate(geom="rect", xmin = c(0, 0.4, 0.7), xmax = c(0.3, 0.6, 0.9), ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.1) + 
    geom_boxplot(aes(y = -5), width = 7, outlier.size = 0.5, position=position_dodge(width=10)) + 
    geom_density(alpha = 0.8) +
    scale_fill_manual(values = RHO_COLORS) + 
    scale_x_continuous(breaks = unique(c(0, 1, rho_breaks)), limits = c(min(df$RhoValues), 1)) +
    scale_y_continuous(breaks = c(0, 5, 10)) +
    xlab("Paarweise Korrelationen") + ylab("") + facet_grid(rows = vars(Group)) +
    geom_vline(xintercept = rho_breaks, col = "black", alpha = 0.5) +
    theme(legend.position = "none", panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), text = element_text(size = FONT_SIZE_PLOT),
          strip.text = element_text(face = "bold"))
  
  grid_text <- textGrob(paste0("Szenario ", read_params$Scenario, 
                               "\np = ", read_params$DiffVars,
                               "\nseed = ", data_obj$Seed), 
                        gp = gpar(fontsize = FONT_SIZE_TEXT))
  grid_legend <- get_shared_boxplot_legend(RhoStrings, RHO_COLORS, font_size = FONT_SIZE_TEXT)
  

  create_tg <- function(x) {
    df <- do.call("rbind", x)[c(1, 4, 7),c(1, 9, 10)]
    colnames(df) <- c("Rho", "RMSR", "ImIntervall")
    
    tg <- tableGrob(df, rows =  NULL, 
                    theme = ttheme_default(base_size = FONT_SIZE_TABLE, padding = unit(c(10, 10), "mm")))
    return(tg)
  }
  
  grid_count <- textGrob(paste0(sum(1:(read_params$DiffVars-1)), " paarw. Korrelationen"),
                         gp=gpar(fontsize = FONT_SIZE_TEXT))
  grid_stats <- grid.arrange(grobs = c(lapply(stats, create_tg), list(grid_count)), 
                             layout_matrix = cbind(c(1,1,1, 2,2,2, 3,3,3, 4)))

  p_all <- grid.arrange(grobs = list(plt, grid_text, grid_legend, grid_stats), 
                        layout_matrix = rbind(c(1, 1, 4), c(1, 1, 4), c(1, 1, 4), c(1, 1, 4),
                                              c(2, 3, 3)))
  return(p_all)
}



get_rho_bounds_plot <- function(data_obj) {
  boundsG1 <- read_rho_bounds(data_obj, 1)
  boundsG2 <- read_rho_bounds(data_obj, 2)
  lbG1 <- boundsG1$LowerBound[lower.tri(boundsG1$LowerBound)]
  ubG1 <- boundsG1$UpperBound[lower.tri(boundsG1$UpperBound)]
  lbG2 <- boundsG2$LowerBound[lower.tri(boundsG2$LowerBound)]
  ubG2 <- boundsG2$UpperBound[lower.tri(boundsG2$UpperBound)]
  
  df <- data.frame(x = rep(1:length(lbG1), 2),
                   "LowerBound" = c(lbG1, lbG2),
                   "UpperBound" = c(ubG1, ubG2),
                   "Group" = rep(c(1,2), each = length(lbG1)))
  
  p_bounds <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = LowerBound, col = "orange")) +
    geom_line(aes(y = UpperBound, col = "blue")) +
    ylim(c(-1, 1)) + facet_wrap(vars(Group), ncol = 2)
  
  statsG1 <- data.frame("Bound" = c("Upper", "Lower"), 
                        rbind(round(summary(ubG1), 3),
                              round(summary(lbG1), 3)))
  statsG2 <- data.frame("Bound" = c("Upper", "Lower"), 
                        rbind(round(summary(ubG2), 3),
                              round(summary(lbG2), 3)))
  
  sc <- data_obj$GenerationParameter[[1]]$Scenario
  p <- data_obj$GenerationParameter[[3]]$DiffVars
  
  p_statsG1 <- tableGrob(statsG1, rows = NULL, theme = ttheme_default(base_size = 8))
  p_statsG2 <- tableGrob(statsG2, rows = NULL, theme = ttheme_default(base_size = 8))
  p_all <- grid.arrange(p_bounds, p_statsG1, p_statsG2, layout_matrix = rbind(c(1, 1),
                                                                              c(1, 1),
                                                                              c(1, 1),
                                                                              c(1, 1),
                                                                              c(2, 3)), top = paste0("Scenario ", sc, ", p = ", p))
  
  return(p_all)
}

read_rho_bounds <- function(data_obj, group = 1) {
  gp <- data_obj$GenerationParameter
  
  name <- paste0("sc",gp[[1]]$Scenario,"G",group,
                 "_v",gp[[3]]$EqualVars,"-",gp[[3]]$DiffVars, 
                 get_change_string(gp[[1]]$ChangeMeans, gp[[1]]$ChangeDispersions, gp[[1]]$ChangeZProbs),
                 "_nbounds",gp[[2]]$NBounds,
                 "_s", data_obj$Seed)
  file <- paste0(PATH_CORRELATION_BOUNDS, name, ".Rds")
  
  if (file.exists(file = file)) {
    corr_bounds <- readRDS(file = file)
  } else {
    cat("Bounds not found")
  }
  
  return(corr_bounds)
}

#===============================================================================
# .------------------------.
# |     OTHER FUNCTIONS    |
# '------------------------'
#===============================================================================

get_variance_boxplots <- function(read_params, correlations, scale = FALSE) {
  RhoStrings <- unlist(lapply(correlations, function(x) paste0("[", sprintf("%.1f", x$MinRho), "-", sprintf("%.1f", x$MaxRho), "]")))

  df <- data.frame(matrix(vector(), nrow = 0, ncol = 3,
                          dimnames = list(c(), c("Varianz", "Rho", "Gruppe"))))
  
  n_rho <- length(correlations)
  stats <- list("All" = vector("list", n_rho), "G1" = vector("list", n_rho), "G2" = vector("list", n_rho))
  
  for (i in 1:n_rho) {
    read_params$MinRho <- correlations[[i]]$MinRho
    read_params$MaxRho <- correlations[[i]]$MaxRho
    read_params$Method <- correlations[[i]]$Method
    
    data_files <- get_data(read_params)
    data <- data_files$data
    data_obj <- data_files$data_obj
    ns <- data_obj$GenerationParameter[[3]][["ns"]]
    if (scale) data <- scale(data)
    
    var_data <- var(data); var_data <- var_data[upper.tri(var_data)]
    var_G1 <- var(data[1:250,]); var_G1 <- var_G1[upper.tri(var_G1)]
    var_G2 <- var(data[251:500,]); var_G2 <- var_G2[upper.tri(var_G2)]
    
    stats[["All"]][[i]] <- data.frame("Rho" = RhoStrings[i], rbind(round(summary(var_data), 0)))
    stats[["G1"]][[i]] <- data.frame("Rho" = RhoStrings[i], rbind(round(summary(var_G1), 0)))
    stats[["G2"]][[i]] <- data.frame("Rho" = RhoStrings[i], rbind(round(summary(var_G2), 0)))
    
    df <- rbind(df, data.frame("Varianz" = c(var_data, var_G1, var_G2), 
                               "Rho" = RhoStrings[i], 
                               "Gruppe" = rep(c("Gesamt", "Gruppe 1", "Gruppe 2"), each = length(var_data))))

  }

  df$Gruppe <- factor(df$Gruppe, levels = c("Gesamt", "Gruppe 1", "Gruppe 2"))
  df$Rho <- factor(df$Rho, levels = rev(RhoStrings))

  plt <- ggplot(df) + geom_boxplot(aes(x = Varianz, fill = Rho)) + 
    scale_fill_manual(values = RHO_COLORS) + 
    xlab("Paarweise Varianzen") + ylab("") + facet_grid(rows = vars(Gruppe)) +
    theme(legend.position = "none", text = element_text(size = FONT_SIZE_PLOT),
          strip.text = element_text(face = "bold"),
          axis.ticks.y = element_blank(), axis.text.y = element_blank())

  grid_text <- textGrob(paste0("Szenario ", read_params$Scenario, 
                               "\np = ", read_params$DiffVars,
                               "\nseed = ", data_obj$Seed), 
                        gp = gpar(fontsize = FONT_SIZE_TEXT))
  grid_legend <- get_shared_boxplot_legend(RhoStrings, RHO_COLORS, font_size = FONT_SIZE_TEXT)
  
  
  create_tg <- function(x) {
    df <- do.call("rbind", x)[,c(1, 3, 5, 6)]
    colnames(df) <- c("Rho", "Un.Qu.", "Mittelw.", "Ob.Qu.")
    tg <- tableGrob(df, rows =  NULL, theme = ttheme_default(base_size = FONT_SIZE_TABLE))
    return(tg)
  }

  p_all <- grid.arrange(grobs = list(plt, grid_text, grid_legend), 
                        layout_matrix = rbind(c(1, 1, 1), c(1, 1, 1), c(1, 1, 1), c(1, 1, 1),
                                              c(2, 3, 3)))

  return(p_all)
}
