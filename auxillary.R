#===============================================================================
# CONSTANTS
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
library(GGally)

PATH_DATA <- "data/"
PATH_DATA_OBJ <- paste0(PATH_DATA, "datasets/")
PATH_CORRELATION_MATRIX <- paste0(PATH_DATA, "correlation_matrices/")
PATH_CORRELATION_BOUNDS <- paste0(PATH_DATA, "correlation_bounds/")
PATH_QUANTILES <- paste0(PATH_CORRELATION_BOUNDS, "quantiles/")

PATH_GRAPHICS <- "graphics/"
PATH_ARI <- "ari/"

TRUE_COLORS <- c("1" = "#ff9500", "2" = "#002e59")
# CLUST_COLORS <- scale_colour_manual(values = c("1" = "#ffb514", "2" = "#004c94"))
RHO_COLORS <- c("[0.0-0.3]" = "#C7E9C0", "[0.4-0.6]" = "#74C476", "[0.7-0.9]" = "#006D2C")
P_COLORS <- c("250" = "#c9df8a", "500" = "#77ab59", "1000" = "#36802d", "2000" =  "#234d20")

THEME_BLANK <- theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank())
THEME_STRIP_BLANK <- theme(strip.background = element_blank(), strip.text = element_blank())

FONT_SIZE_TABLE <- 30
FONT_SIZE_TEXT <- 40
FONT_SIZE_PLOT <- 40

DPI <- 144
GRAPH_HEIGHT <- 1080/72
GRAPH_WIDTH <- 1920/72

#===============================================================================
# INCLUDED FUNCTIONS
#
# -get_data
# -scale2ab
# -get_change_string
# -get_rho_method_string
# -print_lists
# -print_list
#
# -compare_columns
#
# -get_shared_point_legend
# -get_shared_boxplot_legend
# -get_shared_boxplot_legend_p
#
#===============================================================================

#Create needed subfolders if they dont exists
create_folders <- function() {
  dir.create(file.path(PATH_DATA), showWarnings = FALSE)
  dir.create(file.path(PATH_DATA_OBJ), showWarnings = FALSE)
  dir.create(file.path(PATH_CORRELATION_MATRIX), showWarnings = FALSE)
  dir.create(file.path(PATH_CORRELATION_BOUNDS), showWarnings = FALSE)
  dir.create(file.path(PATH_QUANTILES), showWarnings = FALSE)
  
  dir.create(file.path(PATH_GRAPHICS), showWarnings = FALSE)
  dir.create(file.path(PATH_GRAPHICS, "umap"), showWarnings = FALSE)
  dir.create(file.path(PATH_GRAPHICS, "tsne"), showWarnings = FALSE)
  dir.create(file.path(PATH_GRAPHICS, "pca"), showWarnings = FALSE)
  dir.create(file.path(PATH_GRAPHICS, "boxplots"), showWarnings = FALSE)
  
  dir.create(file.path(PATH_ARI), showWarnings = FALSE)
}

#===============================================================================
# .--------------------------.
# |     VARIOUS FUNCTIONS    |
# '--------------------------'
#===============================================================================
get_data <- function(read_params) {
  name <- paste0("sc", read_params$Scenario,
                 "_n", paste(read_params$ns, collapse="-"),
                 "_v", read_params$EqualVars, "-", read_params$DiffVars, 
                 get_change_string(read_params$ChangeMeans, read_params$ChangeDispersions, read_params$ChangeZProbs),
                 "_r", get_rho_method_string(read_params$Method), ifelse(read_params$RespectBounds, "-RB", ""), 
                 read_params$MinRho, "-", read_params$MaxRho,
                 "_s", read_params$Seed)

  data_obj <- readRDS(file=paste0(PATH_DATA_OBJ, name, ".Rds"))
  data <- data_obj$Data

  return(list("data" = data, "data_obj" = data_obj))
}

scale2ab <- function(x, a, b) {
  xmin <- min(x)
  xmax <- max(x)
  x <- (((b - a)*(x - xmin))/(xmax-xmin)) + a
  return(x)
}

get_change_string <- function(change_means = TRUE, change_dispersion = FALSE, change_zprobs = FALSE) {
  ret <- paste0(ifelse(change_means, "T", "F"), ifelse(change_dispersion, "T", "F"), ifelse(change_zprobs, "T", "F"))
  return(ret)
}

get_rho_method_string <- function(rho_method) {
  if (rho_method == "Gauss") {
    return("G")
  } else if (rho_method == "Uniform") {
    return("U")
  } else if (rho_method == "Beta") {
    return("B")
  }
}

print_lists <- function(lists) {
  df <- lapply(lists, function(x) data.frame(lapply(x, paste0, collapse = "-")))
  df <- do.call("rbind", df)
  print(df)
}

print_list <- function(list) {
  df <- data.frame(lapply(list, paste0, collapse = "-"))
  print(df)
}

compare_columns <- function(data, data_obj, columns, group = NULL) {
  if (is.null(group)) {group <- get_group_split(data, data_obj)}
  df <- cbind(data[,columns], "Group" = group)
  p_cols <- ggpairs(df, mapping = aes(colour = Group, alpha = 0.4))
  return(p_cols)
}


get_shared_point_legend <- function(groups, scale_color, font_size = 8){
  # Function adapted from:
  # https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots/13650878#13650878
  groups <- paste0(groups, "    ")
  names(scale_color) <- paste0(names(scale_color), "    ")
  plt <- ggplot(data.frame(x = rep(3, length(groups)), y = rep(4, length(groups)), Gruppe = groups), aes(x = x, y = y, col = Gruppe)) + 
    geom_point(size = 16) + scale_color_manual(values = scale_color) + 
    theme_minimal(base_size = font_size) +
    theme(legend.title=element_text(size=font_size), legend.text=element_text(size=font_size), legend.position = "bottom") 
  # guides(color = guide_legend(override.aes = list(size = 4)))
  
  tmp <- ggplot_gtable(ggplot_build(plt))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

get_shared_boxplot_legend <- function(groups, scale_fill, font_size = 8){
  # Function adapted from:
  # https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots/13650878#13650878
  groups <- factor(paste0(groups, "    "), levels = paste0(groups, "    "))
  names(scale_fill) <- paste0(names(scale_fill), "    ")
  plt <- ggplot(data.frame(x = rep(3, length(groups)), y = rep(4, length(groups)), Rho = groups), aes(x = x, y = y, fill = Rho)) + 
    geom_point(shape = 22, size = 20) + scale_fill_manual(values = scale_fill) + 
    theme_minimal(base_size = font_size) +
    theme(legend.title=element_text(size=font_size), legend.text=element_text(size=font_size), 
          legend.position = "bottom")
  
  # guides(fill = guide_legend(override.aes = list(size = 1)))
  
  tmp <- ggplot_gtable(ggplot_build(plt))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

get_shared_boxplot_legend_p <- function(groups, scale_fill, font_size = 8){
  # Function adapted from:
  # https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots/13650878#13650878
  groups <- factor(paste0(groups, "    "), levels = paste0(groups, "    "))
  names(scale_fill) <- paste0(names(scale_fill), "    ")
  plt <- ggplot(data.frame(x = rep(3, length(groups)), y = rep(4, length(groups)), p = groups), aes(x = x, y = y, fill = p)) + 
    geom_point(shape = 22, size = 20) + scale_fill_manual(values = scale_fill) + 
    theme_minimal(base_size = font_size) +
    theme(legend.title=element_text(size=font_size), legend.text=element_text(size=font_size), 
          legend.position = "bottom")
  
  # guides(fill = guide_legend(override.aes = list(size = 1)))
  
  tmp <- ggplot_gtable(ggplot_build(plt))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}