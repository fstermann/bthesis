# Set Working Directory
# Set this to your desired Path where the files below are located
setwd("R")

# Load External Functions
source("auxillary.R")
create_folders() # Create the needed subfolders if they don't exist

source("data_generation.R")
source("dim_red.R")
source("clustering.R")
source("correlations.R")

# .---------------------------.
# |         SCENARIOS         |
# '---------------------------'
# - Scnenario: Scenario ID, ["String"]
# - MinMeans: Vector of minimal value of mean for the ZINB Distribution for each group, e.g. c(10, 20), [1, oo]
# - MaxMeans: Vector of maximal value of mean for the ZINB Distribution for each group, e.g. c(100, 200), [1, oo]
# - MinDispersions: Vector of minimal value of dispersion for the ZINB Distribution for each group, e.g. c(0.3, 0.4), [0, oo]
# - MaxDispersions: Vector of maximal value of dispersion for the ZINB Distribution for each group, e.g. c(0.7, 0.8), [0, oo]
# - MinZProbs: Vector of minimal value of zero probability for the ZINB Distribution for each group, e.g. c(0.001, 0.0001), [0, 1]
# - MaxZProbs: Vector of maximal value of zero probability for the ZINB Distribution for each group, e.g. c(0.1, 0.08), [0, 1]
# - ChangeMeans: Change the mean for each group in specified "DiffVars" variables
# - ChangeDispersions: Change the dispersion for each group in specified "DiffVars" variables
# - ChangeZProbs: Change the zero probability  for each group in specified "DiffVars" variables

scenarios <- list()
scenarios[[1]] <- list(Scenario = "01", 
                       MinMeans = c(45, 12), MaxMeans = c(293, 112), 
                       MinDispersions = c(0.27, 0.27), MaxDispersions = c(0.47, 0.47), 
                       MinZProbs = c(5.30*10^-7, 4.85*10^-7), MaxZProbs = c(0.01, 2.11*10^-5),
                       ChangeMeans = TRUE, ChangeDispersions = TRUE, ChangeZProbs = TRUE)
scenarios[[2]] <- list(Scenario = "02", 
                       MinMeans = c(27, 6), MaxMeans = c(397, 171), 
                       MinDispersions = c(0.24, 0.23), MaxDispersions = c(0.55, 0.55), 
                       MinZProbs = c(3.65*10^-7, 3.26*10^-7), MaxZProbs = c(0.04, 2.91*10^-2),
                       ChangeMeans = TRUE, ChangeDispersions = TRUE, ChangeZProbs = TRUE)
scenarios[[3]] <- list(Scenario = "03", 
                       MinMeans = c(19, 2), MaxMeans = c(576, 217), 
                       MinDispersions = c(0.18, 0.17), MaxDispersions = c(0.78, 0.82), 
                       MinZProbs = c(2.28*10^-7, 2.18*10^-7), MaxZProbs = c(0.08, 6.11*10^-2),
                       ChangeMeans = TRUE, ChangeDispersions = TRUE, ChangeZProbs = TRUE)

# .------------------------------.
# |         CORRELATIONS         |
# '------------------------------'
# - MinRho: Minimal Correlation, [-1, 1]
# - MaxRho: Maximal Correlation,  [-1, 1]
# - Method: Distribution of the pairwise Correlations, ["Uniform", "Gauss", "Beta"]
# - Iterations: Maximal number of iterations in trying to find a valid Correlation Matrix, [1, oo]
# - RespectBounds: Correct the target Correlation, if it lies out of the Correlation Bounds [TRUE, FALSE]
# - NBounds: Length of the Uniform vector used to calculate the Bounds, [1, oo]

rhos <- list()
rhos[[1]] <- data.frame(MinRho = 0, MaxRho = 0.3, Method = "Gauss", Iterations = 10, RespectBounds = TRUE, NBounds = 100000)
rhos[[2]] <- data.frame(MinRho = 0.4, MaxRho = 0.6, Method = "Gauss", Iterations = 10, RespectBounds = TRUE, NBounds = 100000)
rhos[[3]] <- data.frame(MinRho = 0.7, MaxRho = 0.9, Method = "Gauss", Iterations = 10, RespectBounds = TRUE, NBounds = 100000)

# .--------------------------------------------.
# |         NUMBER OF ROWS AND COLUMNS         |
# '--------------------------------------------'
# - ns: Vector of group sizes, e.g. c(250, 250) (Length has to match parameters in scenarios)
# - EqualVars: Number of Variables with equal Distribution parameters [0, oo]
# - DiffVars: Number of Variables with different Distribution parameters [0, oo] (p = EqualVars + DiffVars)

sizes <- list()
sizes[[1]] <- list(ns = c(250, 250), EqualVars = 0, DiffVars = 250)
sizes[[2]] <- list(ns = c(250, 250), EqualVars = 0, DiffVars = 500)
sizes[[3]] <- list(ns = c(250, 250), EqualVars = 0, DiffVars = 1000)
sizes[[4]] <- list(ns = c(250, 250), EqualVars = 0, DiffVars = 2000)



# Choose the desired scenario, rho and size values
i_sc <- c(1:3)
i_rho <- c(1:3)
i_sz <- c(1:4)

print_lists(scenarios[i_sc])
print_lists(rhos[i_rho])
print_lists(sizes[i_sz])


#===============================================================================
# .-----------------------------------.
# |         GENERATE DATASETS         |
# '-----------------------------------'
#===============================================================================

# If calculated for new datasets without existing correlation bounds, 
# this takes a long time. 
seeds <- 30091997 + seq(0, 90, by = 10)

for (seed in seeds) {
  generate_datasets(scenarios = scenarios[i_sc], rhos = rhos[i_rho], sizes = sizes[i_sz], seed = seed)
}

#===============================================================================
# .-----------------------------------------------.
# |             DEFINE PLOT PARAMETER             |
# '-----------------------------------------------'
#===============================================================================
seeds <- 30091997 + seq(0, 90, by = 10)

all_correlations <- list(list("MinRho" = 0, "MaxRho" = 0.3, "Method" = "Gauss"),
                         list("MinRho" = 0.4, "MaxRho" = 0.6, "Method" = "Gauss"), 
                         list("MinRho" = 0.7, "MaxRho" = 0.9, "Method" = "Gauss"))

all_scenarios <- c("01", "02", "03")
all_ps <- c(250, 500, 1000, 2000)

read_params <- list(Scenario = "01",
                    MinRho = 0.0, MaxRho =  0.3, Method = "Gauss", RespectBounds = TRUE,
                    ns = c(250, 250), EqualVars = 0, DiffVars = 250, 
                    ChangeMeans = TRUE, ChangeDispersions = TRUE, ChangeZProbs = TRUE,
                    Seed = 30091997)

dimred_params <- list(scale_embedding = TRUE, scale_data = FALSE,
                      n_neighbors = 15, min_dist = 0.01, 
                      tsne_iter = 1000, perplexity = 30,
                      pca_center = TRUE, pca_scale = TRUE)

group <- factor(c(rep(1, times = read_params$ns[1]), rep(2, times = read_params$ns[2])))

#===============================================================================
# .--------------------------------------.
# |             THESIS PLOTS             |
# |              AND TABLES              |
# '--------------------------------------'
#===============================================================================

# Save example Rho Density/Boxplot
read_params$DiffVars <- 2000
p_rho <- get_rho_density_boxplots(read_params, all_correlations)
ggsave(paste0(PATH_GRAPHICS, "rho_p2000.png"), p_rho, dpi = DPI, width = GRAPH_WIDTH, height = GRAPH_HEIGHT)


# Scatter Plots for Dimensionaly Reduction and Clustering
for (sc in all_scenarios) {
  read_params$Scenario <- sc
  for (p in all_ps) {
    read_params$DiffVars <- p
    cat(paste0("\r", rep(" ", 50)))
    cat("Producing plots for Scenario: ", sc, ", p: ", p, ", seed: ", read_params$Seed, ", DimRed: ", sep = "")
    
    cat("UMAP")
    dimred_params$scale_embedding <- FALSE
    p_umap <- show_CL_correlations("UMAP", dimred_params, read_params, all_correlations,
                                   group = group, seed = read_params$Seed)
    
    cat("/ PCA ")
    dimred_params$scale_embedding <- TRUE
    p_pca <- show_CL_correlations("PCA", dimred_params, read_params, all_correlations,
                                  group = group, seed = read_params$Seed)
    
    cat("/ t-SNE")
    p_tsne <- show_CL_correlations("t-SNE", dimred_params, read_params, all_correlations, 
                                   group = group, seed = read_params$Seed)
    
    ggsave(paste0(PATH_GRAPHICS, "umap/umap_sc",sc,"_p",p,"_s",read_params$Seed,".pdf"), p_umap, dpi = DPI, width = GRAPH_WIDTH, height = GRAPH_HEIGHT)
    ggsave(paste0(PATH_GRAPHICS, "pca/pca_sc",sc,"_p",p,"_s",read_params$Seed,".pdf"), p_pca, dpi = DPI, width = GRAPH_WIDTH, height = GRAPH_HEIGHT)
    ggsave(paste0(PATH_GRAPHICS, "tsne/tsne_sc",sc,"_p",p,"_s",read_params$Seed,".pdf"), p_tsne, dpi = DPI, width = GRAPH_WIDTH, height = GRAPH_HEIGHT)
  }
}

# Get ARI DataFrames
for (sc in all_scenarios) {
  read_params$Scenario <- sc
  df <- get_ARI_ForSeeds(seeds = seeds, read_params = read_params, dimred_params = dimred_params, 
                         correlations = all_correlations, ps = all_ps, 
                         dimreds = c("UMAP", "t-SNE", "PCA"), clusterings = c("kmeans", "hclust", "mclust"), iterations = 10)
  saveRDS(df, file = paste0(PATH_ARI, "tsne_ari_sc", sc, ".Rds"))
}

df1 <- readRDS(file=paste0("archiv/", PATH_ARI, "ari_sc01.Rds"))
df2 <- readRDS(file=paste0("archiv/", PATH_ARI, "ari_sc02.Rds"))
df3 <- readRDS(file=paste0("archiv/", PATH_ARI, "ari_sc03.Rds"))


# Get ARI Boxplots
# Exclude the "AVG", "Min", "Max"
plots_DR <- get_ARI_boxplots_byDR(df1[1:(length(df1)-3)], df2[1:(length(df2)-3)], df3[1:(length(df3)-3)])
for (dr in c("UMAP", "t-SNE", "PCA")) {
  ggsave(paste0(PATH_GRAPHICS, "boxplots/ari_boxplots_", dr, ".pdf"), plots_DR[[toString(dr)]], dpi = DPI, width = GRAPH_WIDTH, height = GRAPH_HEIGHT)
}

# Get ARI Tables 
df1$Average
df2$Average
df3$Average



# Additional functions which are not presented in the paper.
# One can visualize some methods with this and play around with
# different parameter.

#===============================================================================
# .------------------------------------------.
# |             EXPLORE DATASETS             |
# '------------------------------------------'
#===============================================================================

read_params$Scenario <- "01"
read_params$DiffVars <- 250

data_files <- get_data(read_params)
data <- data_files$data
data_obj <- data_files$data_obj

# --> Parameter
data_obj$Seed
data_obj$ComputationTime

data_obj$VariableParameters
print_list(data_obj$GenerationParameter[[1]])
print_list(data_obj$GenerationParameter[[2]])
print_list(data_obj$GenerationParameter[[3]])

data_obj$Rhos
data_obj$Sigmas

compare_columns(data, data_obj, 1:5, group = group)

#===============================================================================
# .----------------------------------------------.
# |             EXPLORE CORRELATIONS             |
# '----------------------------------------------'
#                 correlations.R
#===============================================================================

# Pariwise correlation stats
get_rho_stats(data, data_obj)

# Pairwise Correlation scatter plot to Target correlation
get_rho_scatter_plot(data, data_obj)

# Visualize correlation bounds
get_rho_bounds_plot(data_obj)


# Pairwise Correlation scatter plot to Target correlation for
# multiple correlations
get_rho_scatter_plots(read_params, all_correlations)
# Additionally split by Group
get_rho_scatter_group_plots(read_params, all_correlations, group = group)

# Visualize distribution for multiple correlations
get_rho_density_boxplots(read_params, all_correlations)


# Visualize the Variance in datasets
# Optionally scale the data beforehand (PCA input)
get_variance_boxplots(read_params, all_correlations, scale = TRUE)


#===============================================================================
# .-----------------------------------------------------.
# |             EXPLORE DIMENSION REDUCTION             |
# '-----------------------------------------------------'
#                       dim_red.R
#===============================================================================

# Dimension Reduction plots for all correlations
show_umap_correlations(read_params, all_correlations, n_neighbors = dimred_params$n_neighbors, min_dist = dimred_params$min_dist, group = group)
show_tsne_correlations(read_params, all_correlations, tsne_iter = dimred_params$tsne_iter, perplexity = dimred_params$perplexity, group = group)
show_pca_correlations(read_params, all_correlations, center = dimred_params$pca_center, scale. = dimred_params$pca_scale, group = group)

# All DimReds in one plot
show_all_dim_red(data, data_obj, group = group)

# All Correlations in one plot
show_all_correlations(read_params, all_correlations, group = group)

# All Scenarios in one plot
show_all_scenarios(read_params, all_scenarios, group = group)

# Try out different UMAP parameter
show_umap_parameters(read_params, n_neighbors_range = c(2, 4, 8, 16, 32, 64), min_dist_range = c(0, 0.25, 0.5, 0.75, 1), group = group)


#===============================================================================
# .----------------------------------------------------.
# |             EXPLORE CLUSTERING METHODS             |
# '----------------------------------------------------'
#                       clustering.R
#===============================================================================

# Plots as produced in section THESIS PLOTS
show_CL_correlations("UMAP", dimred_params, read_params, all_correlations, group = group, seed = read_params$Seed)
