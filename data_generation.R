source("auxillary.R")

library(VGAM)
library(Matrix)
library(parallel)
library(doSNOW)

#===============================================================================
# INCLUDED FUNCTIONS
#
# -generate_datasets
# --generate_corr_data
# ---calculate_data
# 
# -get_rho
# --generate_rho
# --correct_rho
# -get_rho_bounds
# --calculate_rho_bounds
# ---get_quantiles
# ---calculate_quantiles
#
# -calculate_sigma
# -draw_X
# 
# -get_params
#
#
# Most of these methods are adapted from the Package SimCorrMix by Allison Fialkowski:
# https://cran.r-project.org/package=SimCorrMix
# Some calculations for the specific ZINB-Distribution have been adjusted.
#===============================================================================


#===============================================================================
# .-------------------------------.
# |     GENERAL DATA FUNCTIONS    |
# '-------------------------------'
#===============================================================================

generate_datasets <- function(scenarios, rhos, sizes, seed) {
  max_ds <- length(scenarios)*length(rhos)*length(sizes)
  ds <- 0
  for (i in 1:length(scenarios)) {
    for (j in 1:length(rhos)) {
      for (k in 1:length(sizes)) {
        ds <- ds +1
        sc <- scenarios[[i]]
        rh <- rhos[[j]]
        sz <- sizes[[k]]
        
        cat("\n\n[",ds,"/",max_ds,"] Simulating data for params:\n")
        print_list(sc)
        print_list(rh)
        print_list(sz)
        cat("\n")
        
        ret <- generate_corr_data(ns = sz$ns, vars_equal = sz$EqualVars, vars_diff = sz$DiffVars, 
                                      min_means = sc$MinMeans, max_means = sc$MaxMeans, 
                                      min_dispersions = sc$MinDispersions, max_dispersions = sc$MaxDispersions, 
                                      min_zprobs = sc$MinZProbs, max_zprobs = sc$MaxZProbs, 
                                      min_rho = rh$MinRho, max_rho = rh$MaxRho, rho_method = rh$Method, rho_iter = rh$Iterations, 
                                      respect_bounds = rh$RespectBounds,
                                      change_means = sc$ChangeMeans, change_dispersion = sc$ChangeDispersions, change_zprobs = sc$ChangeZProbs,
                                      seed = seed, scenario = sc$Scenario, nbounds = rh$NBounds)

        # Save whole dataset and subset objects + params
        name <- paste0("sc", sc$Scenario,
                       "_n", paste(sz$ns, collapse="-"),
                       "_v", sz$EqualVars, "-", sz$DiffVars, get_change_string(sc$ChangeMeans, sc$ChangeDispersions, sc$ChangeZProbs),
                       "_r", get_rho_method_string(rh$Method), ifelse(rh$RespectBounds, "-RB", ""), rh$MinRho, "-", rh$MaxRho, 
                       "_s", seed)
        cat("Saving data as", name,"\n")
        
        saveRDS(list(Data = ret$Data, VariableParameters = ret$VariableParameters, 
                     GenerationParameter = list(sc, rh, sz),
                     Rho = ret$Rho, Rhos = ret$Rhos, Sigmas = ret$Sigmas,
                     ComputationTime = ret$ComputationTime, Seed = seed), 
                file=paste0(PATH_DATA_OBJ, name, ".Rds"))
      }
    }
  }
}

generate_corr_data <- function(ns, vars_equal, vars_diff, 
                                   min_means, max_means,
                                   min_dispersions, max_dispersions, 
                                   min_zprobs, max_zprobs,
                                   min_rho, max_rho, rho_method = "Gauss", rho_iter = 10, respect_bounds = FALSE,
                                   change_means = TRUE, change_dispersion = TRUE, change_zprobs = TRUE,
                                   seed = 1234, scenario = "01", nbounds = 100000) {
  
  start_time <- Sys.time()
  set.seed(seed)
  
  # Parameter
  n <- sum(ns)
  n_vars <- vars_equal + vars_diff
  n_groups <- length(ns)
  
  # Draw Parameter
  # Means
  means <- get_params(n_equal = ifelse(change_means, vars_equal, n_vars), n_diff = ifelse(change_means, vars_diff, 0), 
                      mins = min_means, maxs = max_means, digits = 0, seed = seed)
  # Dispersion Parameter
  dispersions <- get_params(n_equal = ifelse(change_dispersion, vars_equal, n_vars), n_diff = ifelse(change_dispersion, vars_diff, 0), 
                            mins = min_dispersions, maxs = max_dispersions, digits = 2, seed = seed)
  # Probability Zero Inflated Structural Zeros
  zprobs <- get_params(n_equal = ifelse(change_zprobs, vars_equal, n_vars), n_diff = ifelse(change_zprobs, vars_diff, 0), 
                       mins = min_zprobs, maxs = max_zprobs, digits = 5, seed = seed)
  
  # Correlation Matrix
  Rho <- get_rho(min_rho = min_rho, max_rho = max_rho, n_vars = n_vars, n_iter = rho_iter, method = rho_method, seed = seed)

  # Save Parameter in DataFrame
  params <- data.frame(Var = sprintf("X%s", seq(1:n_vars)))
  for (g in 1:n_groups) {
    params[paste0("MeansG",g)] <- means[[g]]
  } 
  for (g in 1:n_groups) {
    params[paste0("DispG",g)] <- dispersions[[g]]
  } 
  for (g in 1:n_groups) {
    params[paste0("ZProbG",g)] <- zprobs[[g]]
  } 
  
  Rhos <- list()
  Sigmas <- list()
  datasets <- list()
  
  # Generate Datasets for each group
  for (g in 1:n_groups) {
    Rhos[[g]] <- Rho
    # Calculate Correlation Bounds
    Rho_bounds <- get_rho_bounds(ns = ns[g], vars_equal = vars_equal, vars_diff = vars_diff, sc = paste0(scenario, "G", g), 
                                 means = means[[g]], dispersions = dispersions[[g]], zprobs = zprobs[[g]], 
                                 nbounds = nbounds, 
                                 change_means = change_means, change_dispersion = change_dispersion, change_zprobs = change_zprobs, 
                                 seed = seed)
    
    if (respect_bounds) {Rhos[[g]] <- correct_rho(rho = Rhos[[g]], rho_lowerbound = Rho_bounds$LowerBound, rho_upperbound = Rho_bounds$UpperBound)}

    # Calculate Sigma
    Sigmas[[g]] <- calculate_sigma(rho_lowerbound = Rho_bounds$LowerBound, rho_upperbound = Rho_bounds$UpperBound, rho_target = Rhos[[g]])
    
    # Generate data
    datasets[[g]] <- calculate_data(ns = ns[g], n_vars = n_vars,
                                  means = means[[g]], dispersions = dispersions[[g]], zprobs = zprobs[[g]],
                                  rho = Rhos[[g]], Sigma = Sigmas[[g]], seed = seed)
  }
  
  data <- do.call("rbind", datasets)

  end_time <- Sys.time()
  comp_time <- difftime(end_time, start_time, units = "mins")[[1]]
  
  cat("Computation time:",comp_time,"mins\n")
  
  return(list("Data" = data, "VariableParameters" = params, "Rho" = Rho, "Rhos" = Rhos, "Sigmas" = Sigmas, 
              "ComputationTime" = comp_time))
}

calculate_data <- function(ns, n_vars, means, dispersions, zprobs,
                           rho = NULL, Sigma = NULL, seed = 1234) {
  # N and Variables
  n <- sum(ns)

  # Generate Data with NORTA-method
  X <- draw_X(n = n, sigma = Sigma, seed = seed)
  Y <- matrix(1, nrow = n, ncol = n_vars)

  for (i in 1:n_vars) {
    cat("Calculating variable", i, "values \r")
    Y[, i] <- qzinegbin(p = pnorm(X[, i]), munb = means[i], size = dispersions[i], pstr0 = zprobs[i])
  }
  
  cat("\n\n")
  
  Y <- data.frame(Y)
  
  return(Y)
}


#===============================================================================
# .----------------------.
# |     RHO FUNCTIONS    |
# '----------------------'
#===============================================================================

get_rho <- function(min_rho, max_rho, n_vars, n_iter = 10, method = "Uniform", seed = 1234) {
  name <- paste0("r", method, min_rho,"-",max_rho,
                 "_v",n_vars)
  file <- paste0(PATH_CORRELATION_MATRIX, name, ".Rds")
  
  if (file.exists(file=file)) {
    cat("Reading already existing target correlation matrix from",file,"\n")
    rho <- readRDS(file=file)
  } else {
    cat("Calculating new correlation matrix!\n")
    rho <- generate_rho(min_rho = min_rho, max_rho = max_rho, n_vars = n_vars, n_iter = n_iter, seed = seed, method = method)
    saveRDS(rho, file=file)
  }
  
  return(rho)
  # Returns Matrix
}

get_rho_bounds <- function(ns, vars_equal, vars_diff, sc, means, dispersions, zprobs, nbounds = 10000, 
                           change_means = TRUE, change_dispersion = FALSE, change_zprobs = FALSE, seed = 1234) {
  name <- paste0("sc",sc,
                 "_v",vars_equal,"-",vars_diff, get_change_string(change_means, change_dispersion, change_zprobs),
                 "_nbounds",nbounds,
                 "_s",seed)
  file <- paste0(PATH_CORRELATION_BOUNDS, name, ".Rds")
  
  if (file.exists(file = file)) {
    cat("Reading already existing correlation bounds from",file,"\n")
    rho_bounds <- readRDS(file = file)
  } else {
    quantiles <- get_quantiles(sc, vars_equal, vars_diff, 
                               get_change_string(change_means, change_dispersion, change_zprobs),
                               means, dispersions, zprobs, nbounds, seed)
    
    cat("Calculating new correlation bounds...")
    rho_bounds <- calculate_rho_bounds(quantiles, vars_equal + vars_diff, seed)
    saveRDS(rho_bounds, file=file)
  }
  
  return(rho_bounds)
  # Returns list with "UpperBound", "LowerBound", "ComputationTime", "Seed"
}

get_quantiles <- function(scenario, vars_equal, vars_diff, change_string, means, dispersions, zprobs, nbounds, seed) {
  # The name quantiles might be misleading. What's calculated are actually
  # the values of the quantile function, so we get real values of the 
  # ZINB-Distribution. So "quantiles" actually means the "quantile function".
  
  name <- paste0("sc",scenario,
                 "_",change_string,
                 "_nbounds",nbounds,
                 "_s",seed)
  file <- paste0(PATH_QUANTILES, name, ".Rds")
  
  n_vars <- vars_equal + vars_diff
    
  if (file.exists(file = file)) {
    cat("Reading already existing quantiles from", file, "\n")
    r_quantiles <- readRDS(file = file)
    
    if (length(r_quantiles[["U"]]) == n_vars) {
      # Use exisiting quantiles
      quantiles <- r_quantiles
    } else if (length(r_quantiles[["U"]]) > n_vars) {
      # Use part of existing quantiles
      quantiles <- list("U" = r_quantiles[["U"]][1:n_vars],
                        "1-U" = r_quantiles[["1-U"]][1:n_vars],
                        "ComputationTime" = r_quantiles$ComputationTime,
                        "Seed" = seed)
    } else if (length(r_quantiles$U) < n_vars) {
      # Use existing quantiles and calculate the missing ones
      quantiles <- calculate_quantiles(n_vars, means, dispersions, zprobs, 
                                       quantiles = r_quantiles, nbounds = nbounds, seed = seed)
      saveRDS(quantiles, file=file)
    }
  } else {
    cat("Calculating new quantiles!\n")
    quantiles <- calculate_quantiles(n_vars, means, dispersions, zprobs, 
                                     nbounds = nbounds, seed = seed)
    saveRDS(quantiles, file=file)
  }
  
  return(quantiles)
  # Returns list with "U", "1-U", "ComputationTime", "Seed"
}

calculate_quantiles <- function(n_vars, means, dispersions, zprobs, quantiles = NULL, nbounds = 10000, seed = 1234){
  start_time <- Sys.time()
  
  if (!is.null(quantiles)) {
    start_v <- length(quantiles[["U"]]) + 1
  } else {
    start_v <- 1
    quantiles <- list("U" = list(), "1-U" = list(), "ComputationTime" = 0, "Seed" = seed)
  }
  
  # Draw U Vector of size nbounds
  set.seed(seed)
  u <- runif(nbounds, 0, 1)
  
  # Parallelise the computation of quantile function
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  clusterExport(cl, list("qzinegbin", "u", "means", "dispersions", "zprobs"), envir=environment())
  
  # Add progress bar
  invisible(capture.output(pb <- txtProgressBar(max = n_vars-start_v+1, style = 3)))
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Calculate quantile function
  cat("Calculating U quantile functions... \n")
  qu <- foreach(i = start_v:n_vars, .options.snow = opts) %dopar%
    {return(qzinegbin(u, munb = means[i], size = dispersions[i], pstr0 = zprobs[i]))}
  
  cat("\nCalculating 1-U quantile functions... \n")
  q1u <- foreach(i = start_v:n_vars, .options.snow = opts) %dopar%
    {return(qzinegbin(1-u, munb = means[i], size = dispersions[i], pstr0 = zprobs[i]))}
  
  close(pb)
  stopCluster(cl)
  
  # Combine output
  quantiles[["U"]] <- c(quantiles[["U"]], qu)
  quantiles[["1-U"]] <- c(quantiles[["1-U"]], q1u)
  
  end_time <- Sys.time()
  comp_time <- difftime(end_time, start_time, units = "mins")[[1]]
  quantiles$ComputationTime <- quantiles$ComputationTime + comp_time
  
  return(quantiles)
  # Returns list with "U", "1-U", "ComputationTime", "Seed"
}

calculate_rho_bounds <- function(quantiles, n_vars, seed) {
  start_time <- Sys.time()

  lower_bound <- matrix(-1, n_vars, n_vars)
  upper_bound <- matrix(1, n_vars, n_vars)

  # This is the basic method.
  # The used code below produces the same output in less time.
  # for (i in 1:(n_vars-1)) {
  #   for (j in (i+1):n_vars) {
  # 
  #     maxcor <- cor(quantiles[["U"]][[i]], quantiles[["U"]][[j]])
  #     mincor <- cor(quantiles[["U"]][[i]], quantiles[["1-U"]][[j]])
  # 
  #     lower_bound[i, j] <- lower_bound[j, i] <- mincor
  #     upper_bound[i, j] <- upper_bound[j, i] <- maxcor
  #   }
  # }
  
  Matrix_U <- matrix(unlist(quantiles[["U"]]), ncol = n_vars)
  Matrix_1U <- matrix(unlist(quantiles[["1-U"]]), ncol = n_vars)
  
  cat("Lower Bound...")
  upper_bound <- cor(Matrix_U)
  
  cat(" Upper Bound...\n")
  cor_U_1U <- cor(Matrix_U, Matrix_1U)
  lower_bound[upper.tri(lower_bound)] <- cor_U_1U[upper.tri(cor_U_1U)]
  lower_bound[lower.tri(lower_bound)] <- t(lower_bound)[lower.tri(lower_bound)]
  
  end_time <- Sys.time()
  comp_time <- difftime(end_time, start_time, units = "mins")[[1]]
  
  return(list("UpperBound" = upper_bound, "LowerBound" = lower_bound, "ComputationTime" = comp_time, "Seed" = seed))
}

generate_rho <- function(min_rho, max_rho, n_vars, n_iter = 10, seed = 1234, method = "Gauss"){
  is_valid <- function(R) {
    # Check if all pairwise correlations are in Range and R is positive definite
    return(min_rho <= min(R[lower.tri(R)]) & 
             max_rho >= max(R[lower.tri(R)]) &
             sum(eigen(R)$values < 0) == 0) 
  }
  
  set.seed(seed)
  correct <- 0
  for (i in 1:n_iter) {
    for (j in 1:n_iter) {
      cat("Trying to generate Rho (Method:", method, "), iteration", j, " with correction", correct, "!\r")
      
      # Draw pairwise correlations by chosen Method
      if (method == "Gauss") {
        correlations <- scale2ab(rnorm(sum(1:(n_vars-1))), min_rho + correct, max_rho - correct)
      } else if (method == "Uniform") {
        correlations <- runif(sum(1:(n_vars-1)), min_rho + correct, max_rho - correct)
      } else if (method == "Beta") {
        correlations <- scale2ab(rbeta(sum(1:(n_vars-1)), shape1 = 2, shape2 = 5), min_rho + correct, max_rho - correct)
      } 
      
      # Fill matrix
      R <- matrix(1, n_vars, n_vars)
      R[lower.tri(R)] <- correlations
      R <- t(R)
      R[lower.tri(R)] <- correlations
      
      # Correct to nearest postive definite matrix
      R <- as.matrix(nearPD(R, corr = TRUE, keepDiag = TRUE)$mat)
      
      if (is_valid(R)) {
        cat("Valid Rho found after", j, "iterations with correction", correct, "(Method", method, ")!\n\r")
        return(R)
      }
    }
    # Shrink the intervall by value "correct" if no matrix 
    # is found after n_iter iterations.
    correct <- 0.001 * i
  }
  
  cat("No valid Rho found after", i, "iterations(Method", method, ")!\n")
  return(NULL)
}

correct_rho <- function(rho, rho_lowerbound, rho_upperbound, correction = 0.01) {
  cat("Correcting rho to fit in bounds...\n")
  
  # Correct Target Rho if correlations are not in feasible range.
  for (i in 1:(ncol(rho)-1)) {
    for (j in (i+1):(ncol(rho))) {
      lb <- rho_lowerbound[i, j]
      ub <- rho_upperbound[i, j]
      
      if (rho[i, j] > ub) {rho[i, j] <- rho[j, i] <- ub - correction} 
      if (rho[i, j] < lb) {rho[i, j] <- rho[j, i] <- lb + correction}
    }
  }
  
  if (sum(eigen(rho)$values < 0) > 0) {as.matrix(nearPD(rho, corr = TRUE, keepDiag = TRUE)$mat)}
  return(rho)
}


#===============================================================================
# .------------------------.
# |     SIGMA FUNCTIONS    |
# '------------------------'
#===============================================================================

calculate_sigma <- function(rho_lowerbound, rho_upperbound, rho_target) {
  cat("Calculating Sigma...\n")
  # Calculate Sigma for the Normal Distribution from
  # Yahav and Shmueli (2011): https://doi.org/10.1002/asmb.901
  
  a <- -(rho_upperbound * rho_lowerbound)/(rho_upperbound + rho_lowerbound)
  b <- log((rho_upperbound + a)/a, exp(1))
  c <- -a
  
  sigma <- (1/b) * log((rho_target - c)/a, exp(1))
  diag(sigma) <- 1
  
  sigma[is.na(sigma)] <- 1
  
  if (sum(eigen(sigma)$values < 0) != 0) {sigma <- as.matrix(nearPD(sigma, corr = TRUE, keepDiag = TRUE)$mat)}
  return(sigma)
}

draw_X <- function(n, sigma, seed = 1234){
  # Draw a multivariate Standardnormal matrix with prespecified 
  # correlation sigma. The function extracted from the package
  # SimCorrMix: https://cran.r-project.org/package=SimCorrMix
  
  eig <- eigen(sigma, symmetric = TRUE)
  sqrteigval <- diag(sqrt(pmax(eig$values, 0)))
  eigvec <- eig$vectors
  fry <- eigvec %*% sqrteigval
  set.seed(seed)
  X <- matrix(rnorm(ncol(sigma) * n), n)
  X <- scale(X, TRUE, FALSE)
  X <- X %*% svd(X, nu = 0, nv = ncol(X))$v # Added "nv = ncol(X)" to support n < p
  X <- scale(X, FALSE, TRUE)
  X <- fry %*% t(X)
  X <- t(X)
  return(X)
}


#===============================================================================
# .-------------------------.
# |     HELPER FUNCTIONS    |
# '-------------------------'
#===============================================================================

get_params <- function(n_equal, n_diff, mins, maxs, digits, seed = 1234) {
  n_groups <- length(mins)
  n_vars <- n_equal + n_diff
  params <- list()
  
  # Draw Parameter for each group
  for (i in 1:n_groups) {
    set.seed(seed+i)
    
    if (digits == 0) {
      params_group <- sample(mins[i]:maxs[i], n_vars, replace=TRUE)
    } else {
      params_group <- runif(n_vars, mins[i], maxs[i])
    }
    
    if (i == 1) {
      params[[i]] <- params_group
    } else {
      if (n_equal == 0) { 
        # Only Variables with different Parameters
        params[[i]] <- params_group
      } else if (n_diff == 0) { 
        # Only Variables with same Parameters, take all from group 1
        params[[i]] <- params[[1]]
      } else { 
        # Mixed Variables with same and different Parameters, take n_equal from group 1
        params[[i]] <- c(params[[1]][1:n_equal], params_group[(n_equal+1):n_vars])
      }
    }
  }
  
  # Round Parameter to specified digits
  params <- lapply(params, round, digits)

  return(params)
  # Return list with Parameter Vectors for each group
}
