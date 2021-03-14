library(SimCorrMix)
library(VGAM)
library(ggplot2)
library(gridExtra)
library(dplyr)

seed <- 30091997
set.seed(seed)

gen_data <- function(n1, n2, p1, p2, size_range, prob_range, p0_range, rho, seed) {
  # Parameter
  n <- n1 + n2
  p <- p1 + p2
  
  # Size
  min_size <- size_range[1]
  max_size <- size_range[2]
  size <- sample(min_size:max_size, p2, replace=TRUE)
  sizes <- list(c(size, sample(min_size:max_size, p1, replace=TRUE)), c(size, sample(min_size:max_size, p1, replace=TRUE)))
  
  # Probability Negativ Binomial
  min_prob <- prob_range[1]
  max_prob <- prob_range[2]
  prob <- runif(p2, min_prob, max_prob)
  probs <- list(c(prob, runif(p1, min_prob, max_prob)), c(prob, runif(p1, min_prob, max_prob)))
  
  #mu <- size * (1 - prob)/prob
  # Probability Zero Inflated Structural Zeros
  min_p0 <- p0_range[1]
  max_p0 <- p0_range[2]
  p_zinb <- runif(p2, min_p0, max_p0)
  ps_zinb <- list(c(p_zinb, runif(p1, min_p0, max_p0)), c(p_zinb, runif(p1, min_p0, max_p0)))
  
  # Correlation Matrix
  rho_matrix <- matrix(rho, p, p)
  diag(rho_matrix) <- 1
  
  # Check input / Estimate computation time
  validpar(k_nb = p, method = "Polynomial", size = sizes[[1]], prob = probs[[1]], p_zinb = ps_zinb[[1]], rho = rho_matrix)
  validpar(k_nb = p, method = "Polynomial", size = sizes[[2]], prob = probs[[2]], p_zinb = ps_zinb[[2]], rho = rho_matrix)
  validcorr(n = n, k_nb = p, method = "Polynomial", size = sizes[[1]], prob = probs[[1]], p_zinb = ps_zinb[[1]], rho = rho_matrix)
  validcorr(n = n, k_nb = p, method = "Polynomial", size = sizes[[2]], prob = probs[[2]], p_zinb = ps_zinb[[2]], rho = rho_matrix)
  
  # Generate data
  data1 <- corrvar(n = n1, k_nb = p, method = "Polynomial", size = sizes[[1]], prob = probs[[1]], p_zinb = ps_zinb[[1]], rho = rho_matrix, seed = seed)
  data2 <- corrvar(n = n2, k_nb = p, method = "Polynomial", size = sizes[[2]], prob = probs[[2]], p_zinb = ps_zinb[[2]], rho = rho_matrix, seed = seed)
  
  data <- data.frame(do.call(rbind, list(data1$Y_nb, data2$Y_nb)))

  return(list(data = data, data1 = data1, data2 = data2))
}

get_time_estimate <- function(p) {
  # Parameters of a polynomial model, estimated by previous data generation runs
  t <- 0.548182 + p*-0.009764 + (p^2)*0.004375
  cat("Estimated computation time for the dataset: ", round(t*2/60, 2), " hours...", sep = "") 
}

rhos <- c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9)
sizes <- list(c(45, 293), c(27, 397), c(19, 576))
probs <- list(c(0.27, 0.47), c(0.24, 0.55), c(0.18, 0.78))
p0s <- list(c(5.30*10^-7, 0.01), c(3.65*10^-7, 0.04), c(2.28*10^-7, 0.08))
n1 <- n2 <- 500
p1 <- 100
p2 <- 0

for (rho in rhos) {
  get_time_estimate(p1+p2)
  
  ret <- gen_data(n1, n2, p1, p2, c(sizes[[1]][1], sizes[[1]][2]), c(probs[[1]][1], probs[[1]][2]), c(p0s[[1]][1], p0s[[1]][2]), rho, seed)
  #data1 <- gen_data(n1, n2, p1, p2, c(sizes[[i]][1], sizes[[i]][2]), c(probs[[i]][1], probs[[i]][2]), c(p0s[[i]][1], p0s[[i]][2]), 1)
  #data2 <- gen_data(n1, n2, p1, p2, c(sizes[[i]][1], sizes[[i]][2]), c(probs[[i]][1], probs[[i]][2]), c(p0s[[i]][1], p0s[[i]][2]), 2)
  #data <- data.frame(do.call(rbind, list(data1$Y_nb, data2$Y_nb)))
  
  # Save whole dataset and subset objects
  name <- paste0("r",round(rho*100,0),"_n",n1+n2,"_p",p1+p2,"_s",sizes[[i]][1],"-",sizes[[i]][2],"_p",round(probs[[i]][1]*100,0),"-",round(probs[[i]][2]*100,0))
  saveRDS(ret$data, file=paste0("data/datasets/",sprintf("%02d", i),"_",name,".Rda"))
  saveRDS(ret$data1, file=paste0("data/datasets/",sprintf("%02d", i),"_d1_",name,".Rds"))
  saveRDS(ret$data2, file=paste0("data/datasets/",sprintf("%02d", i),"_d2_",name,".Rds"))
}

# Save dataset
saveRDS(data, file="datasets/n500_p250_s10-150_p02-06.Rda")
