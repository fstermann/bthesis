library(SimCorrMix)
library(VGAM)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Parameters
n <- 500
n1 <- round(n * 0.5, 0)
n2 <- n-n1
p <- 25
p1 <- round(p * 0.2, 0)
p2 <- p-p1

min_size <- 10
max_size <- 150
size <- sample(min_size:max_size, p2, replace=TRUE)
sizes <- list(c(size, sample(min_size:max_size, p1, replace=TRUE)), c(size, sample(min_size:max_size, p1, replace=TRUE)))

prob <- runif(p2, 0.2, 0.6)
probs <- list(c(prob, runif(p1, 0.2, 0.6)), c(prob, runif(p1, 0.2, 0.6)))

#mu <- size * (1 - prob)/prob
p_zinb <- runif(p2, 0, 0.1)
ps_zinb <- list(c(p_zinb, runif(p1, 0, 0.1)), c(p_zinb, runif(p1, 0, 0.1)))

k_nb <- p

# Covariance
rho <- matrix(0.35, p, p)
diag(rho) <- 1
rownames(rho) <- colnames(rho) <- sprintf("ZINB%s", seq(1:p))

length(sizes[[1]])
length(probs[[1]])

# Check input
validpar(k_nb = k_nb, method = "Polynomial", size = sizes[[1]], prob = probs[[1]], p_zinb = ps_zinb[[1]], rho = rho)
validpar(k_nb = k_nb, method = "Polynomial", size = sizes[[2]], prob = probs[[2]], p_zinb = ps_zinb[[2]], rho = rho)

# Generate data
data1 <- corrvar(n = n1, k_nb = k_nb, method = "Polynomial", size = sizes[[1]], prob = probs[[1]], p_zinb = ps_zinb[[1]], rho = rho)
data2 <- corrvar(n = n2, k_nb = k_nb, method = "Polynomial", size = sizes[[2]], prob = probs[[2]], p_zinb = ps_zinb[[2]], rho = rho)

data <- data.frame(do.call(rbind, list(data1$Y_nb, data2$Y_nb)))

# Save dataset
saveRDS(data, file="datasets/n500_p25_s10-150_p02-06.Rda")
