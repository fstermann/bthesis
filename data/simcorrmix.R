library(SimCorrMix)
library(VGAM)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Parameter
n <- 500
p <- 10
size <- sample(1:19, p, replace=TRUE)
prob <- runif(p, 0.2, 0.8)
mu <- size * (1 - prob)/prob
p_zinb <- runif(p, 0.15, 0.25)
k_nb <- p

# Covariance
rho <- matrix(0.35, p, p)
diag(rho) <- 1
rownames(rho) <- colnames(rho) <- sprintf("ZINB%s", seq(1:p))

# Check input
validpar(k_nb = k_nb, method = "Polynomial", size = size, mu = mu, p_zinb = p_zinb, rho = rho)

# Generate data
data <- corrvar(n = n, k_nb = k_nb, method = "Polynomial", size = size, mu = mu, p_zinb = p_zinb, rho = rho)
data$Sigma

# Compare plots
get_abs_value <- function(x, y){
  is.integer0 <- function(x){is.integer(x) && length(x) == 0L}
  
  if (is.integer0(x)){x <- 0}
  if (is.integer0(y)){y <- 0}
  
  return(abs(x-y))
}

compare_results <- function(data, i){
  y_sim <- data$Y_nb[,i]
  y_true <- rzinegbin(length(y_sim), size = size[i], munb = mu[i], pstr0 = p_zinb[i])
  d <- data.frame(append(y_sim, y_true), rep(c("Sim", "True"), each = length(y_sim)))
  colnames(d) <- c("Values", "Fill")
  
  # Plots
  p_bars <- ggplot(d, aes(x = factor(Values), fill = Fill)) +
    geom_bar(stat="count", position=position_dodge(), width = 0.6) + xlab("")
  p_density <- ggplot(d, aes(x = Values, fill = Fill)) +
    geom_density(aes(alpha = 0.5)) + xlab("")
  p_all <- grid.arrange(p_bars, p_density, ncol=2)
  
  # Loss
  l_tbl <- d %>% group_by(Fill, Values) %>% summarise(counts = n()) %>% data.frame()
  loss <- 0
  l_n <- max(l_tbl$Values)
  for (j in 1:l_n){
    loss <- loss + get_abs_value(l_tbl[l_tbl["Fill"] == "Sim" & l_tbl["Values"] == j,]$counts, 
                                 l_tbl[l_tbl["Fill"] == "True" & l_tbl["Values"] == j,]$counts)
  }

  return(list(plot = p_all, 
              loss = loss/l_n))
}

get_avg_loss <- function(data){
  avg_loss <- 0
  sgl_loss <- rep(0, ncol(data$Y_nb))
  
  for (p in 1:ncol(data$Y_nb)){
    loss <- 0
    y_sim <- data$Y_nb[,p]
    y_true <- rzinegbin(length(y_sim), size = size[p], munb = mu[p], pstr0 = p_zinb[p])
    d <- data.frame(append(y_sim, y_true), rep(c("Sim", "True"), each = length(y_sim)))
    colnames(d) <- c("Values", "Fill")
    
    l_tbl <- d %>% group_by(Fill, Values) %>% summarise(counts = n()) %>% data.frame()
    l_n <- max(l_tbl$Values)
    for (j in 1:l_n){
      loss <- loss + get_abs_value(l_tbl[l_tbl["Fill"] == "Sim" & l_tbl["Values"] == j,]$counts, 
                                   l_tbl[l_tbl["Fill"] == "True" & l_tbl["Values"] == j,]$counts)
    }
    sgl_loss[p] <- loss/l_n
    avg_loss <- avg_loss + loss/l_n
  }
  
  return(list(AvgLoss = avg_loss/ncol(data$Y_nb), 
              SingleLosses = sgl_loss))
}

# Plot & Loss
results <- compare_results(data, 1)
results$loss
plot(results$plot)

# Losses over all columns
losses <- get_avg_loss(data)
losses$AvgLoss
losses$SingleLosses
