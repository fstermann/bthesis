library(SimCorrMix)
library(VGAM)
library(ggplot2)

n <- 500
p <- 10
size <- sample(1:19, p, replace=TRUE)
prob <- runif(p, 0.2, 0.8)
mu <- size * (1 - prob)/prob
p_zinb <- runif(p, 0.15, 0.25)

k_nb <- length(size)

Rey <- matrix(0.35, k_nb, k_nb)
diag(Rey) <- 1

rownames(Rey) <- colnames(Rey) <- sprintf("NB%s",seq(1:k_nb))

# Check input
validpar(k_nb = k_nb, method = "Polynomial", size = size, prob = prob, mu = NULL, p_zinb = p_zinb, rho = Rey)

data <- corrvar(n = n, k_nb = k_nb, method = "Polynomial", size = size, prob = prob, mu = NULL, p_zinb = p_zinb, rho = Rey)
data2 <- corrvar2(n = n, k_nb = k_nb, method = "Polynomial", size = size, prob = prob, mu = NULL, p_zinb = p_zinb, rho = Rey)

head(data$Y_nb)
summary(data$Y_nb)
data2$Sigma

plot_simpdf_theory(sim_y = data2$Y_nb[,3], Dist = "Negative_Binomial",
                   params = c(p_zinb), cont_var = FALSE, col_width = 0.25)

plot_simtheory(sim_y = data2$Y_nb[,3], Dist = "Negative_Binomial",
                   params = c(p_zinb), cont_var = FALSE, binwidth = 0.25)

rzinegbin(n, size = size[3], munb = mu[3], pstr0 = p_zinb[3])

compare_results <- function(data, i){
  y_sim <- data$Y_nb[,i]
  y_true <- rzinegbin(n, size = size[i], munb = mu[i], pstr0 = p_zinb[i])
  d <- data.frame(append(y_sim, y_true), rep(c("Sim", "True"), each = length(y_sim)))
  colnames(d) <- c("Values", "Fill")
  
  ggplot(d, aes(x = Values, fill = Fill)) +
    #geom_bar(stat="count", position=position_dodge(), width = 0.6) +
    geom_density(aes(alpha = 0.5))
}

compare_results(data2, 4)
