library(VGAM)
library(ggplot2)
library(gridExtra)
library(fitdistrplus)

data <- readRDS(file="datasets/n500_p25_s10-150_p02-06.Rda")
  
# Compare plots
get_abs_value <- function(x, y){
  is.integer0 <- function(x){is.integer(x) && length(x) == 0L}
  
  if (is.integer0(x)){x <- 0}
  if (is.integer0(y)){y <- 0}
  
  return(abs(x-y))
}

dens <- density(data[,1])
dens[1:2]

xtemp <- sort(data$X1)
dtemp <- data.frame(x = data$X1, y = density(data$X1, n = 500, from = 0, to = max(data$X1))$y, group = rep(c(1,2), each = 250))
ggplot(dtemp, aes(x = x)) +
  geom_histogram(aes(y = ..density..), binwidth = 1) +
  geom_line(aes(y = y, x = density(data$X1, n = 500, from = 0, to = max(data$X1))$x)) +
  facet_grid(vars(group))

#ggplot(data, aes(x = factor(X1))) +
#  geom_bar(stat="count", position=position_dodge(), width = 0.6) + xlab("") +
#  scale_x_discrete(breaks = seq(0, max(data$X1), by = 5))

dzinegbin(c(2, 2, 2), size = sizes[[1]][i], 
          munb = sizes[[1]][i] * (1 - probs[[1]][i])/probs[[1]][i], 
          pstr0 = ps_zinb[[1]][i])

seq(0, length.out=250, by=ceiling(400/250))

compare_results <- function(data, i){
  min_x <- -1
  d <- data[,i]
  d1 <- d[1:n1]
  d2 <- d[(n1+1):n]
  
  y1_sim <- density(d1, n=250, from=min_x, to=max(d1))$y
  x1_sim <- density(d1, n=250, from=min_x, to=max(d1))$x
  y2_sim <- density(d2, n=250, from=min_x, to=max(d2))$y
  x2_sim <- density(d2, n=250, from=min_x, to=max(d2))$x
  
  x1_true <- seq(min_x, length.out=250, by=ceiling(max(d1)/250))
  y1_true <- dzinegbin(x1_true, size = sizes[[1]][i], 
                       munb = sizes[[1]][i] * (1 - probs[[1]][i])/probs[[1]][i], 
                       pstr0 = ps_zinb[[1]][i])
  
  x2_true <- seq(min_x, length.out=250, by=ceiling(max(d2)/250))
  y2_true <- dzinegbin(x2_true, size = sizes[[2]][i], 
                       munb = sizes[[2]][i] * (1 - probs[[2]][i])/probs[[2]][i], 
                       pstr0 = ps_zinb[[2]][i])
  
  df <- data.frame(d = c(d1, d2), 
                   x_sim = c(x1_sim, x2_sim),
                   y_sim = c(y1_sim, y2_sim),
                   x_true = c(x1_true, x2_true),
                   y_true = c(y1_true, y2_true),
                   group = c(rep(1, length(d1)), rep(2, length(d2))))
  #colnames(df) <- c("d", "x_sim", "y_sim", "x_true", "y_true", "group")
  
  df_text <- data.frame(x = c(Inf, Inf), #c(max(d)*0.8, max(d)*0.8),
                        y = c(Inf, Inf), #c(max(df$y_true)*0.8, max(df$y_true)*0.8),
                        group = c(1, 2),
                        label = c(paste0("size = ", sizes[[1]][i], "\n", "prob = ", round(probs[[1]][i], 2)), 
                                  paste0("size = ", sizes[[2]][i], "\n", "prob = ",  round(probs[[2]][i], 2))))
  # 
  # df <- data.frame(c(x1_sim, x2_sim, 0:(length(y1_true)-1), 0:(length(y2_true)-1)), 
  #                 c(y1_sim, y2_sim, y1_true, y2_true), 
  #                 c(rep("Sim1", length(y1_sim)), rep("Sim2", length(y2_sim)), rep("True1", length(y1_true)), rep("True2", length(y2_true))),
  #                 c(rep(1, length(y1_sim)), rep(2, length(y2_sim)), rep(1, length(y1_true)), rep(2, length(y2_true))))
  # colnames(df) <- c("x", "y", "fill", "group")
  
  # Plots
  #p_bars <- ggplot(d, aes(x = factor(Values), fill = Fill)) +
  #  geom_bar(stat="count", position=position_dodge(), width = 0.6) + xlab("")
  p_density <- ggplot(df, aes(x = d)) +
    geom_histogram(aes(y = ..density..), binwidth = 1) +
    #geom_area(stat = "identity", alpha = 0.5) +
    geom_line(aes(x = x_true, y = y_true), col = "green") +
    geom_area(aes(x = x_true, y = y_true), fill = "green", alpha = 0.3) +
    geom_line(aes(x = x_sim, y = y_sim), col = "yellow") +
    geom_area(aes(x = x_sim, y = y_sim), fill = "yellow", alpha = 0.3) +
    xlab("") + xlim(c(min_x, max(d))) +
    geom_text(data = df_text, aes(x = Inf, y = Inf, label = label), hjust = "inward", vjust = "inward") +
    facet_grid(vars(group))
  #p_all <- grid.arrange(p_bars, p_density, ncol=2)
  
  # Loss
  #l_tbl <- d %>% group_by(Fill, Values) %>% summarise(counts = n()) %>% data.frame()
  #loss <- 0
  #l_n <- max(l_tbl$Values)
  #for (j in 1:l_n){
  #  loss <- loss + get_abs_value(l_tbl[l_tbl["Fill"] == "Sim" & l_tbl["Values"] == j,]$counts, 
  #                               l_tbl[l_tbl["Fill"] == "True" & l_tbl["Values"] == j,]$counts)
  #}
  
  return(list(plot = p_density, 
              loss = 1)) #loss/l_n))
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
results <- compare_results(data, 21)
results$loss
plot(results$plot)

# Losses over all columns
losses <- get_avg_loss(data)
losses$AvgLoss
losses$SingleLosses
