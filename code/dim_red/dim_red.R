library(uwot)
library(ggplot2)

data <- readRDS(file="datasets/n500_p25_s10-150_p02-06.Rda")

# UMAP
umap_raw <- umap(data)
umap_projection <- data.frame(x = umap_raw[,1], 
                              y= umap_raw[,2], 
                              index = factor(rep(c(1, 2), each = 250)))
ggplot(umap_projection, aes(x = x, y = y, col = index)) +
  geom_point()

saveRDS(umap_projection, file="dim_red/01_umap.Rda")

df <- data.frame(x = dim_red[,1], y= dim_red[,2], index = factor(rep(c(1, 2), each = 250)))
ggplot(df, aes(x = x, y = y, col = index)) +
  geom_point()


dim_red <- prcomp(data, center = TRUE, scale. = TRUE)
df <- data.frame(PC1 = dim_red$x[,1], PC2 = dim_red$x[,2], index = factor(rep(c(1, 2), each = 250)))
ggplot(df, aes(x = PC1, y = PC2, col = index)) +
  geom_point()

ggplot(data.frame(x = data$X1, y = data$X2, index = factor(rep(c(1, 2), each = 250))), aes(x = x, y = y, col = index)) +
  geom_point()
