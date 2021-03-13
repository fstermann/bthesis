library(gridExtra)
library(ggplot2)

data <- readRDS(file="dim_red/01_umap.Rda")

# kmeans
umap_kmeans <- kmeans(data[,c("x", "y")], centers = 2)
data$kmeans <- factor(umap_kmeans$cluster)

p_cluster <- ggplot(data, aes(x = x, y = y, col = kmeans)) +
  geom_point()
p_true <- ggplot(data, aes(x = x, y = y, col = index)) +
  geom_point()
p_all <- grid.arrange(p_cluster, p_true, ncol=2)
p_all
