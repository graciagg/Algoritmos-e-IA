rm(list=ls())

path <- "/Users/vic/Library/CloudStorage/GoogleDrive-vdelaopascual@gmail.com/Mi unidad/MU en Bioinformática (UNIR 2023)/Actividades (AAAA-MM-DD)/Actividad 1_junio2024 update"

setwd(path)

df <- read.csv("Dataset expresión genes.csv")


library(dplyr)
df_genes <- df %>% dplyr::select(starts_with("AQ_"))
str(df_genes)


is.na(colSums(df_genes)) # ver si hay missing
df_genes_scale <- scale(df_genes)  # Normalización z-score



#### ---- Clustering no jerárquico con kmeans ----
library(factoextra)
library(stats)


kmeans.result <- kmeans(df_genes_scale, centers = 2, iter.max = 100, nstart = 25)
fviz_cluster(kmeans.result, df_genes_scale, xlab = '', ylab = '') +
  ggtitle("Cluster plot, centers = 2", subtitle = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10)))


kmeans.result <- kmeans(df_genes_scale, centers = 3, iter.max = 100, nstart = 25)
fviz_cluster(kmeans.result, df_genes_scale, xlab = '', ylab = '') +
  ggtitle("Cluster plot, centers = 3", subtitle = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10)))



# n optimo de clusters
fviz_nbclust(df_genes_scale, kmeans, method = "wss") +
  ggtitle("Optimal number of clusters", subtitle = "") +
  theme_classic()