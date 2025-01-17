#Carga de datos y establecimiento de pwd#
rm(list=ls())

path <- "/Users/vic/Library/CloudStorage/GoogleDrive-vdelaopascual@gmail.com/Mi unidad/MU en Bioinformática (UNIR 2023)/Actividades (AAAA-MM-DD)/Actividad 1_junio2024 update"

setwd(path)

df <- read.csv("Dataset expresión genes.csv")


library(dplyr)
df_genes <- df %>% dplyr::select(starts_with("AQ_"))
str(df_genes)


is.na(colSums(df_genes)) # ver si hay missing
df_genes_scale <- scale(df_genes)  # Normalización z-score


#### ---- Clustering jerarquico divisivo (de arriba a abajo): datos con muchas observaciones ----
library(ggplot2) 
library(factoextra)
library(cluster)

# Implementación del clustering divisivo
diana_euclidean <- diana(df_genes_scale, metric = "euclidean", stand = FALSE) # ideal para datos donde las distancias más pequeñas
diana_manhattan <- diana(df_genes_scale, metric = "manhattan", stand = FALSE) # ideal para datos con diferentes escalas o datos categóricos

#Graficamos ambos tipos de clustering
colors <- rainbow(5)
clust_diana_euclidean <- fviz_dend(diana_euclidean, 
                                   cex = 0.5, 
                                   k = 5,
                                   palette = colors, 
                                   main = 'Euclidean',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()


colors <- rainbow(5)
clust_diana_manhattan <- fviz_dend(diana_manhattan, 
                                   cex = 0.5, 
                                   k = 5,
                                   palette = colors, 
                                   main = 'Manhattan',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()


grid.arrange(clust_diana_euclidean, clust_diana_manhattan, nrow = 2)









rm(list=ls())

# Cargar librerías
library(ggplot2)
library(cluster)

# Ejemplo de coordenadas (representación de plantas, clientes, tiendas, etc.)
data <- data.frame(
  x = c(2, 4, 6, 10),  # Coordenadas x (pueden ser precios, calificaciones, etc.)
  y = c(3, 4, 7, 10)   # Coordenadas y
)

# Calcular la matriz de distancias (usando Euclidiana)
dist_matrix <- dist(data)
dist_matrix

# Realizar el clustering con diferentes métodos de enlace
hclust_single <- hclust(dist_matrix, method = "single")
hclust_complete <- hclust(dist_matrix, method = "complete")
hclust_average <- hclust(dist_matrix, method = "average")
hclust_ward <- hclust(dist_matrix, method = "ward.D2")

hclust_single$order
hclust_complete$order
hclust_average$order
hclust_ward$order

# Graficar resultados
par(mfrow = c(2, 2))  # Dividir la ventana en 2x2 para graficar múltiples
plot(hclust_single, main = "Single Linkage")
plot(hclust_complete, main = "Complete Linkage")
plot(hclust_average, main = "Average Linkage")
plot(hclust_ward, main = "Ward's Linkage")
