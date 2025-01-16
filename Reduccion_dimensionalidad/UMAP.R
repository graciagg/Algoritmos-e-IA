######################
# Lineal L Embedding #
######################

# Directorio en donde tendremos nuestro script ISOMAP.R
setwd("~/Desktop/UNIR/scripts")

# Carga de librerias
library(uwot)
library(ggplot2)

# Lectura de datos
data.raw <- read.csv('data/data.csv')
labels.raw<- read.csv('data/labels.csv')

# Guardado en un dataframe de los 500 primeros genes 
data <- sapply(data.raw[2:501], as.numeric)

# Se dejan los parametros a excepcion del n_neighbours

# Funcion umap()
#     x: reduccion de la dimensionalidad
#     n_neighbours: entero qeu indica el numero de vecinos cercanos
#     n_componentes: entero que determina el tamaño del espacio de salida
#     metric: define la distancia entre puntos
#     min_dist: distancia minima permitida entre puntos
#     scale: tipos de escalado
#     verbose: tiempo hasta que se complete el calculo

#     Y resultado




umap.results <- umap(data, n_neighbors=0.2 * nrow(data),
                     n_components = 2, min_dist = 0.01, local_connectivity=1, ret_model = TRUE)

umap.df <- data.frame(umap.results$embedding)

# Graficamos
ggplot(umap.df, aes(x = X1, y = X2, color = labels.raw$Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método UMAP Types of Cancer", x = "X1", y = "X2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))
