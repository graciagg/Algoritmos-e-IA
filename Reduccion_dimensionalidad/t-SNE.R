###############################################
# T-distributed Stochastic Neighbor Embedding #
###############################################

# Directorio en donde tendremos nuestro script t-SNE.R
setwd("~/Desktop/UNIR/scripts")

# Carga de librerias
library(ggplot2)
library(Rtsne)

# Seteamos la semilla para que sea replicable el algoritmo
set.seed(1234)

# Lectura de datos
data.raw <- read.csv('data/data.csv')
labels.raw <- read.csv('data/labels.csv')

# Guardado en un dataframe de los 500 primeros genes 
data <- sapply(data.raw[2:501], as.numeric)

# Algoritmo
# funcion Rtsne()
#   X: datos sobre los que reduciremos la dimensionalidad
#   dims: tamaño final del conjunto de datos (mejor <=3) por eficiencia
#   num_threads: hilos a utilizar (procesadores). No hace falta usarlo de momento
# 
#   Variable Y con matriz del t-SNE

tsne <- Rtsne(X=data)
tsne_result <- data.frame(tsne$Y)

# Graficamos
ggplot(tsne_result, aes(x = X1, y = X2, color = labels.raw$Class)) +
geom_point(size = 3) +
scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
labs(title = "Método t-SNE", x = "PC1", y = "PC2", color = "Grupo") +
theme_classic() +
theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))
