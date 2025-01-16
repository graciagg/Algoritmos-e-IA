######################
# Lineal L Embedding #
######################

# Directorio en donde tendremos nuestro script ISOMAP.R
setwd("~/Desktop/UNIR/scripts")

# Carga de librerias
library(RDRToolbox)
library(lle)
library(ggplot2)

# Lectura de datos
data.raw <- read.csv('data/data.csv')
labels.raw<- read.csv('data/labels.csv')

# Guardado en un dataframe de los 500 primeros genes 
data <- sapply(data.raw[2:501], as.numeric)

# Algoritmo
# Funcion LLE()
#   x: matriz sobre la cual se va a reducir la dimensionalidad
#   dim: numero de dimensiones de salida
#   k: numero de vecinos cercanos. Puede aproximarse con la funcion calc_k()

#   la salida será un dataframe con dimensión = dim

# Funcion calc_k()
#   x: matriz sobre la cual se va a reducir la dimensionalidad
#   k_min, k_max: intervalo para la busqueda de la k optima
#   plotres: representará los resultados de cada k
#   paralel: utilización del calculo paralelo del ordenador
#   cpus: numero de procesadores que se quiera utilizar (max nº procesadores maquina)


lle.results <- LLE(data, dim = 2, k=130)

lle.df <- data.frame(lle.results)

# Graficamos
ggplot(lle.df, aes(x = X1, y = X2, color = labels.raw$Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método LLE Types of Cancer", x = "PC1", y = "PC2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))

