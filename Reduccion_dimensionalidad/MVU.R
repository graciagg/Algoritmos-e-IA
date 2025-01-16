######################
# Lineal L Embedding #
######################

# Directorio en donde tendremos nuestro script ISOMAP.R
setwd("~/Desktop/UNIR/scripts")

# Carga de librerias
library(Rdimtools)
library(ggplot2)

# Lectura de datos
data.raw <- read.csv('data/data.csv')
labels.raw<- read.csv('data/labels.csv')

# Guardado en un dataframe de los 500 primeros genes 
data <- sapply(data.raw[1:100, 2:501], as.numeric)

# Definimos k a través del argumento type de la funcion do.mvu
# Funcion do.mvu()
#   x: matriz donde reduciremos la dimensionalidad
#   ndim: dimension final de los datos
#   type: cantidad de puntos vecinos
#   preprocess: preprocesamiento de datos

# Variable Y

# A más muestras más tiempo de ejecución


mvu.results <- do.mvu(data, ndim=2, type=c("proportion", 0.1))

mvu.df <- data.frame(mvu.results$Y)

# Graficamos
ggplot(mvu.df, aes(x = X1, y = X2, color = labels.raw[1:100,]$Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método MVU Types of Cancer", x = "X1", y = "X2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))
