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
data <- sapply(data.raw[2:501], as.numeric)

# Definimos los parametros k vecinos mas cercanos como el 20% de la muestra
# Funcion do.lapeig()
#   X: matriz sobre la que se reducirá la dimensionalidad
#   Ndim: dimensiones de los datos finales
#   type: forma de generar vecinos
#   preprocess: 6 tipos de preprocesamiento
#   weighted: si queremos darle pesos

#   Y: nuevo espacio con menor dimensionalidad
#   eigenvals: vector con los valores propios de las columnas de Y



le.results <- do.lapeig(data, type=c("proportion", 0.2), weighted=FALSE)

le.df <- data.frame(le.results$Y)

# Graficamos
ggplot(le.df, aes(x = X1, y = X2, color = labels.raw$Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método LE Types of Cancer", x = "X1", y = "X2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))
