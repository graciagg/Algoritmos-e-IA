#####################
# Isometric Mapping #
#####################

# Directorio en donde tendremos nuestro script ISOMAP.R
setwd("~/Desktop/UNIR/scripts")

# Carga de librerias
# En caso de no tener descargado RDRToolbox


library(RDRToolbox)
library(ggplot2)

# Lectura de datos
data.raw <- read.csv('data/data.csv')
labels.raw<- read.csv('data/labels.csv')

# Guardado en un dataframe de los 500 primeros genes 
data <- sapply(data.raw[2:501], as.numeric) #no ponemos dataframe!

# Algoritmo

# Funcion Isomap()
#   data -> datos (matriz) sobre los que haremos reduccion de dimensionalidad
#   dim -> dimensiones de las columnas del espacio reducido
#   k -> numero de vecinos cercanos a cada punto. A mayor k mayor computacion
#   potResiduals -> devuelve la varianza explicada por las diferentes dimensiones

#   Si se ha elegido una unica dimension devuelve una matriz
#   Si se ha elegido un vector de dimensiones devolvera una matriz por cada elemento del vector

# Calculamos isomap de 1 a 10 dimensiones y con 5 vecinos (y 30?)
isomap.results = Isomap(data=data, dims=1:10, k=5, plotResiduals=TRUE)

# Dataframe con los puntos que queremos dibujar en el plano 2D
#     (elegiriamos otro si queremos otra dimension)
isomap.df <- data.frame(isomap.results$dim2) 

# Graficamos
ggplot(isomap.df, aes(x = X1, y = X2, color = labels.raw$Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Isomap Types of Cancer", x = 'dim1', y = 'dim2', color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))
