############################
# Multidimensional Scaling #
#############################

# Directorio en donde tendremos nuestro script MSD.R
setwd("C:/Users/graci/Documents/BIOINFORMATICA_UNIR/Algoritmos e Inteligencia Artificial/Tema 1 Ejercicios/Breast Cancer")

# Carga de librerias
library(stats)   # librería para el MDS
library(ggplot2) # librería para hacer la representación gráfica

# Lectura de datos
data.raw <- read.csv('data.csv')
labels.raw <- read.csv('labels.csv')

# Guardado en un dataframe de los 500 primeros genes 
data <- data.frame(sapply(data.raw[2:501], as.numeric))

# Algoritmo #

# Funcion cmdscale()
#   d: matriz de distancias (usaremos la funcion dist)
#   k: numero que indica el tamaño final de los datos (max num de variables)
#   eig: si calculamos autovalores de las variables. Nos sirve para el calculo 
#        de la varianza explicada, es decir, para coger las columnas de mayor
#        variabilidad
#   x.ret: para devolver la matriz de distancias que calcule el algoritmo

#   points: dataframe de tamaño k que representa las nuevas coordenadas
#   eig: vector con los autovalores para elegir el numero de dimensiones


# Utilizamos la funcion dist para calcular la matriz de distancias euclideas
# matriz NxN de distancias entre todos los puntos
distances <- dist(data, method = 'euclidean')

# Utilizamos la funcon cmdscale para realizar el MSD
mds.results <- cmdscale(distances, eig=TRUE, k=2, x.ret=TRUE)

# Calculamos la varianza explicada
varianza.explicada <- mds.results$eig/sum(mds.results$eig) * 100

# Sacamos en un dataframe los puntos del mds
mds.df <- data.frame(mds.results$points)


# Grafico
ggplot(mds.df, aes(x=X1, y=X2, color=labels.raw$Class)) +
  geom_point(size=3) + 
  scale_color_manual(values=c("red", "blue", "green", "orange", "purple")) +
  labs(title="MDS Different types of Cancer", x="Dimension 1 (X1)", y="Dimension 2 (X2)", color = "Grupo") +
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))
