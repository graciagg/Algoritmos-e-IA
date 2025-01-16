######################
# Lineal L Embedding #
######################

# Directorio en donde tendremos nuestro script ISOMAP.R
setwd("~/Desktop/UNIR/scripts")

# Carga de librerias
library(ica)
library(ggplot2)

set.seed(1234)

# Lectura de datos
data.raw <- read.csv('data/data.csv')
labels.raw<- read.csv('data/labels.csv')

# Guardado en un dataframe de los 500 primeros genes 
data <- sapply(data.raw[2:501], as.numeric)

ica.results <- ica(data, nc=3, method = "fast")

ica.df <- data.frame(ica.results$S)

# Graficamos
ggplot(ica.df, aes(x = X1, y = X2, color = labels.raw$Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método ICA Types of Cancer", x = "X1", y = "X2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))
