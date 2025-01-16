rm(list=ls())

path <- "/Users/vic/Library/CloudStorage/GoogleDrive-vdelaopascual@gmail.com/Mi unidad/MU en Bioinformática (UNIR 2023)/Actividades (AAAA-MM-DD)/Actividad 1_junio2024 update"

setwd(path)

df <- read.csv("Dataset expresión genes.csv")




#### ---- X = GENES Y= TUMOR; training testing ----
library(dplyr)
library(caret)

df_genes_tumor <- df %>% dplyr::select(starts_with("AQ_"), tumor) # lo pongo dplyr::select para que lo detecte dplyr ya que luego con MASS tiene la misma función
df_genes_tumor$tumor <- as.factor(df_genes_tumor$tumor)

set.seed(1995)
train_index <- createDataPartition(df_genes_tumor$tumor, p = 0.8, list = FALSE)
train_index

training_tumor <- df_genes_tumor[train_index, ]
testing_tumor <- df_genes_tumor[-train_index, ]

table(training_tumor$tumor)
table(testing_tumor$tumor)

str(training_tumor)


numerical_columns <- training_tumor[, sapply(training_tumor, is.numeric)]
scaled_data <- scale(numerical_columns)
training_tumor <- cbind(scaled_data, tumor = training_tumor$tumor)

numerical_columns <- testing_tumor[, sapply(testing_tumor, is.numeric)]
scaled_data <- scale(numerical_columns)
testing_tumor <- cbind(scaled_data, tumor = testing_tumor$tumor)

training_tumor <- as.data.frame(training_tumor)
testing_tumor <- as.data.frame(testing_tumor)

# Crear la formula sumando cada gen
genes <- colnames(training_tumor[1:46])
formula_tumor <- as.formula(paste("tumor ~", paste(genes, collapse = "+")))
formula_tumor
training_tumor$tumor <- as.factor(training_tumor$tumor)



#### ---- X = GENES Y= TRAT; training testing ----

df_genes_trat <- df %>% dplyr::select(starts_with("AQ_"), trat) # lo pongo dplyr::select para que lo detecte dplyr ya que luego con MASS tiene la misma función
df_genes_trat$trat <- as.factor(df_genes_trat$trat)

set.seed(1995)
train_index <- createDataPartition(df_genes_trat$trat, p = 0.8, list = FALSE)
train_index

training_trat <- df_genes_trat[train_index, ]
testing_trat <- df_genes_trat[-train_index, ]

table(training_trat$trat)
table(testing_trat$trat)

str(training_trat)


numerical_columns <- training_trat[, sapply(training_trat, is.numeric)]
scaled_data <- scale(numerical_columns)
training_trat <- cbind(scaled_data, trat = training_trat$trat)

numerical_columns <- testing_trat[, sapply(testing_trat, is.numeric)]
scaled_data <- scale(numerical_columns)
testing_trat <- cbind(scaled_data, trat = testing_trat$trat)

training_trat <- as.data.frame(training_trat)
testing_trat <- as.data.frame(testing_trat)

# Crear la formula sumando cada gen
genes <- colnames(training_trat[1:10])
formula_trat <- as.formula(paste("trat ~", paste(genes, collapse = "+")))
formula_trat
training_trat$trat <- as.factor(training_trat$trat)




#### ---- FDA (para datos no paramétricos, tarda mucho en correr...) ----
library(mda)

# Ajustar el modelo RDA en el entrenamiento
fda_model <- fda(formula_tumor, data = training_tumor, method = bruto)
fda_model

fda_pred <- predict(fda_model, newdata = training_tumor)


# Realizar predicciones sobre el conjunto de prueba
fda_predictions <- predict(fda_model, newdata = testing_tumor)

# Obtener la predicción (predicciones de la clase)
predicted_classes <- rda_predictions$class
predicted_classes
length(predicted_classes)

# Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(testing_tumor$tumor)
true_classes
length(true_classes)


# Crear la matriz de confusión (trat predicho testing vs. trat real testing)
confusion <- confusionMatrix(predicted_classes, true_classes)
print(confusion)



























# simulamos más datos (no importa para los bioinfos)
library(synthpop)
library(dplyr)

## simulacion nº1
set.seed(1995) 
df_sintetico <- syn(df)  # Crea un objeto con los datos sintéticos
df_sintetico <- df_sintetico$syn # Obtener el dataframe de pacientes simulados

# Verifica los datos
str(df_sintetico)
str(df)
table(df_sintetico$tumor)
max(df$id)

df_sintetico$id <- 69:(69 + nrow(df_sintetico) - 1)
all(names(df) == names(df_sintetico))  # asegurar de que ambas bases tengan las mismas columnas = Esto debe devolver TRUE
df <- bind_rows(df, df_sintetico) # combinamos
any(duplicated(df$id))  # duplicados en id = debe ser FALSE



## simulacion nº2
set.seed(1995) 
df_sintetico <- syn(df)  # Crea un objeto con los datos sintéticos
df_sintetico <- df_sintetico$syn # Obtener el dataframe de pacientes simulados

# Verifica los datos
str(df_sintetico)
str(df)
table(df_sintetico$tumor)
max(df$id)

df_sintetico$id <- 134:(134 + nrow(df_sintetico) - 1)
all(names(df) == names(df_sintetico))  # asegurar de que ambas bases tengan las mismas columnas = Esto debe devolver TRUE
df <- bind_rows(df, df_sintetico) # combinamos
any(duplicated(df$id))  # duplicados en id = debe ser FALSE


## simulacion nº3
set.seed(1995) 
df_sintetico <- syn(df)  # Crea un objeto con los datos sintéticos
df_sintetico <- df_sintetico$syn # Obtener el dataframe de pacientes simulados

# Verifica los datos
str(df_sintetico)
str(df)
table(df_sintetico$tumor)
max(df$id)

df_sintetico$id <- 264:(264 + nrow(df_sintetico) - 1)
all(names(df) == names(df_sintetico))  # asegurar de que ambas bases tengan las mismas columnas = Esto debe devolver TRUE
df <- bind_rows(df, df_sintetico) # combinamos
any(duplicated(df$id))  # duplicados en id = debe ser FALSE



