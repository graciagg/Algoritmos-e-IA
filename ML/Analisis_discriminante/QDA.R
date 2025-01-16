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




#### ---- QDA (no lineal) ----

# Ajustar el modelo QDA en el entrenamiento
qda_model <- qda(formula_trat, data = training_trat)
qda_model


# Realizar predicciones sobre el conjunto de prueba
qda_predictions <- predict(qda_model, newdata = testing_trat)

# Obtener la predicción (predicciones de la clase)
predicted_classes <- qda_predictions$class
predicted_classes
length(predicted_classes)

# Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(testing_trat$trat)
true_classes
length(true_classes)

# Crear la matriz de confusión (trat predicho testing vs. trat real testing)
confusion <- confusionMatrix(predicted_classes, true_classes)
print(confusion)

# Acuracy: El modelo tiene una precisión del 50%. Esto significa que el 50% de las predicciones fueron correctas (el modelo acertó en la clasificación).
# Sensitivity: La sensibilidad mide cuántas veces el modelo predijo correctamente la clase 1 entre las veces que realmente era clase 1. En este caso es 0.0, lo que significa que el modelo no detectó ninguna instancia de clase 1 correctamente (ningún verdadero positivo para la clase 1).
# Specificity: La especificidad mide cuántas veces el modelo predijo correctamente que no era clase 1. Aquí es 1.0, lo que significa que el modelo predijo correctamente todas las veces que no era clase 1 (es decir, todas las veces que realmente era clase 2).
# Pos Pred Value: Este valor indica el porcentaje de veces que el modelo predijo clase 1 correctamente cuando hizo esa predicción. Dado que no hay verdaderos positivos (TP = 0), este valor es NaN (no disponible).
# Neg Pred Value: Mide cuántas veces el modelo predijo correctamente que no era clase 1 cuando realmente no era clase 1. En este caso es 0.5, lo que significa que el modelo predijo correctamente 50% de las veces que no era clase 1.





