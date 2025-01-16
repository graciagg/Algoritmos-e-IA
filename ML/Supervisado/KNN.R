#### ---- Preparamos la base de datos para modelos de clasificación ----

rm(list=ls())
#path <- '/Users/vic/Library/CloudStorage/GoogleDrive-vdelaopascual@gmail.com/Mi unidad/MU en Bioinformática (UNIR 2023)/Presentaciones Algoritmos e Inteligencia Artificial (AAAA-MM-DD)/Ejemplo clase expresion genica'
setwd("~/BIOINFORMATICA_UNIR/Algoritmos e Inteligencia Artificial/Tema 6 Aprendizaje supervisado")
df <- read.csv("dataset_expresiongenes_cancer.csv")
str(df)
table(df$primaryormetastasis)
data <- subset(df, select = 3:101)
str(data)


library(glmnet) # ElasticNet regression
library(tidyverse)
library(caret) # ML
library(rpart) # DT
library(rpart.plot) # Decision Tree plot
library(rattle) # DT plot
library(pROC) # ROC
library(PRROC) # PR-Curve
library(MASS) # LDA
library(klaR) # RDA
library(gridExtra) # juntar los gráficos

# Machine learning methods: http://topepo.github.io/caret/train-models-by-tag.html
# in method="XXXX"
# in metric="XXXX" -> "RMSE" para regression y "Accuracy" para classification
# Apoyo: https://rpubs.com/nomarpicasso/1150167
names(getModelInfo())
modelLookup(model = "knn")
modelLookup(model = "svmLinear")
modelLookup(model = "svmRadial")
modelLookup(model = "svmPoly")
modelLookup(model = "rpart")








genes <- names(df[4:19197])


# Preparar los datos para el modelo LASSO (para q elimine variables q no sean tan importantes y seleccione genes q contribuyan más a la clasificación)
# Todo esto lo hago para reducir genes, pq tengo 19.000 y pico
x <- as.matrix(df[, genes])
y <- factor(df$primaryormetastasis)

# Ajustar el modelo LASSO
set.seed(1995)
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)
selected_genes <- coef(lasso_model, s = "lambda.min")
selected_genes <- as.matrix(selected_genes) # Convertir a matriz densa si es necesario (esto convierte el formato disperso a un formato manejable)
selected_genes_df <- as.data.frame(selected_genes) # Convertir la matriz a data frame para facilitar su manipulación
non_zero_indices <- selected_genes_df[selected_genes_df$s1 != 0, , drop = FALSE] # Filtrar los coeficientes no cero
dim(non_zero_indices)
non_zero_indices

names <- rownames(non_zero_indices)[2:180]
names

data <- df %>% dplyr:: select(primaryormetastasis, names)
rows <- df$code
rownames(data) <- rows
str(data)

write_csv(data, "~/BIOINFORMATICA_UNIR/Algoritmos e Inteligencia Artificial/Tema 6 Aprendizaje supervisado/stata.csv")
names <- colnames(data)[-1]
data <- data %>% dplyr::select(names, primaryormetastasis)


# Dividir el conjunto de datos en conjuntos de entrenamiento y prueba
set.seed(1995)
trainIndex <- createDataPartition(data$primaryormetastasis, p = 0.8, list = FALSE)
data$primaryormetastasis <- as.factor(data$primaryormetastasis)
trainData <- data[trainIndex,]
testData <- data[-trainIndex,]


#### ---- kNN ----
# Crear un modelo de k-NN utilizando el paquete caret
knnModel <- train(primaryormetastasis ~ .,
                  data = trainData,
                  method = "knn",
                  trControl = trainControl(method = "cv", number = 10),
                  preProcess = c("center", "scale"),
                  tuneLength = 30) #Que me pruebe una longitud de 30 parámetros de K. No tienen por qué ir del 1 al 30. Él empieza por 5 y va subiendo a números impares
knnModel

plot(knnModel)

# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(knnModel, newdata = testData )
predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions, testData$primaryormetastasis)


# Obtener probabilidades
probabilities_knn <- predict(knnModel, newdata = testData, type = "prob")
probabilities_knn