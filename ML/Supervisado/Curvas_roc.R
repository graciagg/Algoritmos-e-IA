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




#### ---- Support Vector Machine ----
# Crear un modelo de SVM lineal utilizando el paquete caret
# parámetro C "cost" por defecto es 1, pero puedes tunearlo. Controla la flexibilidad del modelo para encontrar un equilibrio entre un margen amplio y la clasificación correcta de las muestras
svmModelLineal <- train(primaryormetastasis ~.,
                        data = trainData,
                        method = "svmLinear",
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneGrid = expand.grid(C = seq(0, 2, length = 20)), #C grande lleva al sobreajuste, C pequeño al infraajuste
                        prob.model = TRUE) 
svmModelLineal

plot(svmModelLineal)

# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(svmModelLineal, newdata = testData )
predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions, testData$primaryormetastasis)

# SVM lineal
probabilities_svm_linear <- predict(svmModelLineal, newdata = testData, type = "prob")
probabilities_svm_linear




# Crear un modelo de SVM tipo kernel utilizando el paquete caret
# no hace falta tunear el parámetro C "cost" 
svmModelKernel <- train(primaryormetastasis ~.,
                        data = trainData,
                        method = "svmRadial",
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneLength = 10,
                        prob.model = TRUE) 
svmModelKernel

plot(svmModelKernel)


# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(svmModelKernel, newdata = testData )
predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions, testData$primaryormetastasis)

# SVM kernel
probabilities_svm_kernel <- predict(svmModelKernel, newdata = testData, type = "prob")
probabilities_svm_kernel



# Crear un modelo de SVM tipo kernel polynomial utilizando el paquete caret
# no hace falta tunear el parámetro C "cost" 
svmModelKernelPolynomial <- train(primaryormetastasis ~.,
                                  data = trainData,
                                  method = "svmPoly",
                                  trControl = trainControl(method = "cv", number = 10),
                                  preProcess = c("center", "scale"),
                                  tuneLength = 5,
                                  prob.model = TRUE) 
svmModelKernelPolynomial

plot(svmModelKernelPolynomial)


# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(svmModelKernelPolynomial, newdata = testData )
predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions, testData$primaryormetastasis)

# SVM kernel
probabilities_svm_kernelpol <- predict(svmModelKernelPolynomial, newdata = testData, type = "prob")
probabilities_svm_kernelpol



#### ---- Decission Tree ----
# Crear un modelo de DT utilizando el paquete caret
dtModel <- train(primaryormetastasis ~.,
                 data = trainData,
                 method = "rpart",
                 trControl = trainControl(method = "cv", number = 10),
                 preProcess = c("center", "scale"),
                 tuneLength = 10)
dtModel
plot(dtModel)

fancyRpartPlot(dtModel$finalModel, type=4)


# Evaluar el modelo con el conjunto de prueba
predictions_raw <- predict(dtModel, newdata = testData, type = "raw") # raw = clases
predictions_raw


# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions_raw, testData$primaryormetastasis)

# Obtener probabilidades
probabilities_dt <- predict(dtModel, newdata = testData, type = "prob")







#### ---- LDA (lineal) ----
formula <- as.formula(paste("primaryormetastasis ~", paste(names, collapse = "+")))
formula

# Ajustar el modelo LDA en el entrenamiento
lda_model <- lda(formula, data = trainData)
lda_model$scaling # contribuciones/coeficientes

# Realizar predicciones sobre el conjunto de prueba
lda_predictions <- predict(lda_model, newdata = testData)
lda_predictions$x

predicted_classes <- lda_predictions$class # Obtener la predicción (predicciones de la clase)
true_classes <- as.factor(testData$primaryormetastasis) # Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
confusionMatrix(predicted_classes, true_classes) # Crear la matriz de confusión (tumor predicho testing vs. tumor real testing)
probabilities_lda <- predict(lda_model, newdata = testData, type = "prob") # Obtener probabilidades






#### ---- RDA (regularizado) ----
# Ajustar el modelo RDA en el entrenamiento
rda_model <- rda(formula, data = trainData)

# Realizar predicciones sobre el conjunto de prueba
rda_predictions <- predict(rda_model, newdata = testData)

predicted_classes <- rda_predictions$class # Obtener la predicción (predicciones de la clase)
true_classes <- as.factor(testData$primaryormetastasis) # Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
confusionMatrix(predicted_classes, true_classes) # Crear la matriz de confusión (tumor predicho testing vs. tumor real testing)
probabilities_rda <- predict(rda_model, newdata = testData, type = "prob") # Obtener probabilidades







#### ---- Curvas ROC ----
roc_knn <- roc(testData$primaryormetastasis, probabilities_knn[,2]) # Cambia [,2] según la clase positiva
auc_knn <- auc(roc_knn)
cat("AUC k-NN:", auc_knn, "\n")

roc_svm_linear <- roc(testData$primaryormetastasis, probabilities_svm_linear[,2])
auc_svm_linear <- auc(roc_svm_linear)
cat("AUC SVM Lineal:", auc_svm_linear, "\n")

roc_svm_kernel <- roc(testData$primaryormetastasis, probabilities_svm_kernel[,2])
auc_svm_kernel <- auc(roc_svm_kernel)
cat("AUC SVM Kernel:", auc_svm_kernel, "\n")


roc_svm_kernelpol <- roc(testData$primaryormetastasis, probabilities_svm_kernelpol[,2])
auc_svm_kernelpol <- auc(roc_svm_kernelpol)
cat("AUC SVM Kernel Polynomial:", auc_svm_kernelpol, "\n")

roc_dt <- roc(testData$primaryormetastasis, probabilities_dt[,2])
auc_dt <- auc(roc_dt)
cat("AUC Árbol de decisión:", auc_dt, "\n")

roc_lda <- roc(testData$primaryormetastasis, probabilities_lda$posterior[, 2])
auc_lda <- auc(roc_lda)
cat("AUC LDA:", auc_lda, "\n")

roc_rda <- roc(testData$primaryormetastasis, probabilities_rda$posterior[, 2])
auc_rda <- auc(roc_rda)
cat("AUC QDA:", auc_rda, "\n")





# Asegúrate de que la primera curva ROC sea la base
plot(roc_knn, col = "blue", main = "Curvas ROC", lwd = 2)
plot(roc_svm_linear, col = "red", add = TRUE, lwd = 2)
plot(roc_svm_kernel, col = "green", add = TRUE, lwd = 2)
plot(roc_svm_kernelpol, col = "orange", add = TRUE, lwd = 2)
plot(roc_dt, col = "purple", add = TRUE, lwd = 2)
plot(roc_lda, col = "pink", add = TRUE, lwd = 2)
plot(roc_rda, col = "yellow", add = TRUE, lwd = 2)


# Agregar leyenda
knn_legend <- paste("AUC k-NN:", round(auc_knn, 2))  # Redondeamos a 2 decimales, si es necesario
svm_legend <- paste("AUC SVM Lineal:", round(auc_svm_linear, 2))  # Redondeamos a 2 decimales, si es necesario
svmk_legend <- paste("AUC SVM Kernel:", round(auc_svm_kernel, 2))  # Redondeamos a 2 decimales, si es necesario
svmp_legend <- paste("AUC SVM Kernel Polynomial:", round(auc_svm_kernelpol, 2))  # Redondeamos a 2 decimales, si es necesario
dt_legend <- paste("AUC Decission Tree:", round(auc_dt, 2))  # Redondeamos a 2 decimales, si es necesario
lda_legend <- paste("AUC LDA:", round(auc_lda, 2))  # Redondeamos a 2 decimales, si es necesario
rda_legend <- paste("AUC RDA:", round(auc_rda, 2))  # Redondeamos a 2 decimales, si es necesario

legend("bottomright", legend = c(knn_legend, svm_legend, svmk_legend, svmp_legend, dt_legend, lda_legend, rda_legend),
       col = c("blue", "red", "green", "orange", "purple" ,"pink" ,"yellow" ), lwd = 2)