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


######################## MODELOS DE REGRESION LINEAL

#### ---- Preparamos la base de datos para modelos de regresión ----

path2 <- '/Users/vic/Library/CloudStorage/GoogleDrive-vdelaopascual@gmail.com/Mi unidad/MU en Bioinformática (UNIR 2023)/Presentaciones Algoritmos e Inteligencia Artificial (AAAA-MM-DD)/Otras bases de datos/Framingham Dataset_dep.csv'
df2 <- read.csv(path2)

str(df2)
anyNA(df2)
df2 <- na.omit(df2)


# Eliminar 'code' y configurar las variables
df2 <- df2[, !names(df2) %in% c("code")]
df2 <- df2[df2$antihypertensive_use != "", ]
df2$sex <- as.factor(df2$sex)
df2$current_smoker <- as.factor(df2$current_smoker)
df2$diabetes_status <- as.factor(df2$diabetes_status)
df2$antihypertensive_use <- as.factor(df2$antihypertensive_use)
df2$prev_coronary_hd <- as.factor(df2$prev_coronary_hd)
df2$prev_angina <- as.factor(df2$prev_angina)
df2$prev_myocardial_infarction <- as.factor(df2$prev_myocardial_infarction)
df2$prev_stroke <- as.factor(df2$prev_stroke)
df2$prev_hypertension <- as.factor(df2$prev_hypertension)
df2$exam_cycle <- as.factor(df2$exam_cycle)
df2$death_status <- as.factor(df2$death_status)
df2$followup_angina <- as.factor(df2$followup_angina)
df2$hosp_myocardial_infarction <- as.factor(df2$hosp_myocardial_infarction)
df2$hosp_mi_or_fatal_chd <- as.factor(df2$hosp_mi_or_fatal_chd)
df2$followup_chd <- as.factor(df2$followup_chd)
df2$followup_stroke <- as.factor(df2$followup_stroke)
df2$followup_cvd <- as.factor(df2$followup_cvd)
df2$followup_hypertension <- as.factor(df2$followup_hypertension)
str(df2)
dim(df2)
df2 <- df2[, -c(18, 29:36)]


# Separar conjuntos de entrenamiento y prueba
set.seed(1995)
trainIndex <- createDataPartition(df2$total_cholesterol, p = 0.8, list = FALSE)
trainData_2 <- df2[trainIndex, ]
testData_2 <- df2[-trainIndex, ]




#### ---- kNN regression ----
knnModel_reg <- train(total_cholesterol ~.,
                      data = trainData_2,
                      method = "knn",
                      trControl = trainControl(method = "cv", number = 10),
                      preProcess = c("center", "scale"),
                      tuneLength = 30)
knnModel_reg
plot(knnModel_reg)

# Predicciones y métricas
knn_preds <- predict(knnModel_reg, newdata = testData_2)
postResample(knn_preds, testData_2$total_cholesterol)
knn_rmse <- sqrt(mean((knn_preds - testData_2$total_cholesterol)^2))  # RMSE para kNN

knn_preds_df <- data.frame(Real = testData_2$total_cholesterol,
                           Predicted = knn_preds)

knnModel_graph <- ggplot(knn_preds_df, aes(x = Real, y = Predicted)) +
  geom_point(alpha = 0.6) +  # Puntos de la dispersión
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Línea de igualdad
  labs(title = paste("Modelo kNN - RMSE:", round(knn_rmse, 2)),
       x = "Valores Reales", 
       y = "Predicciones") +
  theme_minimal()

knnModel_graph



#### ---- SVM regression ----
svmModelLineal_reg <- train(total_cholesterol ~.,
                            data = trainData_2,
                            method = "svmLinear",
                            trControl = trainControl(method = "cv", number = 10),
                            preProcess = c("center", "scale"),
                            tuneGrid = expand.grid(C = seq(0, 2, length = 10))) #C grande lleva al sobreajuste, C pequeño al infraajuste
svmModelLineal_reg
plot(svmModelLineal_reg)

# Predicciones y métricas
svmModelLineal_preds <- predict(svmModelLineal_reg, newdata = testData_2)
postResample(svmModelLineal_preds, testData_2$total_cholesterol)
svm_rmse <- sqrt(mean((svmModelLineal_preds - testData_2$total_cholesterol)^2))  # RMSE para SVM

svmModelLineal_preds_df <- data.frame(Real = testData_2$total_cholesterol,
                                      Predicted = svmModelLineal_preds)

svmModelLineal_reg_graph <- ggplot(svmModelLineal_preds_df, aes(x = Real, y = Predicted)) +
  geom_point(alpha = 0.6) +  # Puntos de la dispersión
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Línea de igualdad
  labs(title = paste("Modelo SVM lineal - RMSE:", round(svm_rmse, 2)),
       x = "Valores Reales", 
       y = "Predicciones") +
  theme_minimal()






svmModelRadial_reg <- train(total_cholesterol ~.,
                            data = trainData_2,
                            method = "svmRadial",
                            trControl = trainControl(method = "cv", number = 10),
                            preProcess = c("center", "scale"),
                            tuneLength = 10)
svmModelRadial_reg
plot(svmModelRadial_reg)

# Predicciones y métricas
svmModelRadial_preds <- predict(svmModelRadial_reg, newdata = testData_2)
postResample(svmModelRadial_preds, testData_2$total_cholesterol)
svm_kernel_rmse <- sqrt(mean((svmModelRadial_preds - testData_2$total_cholesterol)^2))  # RMSE para SVM Kernel

svmModelRadial_preds_df <- data.frame(Real = testData_2$total_cholesterol,
                                      Predicted = svmModelRadial_preds)

svmModelRadial_reg_graph <- ggplot(svmModelRadial_preds_df, aes(x = Real, y = Predicted)) +
  geom_point(alpha = 0.6) +  # Puntos de la dispersión
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Línea de igualdad
  labs(title = paste("Modelo SVM Kernel - RMSE:", round(svm_kernel_rmse, 2)),
       x = "Valores Reales", 
       y = "Predicciones") +
  theme_minimal()



#### ---- Decission Tree regression ----
dtModel_reg <- train(total_cholesterol ~.,
                     data = trainData_2,
                     method = "rpart",
                     trControl = trainControl(method = "cv", number = 10),
                     preProcess = c("center", "scale"),
                     tuneGrid = expand.grid(cp = seq(0.0001, 0.5, by = 0.005)))  # Más valores para cp
dtModel_reg
plot(dtModel_reg)

# Predicciones y métricas
dtModel_preds <- predict(dtModel_reg, newdata = testData_2)
postResample(dtModel_preds, testData_2$total_cholesterol)
dt_rmse <- sqrt(mean((dtModel_preds - testData_2$total_cholesterol)^2))  # RMSE para Decision Tree

dtModel_preds_df <- data.frame(Real = testData_2$total_cholesterol,
                               Predicted = dtModel_preds)

dt_reg_graph <- ggplot(dtModel_preds_df, aes(x = Real, y = Predicted)) +
  geom_point(alpha = 0.6) +  # Puntos de la dispersión
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Línea de igualdad
  labs(title = paste("Modelo Decision Tree - RMSE:", round(dt_rmse, 2)),
       x = "Valores Reales", 
       y = "Predicciones") +
  theme_minimal()


graphs_ml_reg <- list(knnModel_graph, svmModelLineal_reg_graph, svmModelRadial_reg_graph, dt_reg_graph)
grid.arrange(grobs = graphs_ml_reg, ncol = 2 )

# Importance plots
plot(varImp(knnModel_reg, scale = TRUE)) # importancia de las variables
plot(varImp(svmModelLineal_reg, scale = TRUE)) # importancia de las variables
plot(varImp(svmModelRadial_reg, scale = TRUE)) # importancia de las variables
plot(varImp(dtModel_reg, scale = TRUE)) # importancia de las variables


