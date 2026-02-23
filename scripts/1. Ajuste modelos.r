
# Cargar de todas las librerías necesarias para el análisis
library(tidymodels)
library(ggplot2)
library(dplyr)
library(factoextra)
library(glmnet)
library(broom)
library(rsample)
library(discrim)
library(vip)
library(yardstick)
library(patchwork)
library(workflows)
library(shapviz)
library(tidyr)
library(klaR)
library(ranger)
library(kernlab)
library(ggbeeswarm)
library(boot)

# Establecer preferencias para el ecosistema tidymodels
tidymodels_prefer()

#semilla utilizada set.seed(123)


# Lectura de los datos desde archivo de texto
Dent <- read.table("Dientes.txt", header = TRUE)

# Exploración inicial de los datos
head(Dent, 20)
summary(Dent)
str(Dent)

# Conversión de variables categóricas a factores
Dent$Ind <- as.factor(Dent$Ind)
Dent$Sexo <- as.factor(Dent$Sexo)


# Filtrar datos para caninos superiores (CS)
CS <- Dent %>%
  filter(Tipo_Diente == "CS")

# Exploración de los datos de CS
summary(CS[, c("MD", "BL", "MDCu", "BLCu")])
table(CS$Sexo)
length(unique(CS$Ind))

# Guardar datos de CS para análisis posteriores
write.table(
  CS,
  file = "caninos_superiores.txt",
  sep = "\t",
  row.names = FALSE
)

# Función para dividir datos manteniendo la integridad por individuo y sexo
split_dentales <- function(data, prop = 0.8) {
  
  # Verificar columnas requeridas
  required_cols <- c("Ind", "Sexo")
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    stop(paste("Faltan columnas:", 
               paste(missing, collapse = ", ")))
  }
  
  set.seed(123)  # Para reproducibilidad
  
  # División estratificada por sexo y agrupada por individuo
  split_obj <- group_initial_split(
    data,
    prop = prop,
    group = Ind,     # Para evitar pseudorreplicación
    strata = Sexo    # Para mantener balance por sexo
  )
  
  # Devolver lista con conjuntos de entrenamiento y prueba
  list(
    train = training(split_obj),
    test  = testing(split_obj),
    split = split_obj
  )
}

# Aplicar la función de división a los datos de caninos superiores
cs_split <- split_dentales(CS)
cs_train <- cs_split$train  # conjunto de entrenamiento
cs_test <- cs_split$test    # conjunto de prueba

# Verificar la división
summary(cs_train)
summary(cs_test)

# Verificar que no haya individuos repetidos entre conjuntos
length(intersect(unique(cs_train$Ind), unique(cs_test$Ind)))

# Identificar casos con todos los valores NA
cs_train %>%
  filter(is.na(MD) & is.na(BL) & is.na(MDCu) & is.na(BLCu))

# Eliminar casos con todos los valores NA
cs_train <- cs_train %>%
  filter(!(is.na(MD) & is.na(BL) & is.na(MDCu) & is.na(BLCu)))

# Verificar que se hayan eliminado los casos NA
cs_train %>%
  filter(is.na(MD) & is.na(BL) & is.na(MDCu) & is.na(BLCu))

# Crear receta para preprocesamiento de datos
cs_recipe <- recipe(Sexo ~ MD + BL + MDCu + BLCu, data = cs_train) %>%
  step_impute_knn(all_predictors()) %>%  # Imputación KNN para valores faltantes
  step_normalize(all_predictors())        # Estandarización de predictores

# Verificar la receta preparada
prep(cs_recipe)

# Configurar validación cruzada estratificada y agrupada
set.seed(123)  # Para reproducibilidad

cv_cs <- group_vfold_cv(
  cs_train,
  v = 10,        # 10 folds
  repeats = 5,   # 5 repeticiones
  group = Ind,   # Agrupar por individuo
  strata = Sexo  # Estratificar por sexo
)

# Definir métricas para evaluar el rendimiento de los modelos
metricas <- metric_set(
  accuracy,       # Precisión general
  sens,           # Sensibilidad (tasa de verdaderos positivos)
  spec,           # Especificidad (tasa de verdaderos negativos)
  bal_accuracy,   # Balanced Accuracy (promedio de sensibilidad y especificidad)
  roc_auc         # Área bajo la curva ROC
)

# Configurar control para resampling
ctrl <- control_resamples(save_pred = TRUE)

# Crear workflow básico con la receta definida
cs_wf <- workflow() %>%
  add_recipe(cs_recipe)

# Especificar modelo de regresión logística
log_spec <- logistic_reg(mode = "classification") %>%
  set_engine("glm")

# Crear workflow completo
log_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(log_spec)

# Evaluar con validación cruzada
log_res <- fit_resamples(
  log_wf,
  resamples = cv_cs,
  control = ctrl,
  metrics = metricas
)

# Mostrar métricas
collect_metrics(log_res)

# Especificar modelo Elastic Net
enet_spec <- logistic_reg(
  penalty = tune(),  # Parámetro de regularización a optimizar
  mixture = tune()   # Tipo de regularización (0=ridge, 1=lasso)
) %>%
  set_engine("glmnet") %>%  # Motor de ajuste
  set_mode("classification")

# Crear workflow
enet_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(enet_spec)

# Definir grilla de hiperparámetros
enet_grid <- grid_regular(
  penalty(range = c(-4, 0)),  # Valores de lambda en escala log
  mixture(range = c(0, 1)),   # Valores de alpha (0=ridge, 1=lasso)
  levels = 10                 # 10x10 combinaciones (100 modelos)
)

# Optimizar hiperparámetros con validación cruzada
set.seed(123)
enet_res <- tune_grid(
  enet_wf,
  resamples = cv_cs,
  grid = enet_grid,
  control = ctrl,
  metrics = metricas
)

# Seleccionar la mejor configuración
best_enet <- select_best(enet_res, metric = "bal_accuracy")

# Especificar modelo LDA
lda_spec <- discrim_linear() %>%
  set_engine("MASS")

# Crear workflow
lda_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(lda_spec)

# Evaluar con validación cruzada
lda_res <- fit_resamples(
  lda_wf,
  resamples = cv_cs,
  metrics = metricas,
  control = ctrl
)

# Mostrar métricas
collect_metrics(lda_res)

# Especificar modelo Naive Bayes
nb_spec <- naive_Bayes() %>%
  set_engine("klaR")  # Usa paquete klaR para el ajuste

# Crear workflow
nb_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(nb_spec)

# Evaluar con validación cruzada
nb_res <- fit_resamples(
  nb_wf,
  resamples = cv_cs,
  metrics = metricas,
  control = ctrl
)

# Mostrar métricas
collect_metrics(nb_res)

# Especificar modelo Random Forest
rf_spec <- rand_forest(
  trees = 1000,    # Número de árboles
  mtry = tune(),   # Número de predictores a considerar en cada división
  min_n = tune()   # Tamaño mínimo de nodo
) %>%
  set_engine("ranger") %>%
  set_mode("classification")

# Crear workflow
rf_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(rf_spec)

# Definir grilla de hiperparámetros
rf_grid <- grid_regular(
  mtry(range = c(1, 4)),    # 1-4 predictores por división
  min_n(range = c(2, 20)),  # Tamaño mínimo de nodo
  levels = 6                # 6x6 combinaciones (~36 modelos)
)

# Optimizar hiperparámetros
set.seed(123)
rf_res <- tune_grid(
  rf_wf,
  resamples = cv_cs,
  grid = rf_grid,
  metrics = metricas,
  control = ctrl
)

# Seleccionar la mejor configuración
best_rf <- select_best(rf_res, metric = "bal_accuracy")

# Especificar modelo SVM lineal
svm_lin <- svm_linear(cost = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# Crear workflow
svm_lin_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(svm_lin)

# Definir grilla de hiperparámetros
svm_lin_grid <- grid_regular(
  cost(range = c(-3, 2)),  # Valores de costo en escala log
  levels = 10
)

# Optimizar hiperparámetros
set.seed(123)
svm_lin_res <- tune_grid(
  svm_lin_wf,
  resamples = cv_cs,
  grid = svm_lin_grid,
  metrics = metricas,
  control = ctrl
)

# Seleccionar la mejor configuración
best_svm_lin <- select_best(svm_lin_res, metric = "bal_accuracy")

# Especificar modelo SVM con kernel radial
svm_rad <- svm_rbf(
  cost = tune(),
  rbf_sigma = tune()  # Parámetro sigma del kernel RBF
) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# Crear workflow
svm_rad_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(svm_rad)

# Definir grilla de hiperparámetros
svm_rad_grid <- grid_regular(
  cost(range = c(-2, 2)),        # Valores de costo (0.01 a 100)
  rbf_sigma(range = c(-3, -1)),  # Valores de sigma (rango razonable)
  levels = 6
)

# Optimizar hiperparámetros
set.seed(123)
svm_rad_res <- tune_grid(
  svm_rad_wf,
  resamples = cv_cs,
  grid = svm_rad_grid,
  metrics = metricas,
  control = ctrl
)

# Seleccionar la mejor configuración
best_svm_rad <- select_best(svm_rad_res, metric = "bal_accuracy")

# Especificar modelo XGBoost
xgb_spec <- boost_tree(
  trees = 800,             # Número de árboles
  learn_rate = tune(),     # Tasa de aprendizaje
  mtry = tune(),           # Número de predictores por división
  tree_depth = tune(),     # Profundidad máxima de los árboles
  min_n = tune(),          # Tamaño mínimo de nodo
  loss_reduction = tune()  # Reducción mínima de pérdida
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# Crear workflow
xgb_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(xgb_spec)

# Definir grilla de hiperparámetros (muestreo espacial)
xgb_grid <- grid_space_filling(
  learn_rate(range = c(-3, -1)),      # 0.001 a 0.1
  mtry(range = c(1, 4)),
  tree_depth(range = c(2L, 6L)),      # Profundidad moderada
  min_n(range = c(2L, 15L)),
  loss_reduction(range = c(-5, -1)),  # 1e-5 a 0.1
  size = 25                           # 25 combinaciones
)

# Optimizar hiperparámetros
set.seed(123)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = cv_cs,
  grid = xgb_grid,
  metrics = metricas,
  control = ctrl
)

# Seleccionar la mejor configuración
best_xgb <- select_best(xgb_res, metric = "bal_accuracy")

# Finalizar workflow y entrenar modelo final
xgb_final_wf <- finalize_workflow(xgb_wf, best_xgb)
xgb_fit <- fit(xgb_final_wf, data = cs_train)

# Función auxiliar para modelos sin tuning
extraer_simple <- function(res_obj, nombre_modelo) {
  collect_metrics(res_obj, summarize = FALSE) %>%
    dplyr::filter(.metric == "bal_accuracy") %>%
    dplyr::summarise(
      mean_BA = mean(.estimate),
      sd_BA   = sd(.estimate)
    ) %>%
    dplyr::mutate(modelo = nombre_modelo)
}

# Función auxiliar para modelos con tuning
extraer_tuneado <- function(res_obj, nombre_modelo) {
  # Identificar mejor configuración
  best_id <- select_best(res_obj, metric = "bal_accuracy")$.config
  
  collect_metrics(res_obj, summarize = FALSE) %>%
    dplyr::filter(.metric == "bal_accuracy",
                  .config == best_id) %>%
    dplyr::summarise(
      mean_BA = mean(.estimate),
      sd_BA   = sd(.estimate)
    ) %>%
    dplyr::mutate(modelo = nombre_modelo)
}

# Recolección robusta de métricas
comparacion_modelos <- bind_rows(
  extraer_simple(log_res, "Logistic"),
  extraer_tuneado(enet_res, "Elastic Net"),
  extraer_simple(lda_res, "LDA"),
  extraer_simple(nb_res, "Naive Bayes"),
  extraer_tuneado(rf_res, "Random Forest"),
  extraer_tuneado(svm_lin_res, "SVM Linear"),
  extraer_tuneado(svm_rad_res, "SVM Radial"),
  extraer_tuneado(xgb_res, "XGBoost")
) %>%
  dplyr::arrange(desc(mean_BA))

# Generar predicciones (clase + probabilidades)
xgb_pred <- predict(xgb_fit, cs_test, type = "prob") %>%
  bind_cols(predict(xgb_fit, cs_test)) %>%
  bind_cols(dplyr::select(cs_test, Sexo))

# Asegurar que Sexo sea factor
xgb_pred$Sexo <- as.factor(xgb_pred$Sexo)

# Definir métricas principales
metricas_clase <- metric_set(
  accuracy,
  sens,
  spec,
  bal_accuracy
)

# Calcular métricas
metricas_clase(xgb_pred, truth = Sexo, estimate = .pred_class)

# Calcular AUC-ROC
roc_auc(xgb_pred, truth = Sexo, .pred_F)

# Generar y graficar la curva ROC
roc_curve(xgb_pred, truth = Sexo, .pred_F) %>%
  autoplot()



# Función para extraer métricas en cada remuestreo
boot_metricas <- function(data, indices) {
  d <- data[indices, ]
  c(
    bal_acc = bal_accuracy(d, truth = Sexo, estimate = .pred_class)$.estimate,
    sens    = sens(d, truth = Sexo, estimate = .pred_class)$.estimate,
    spec    = spec(d, truth = Sexo, estimate = .pred_class)$.estimate,
    acc     = accuracy(d, truth = Sexo, estimate = .pred_class)$.estimate
  )
}

# Realizar bootstrap
set.seed(123)
boot_res <- boot(
  data = xgb_pred,
  statistic = boot_metricas,
  R = 2000
)

# Calcular intervalos de confianza al 95%
boot.ci(boot_res, type = "perc", index = 1)  # Balanced Accuracy
boot.ci(boot_res, type = "perc", index = 2)  # Sensibilidad
boot.ci(boot_res, type = "perc", index = 3)  # Especificidad
boot.ci(boot_res, type = "perc", index = 4)  # Accuracy

# Crear tabla con sensibilidad y especificidad por umbral
threshold_df <- roc_curve(xgb_pred, truth = Sexo, .pred_F) %>%
  dplyr::select(.threshold, sensitivity, specificity)

# Graficar sensibilidad y especificidad por umbral
ggplot(threshold_df, aes(x = .threshold)) +
  geom_line(aes(y = sensitivity, color = "Sensibilidad"), size = 1) +
  geom_line(aes(y = specificity, color = "Especificidad"), size = 1) +
  scale_color_manual(values = c("Sensibilidad" = "red",
                                "Especificidad" = "blue")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Umbral",
    y = "Valor",
    color = "",
    title = "Sensibilidad y Especificidad según Umbral"
  )

# Umbral más balanceado (igualando sensibilidad y especificidad)
best_threshold <- threshold_df %>%
  dplyr::mutate(diff = abs(sensitivity - specificity)) %>%
  dplyr::arrange(diff) %>%
  dplyr::slice(1)

# Umbral de Youden J (maximizando sensibilidad + especificidad - 1)
best_youden <- threshold_df %>%
  mutate(J = sensitivity + specificity - 1) %>%
  arrange(desc(J)) %>%
  dplyr::slice(1)

# Generar matriz de confusión
conf_mat(xgb_pred, truth = Sexo, estimate = .pred_class)

# Calcular importancia de variables
xgb_importance <- xgb_fit %>%
  extract_fit_parsnip() %>%
  vip::vi()

# Graficar importancia de variables
vip(
  xgb_fit,
  num_features = 10,
  geom = "col",
  aesthetics = list(fill = "steelblue")
)

# Extraer modelo XGBoost
xgb_core <- extract_fit_parsnip(xgb_fit)$fit

# Obtener datos procesados por la receta
rec <- extract_recipe(xgb_fit)
X_df <- bake(rec, new_data = cs_train) %>%
  dplyr::select(-Sexo)

# Convertir a matriz numérica
X_matrix <- as.matrix(X_df)

# Calcular valores SHAP
sv <- shapviz(
  xgb_core,
  X = X_matrix,
  X_pred = X_matrix,
  baseline = "auto"
)

# Gráfico de importancia SHAP
sv_importance(sv, kind = "beeswarm")

# Gráficos de dependencia SHAP para variables clave
vars <- c("MDCu", "BLCu", "BL", "MD")
plots <- lapply(vars, function(v){
  sv_dependence(sv, v) +
    ggtitle(paste("Dependencia SHAP -", v)) +
    theme_minimal(base_size = 12)
})
wrap_plots(plots, ncol = 2)

# Filtrar datos para caninos inferiores (CI)
CI <- Dent %>%
  filter(Tipo_Diente == "CI")

# Exploración de los datos de CI
summary(CI[, c("MD", "BL", "MDCu", "BLCu")])
table(CI$Sexo)
length(unique(CI$Ind))

# Guardar datos de CI para análisis posteriores
write.table(
  CI,
  file = "caninos_inferiores.txt",
  sep = "\t",
  row.names = FALSE
)

# Dividir datos usando la función definida
ci_split <- split_dentales(CI)
ci_train <- ci_split$train
ci_test <- ci_split$test

# Verificar la división
summary(ci_train)
summary(ci_test)
length(intersect(unique(ci_train$Ind), unique(ci_test$Ind)))

# Receta para CI
ci_recipe <- recipe(Sexo ~ MD + BL + MDCu + BLCu, data = ci_train) %>%
  step_impute_knn(all_predictors()) %>%
  step_normalize(all_predictors())

# Configuración de validación cruzada para CI
set.seed(123)
cv_ci <- group_vfold_cv(
  ci_train,
  v = 10,
  repeats = 5,
  group = Ind,
  strata = Sexo
)

# Métricas para CI
metricas <- metric_set(
  accuracy,
  sens,
  spec,
  bal_accuracy,
  roc_auc
)
ctrl <- control_resamples(save_pred = TRUE)

# Especificar modelo de regresión logística
log_spec_CI <- logistic_reg(mode = "classification") %>%
  set_engine("glm")

# Crear workflow completo
log_wf_CI <- workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(log_spec_CI)

# Evaluar con validación cruzada
log_res_CI <- fit_resamples(
  log_wf_CI,
  resamples = cv_ci,
  control = ctrl,
  metrics = metricas
)

# Mostrar métricas
collect_metrics(log_res_CI)

# Especificar modelo LDA
lda_spec <- discrim_linear() %>%
  set_engine("MASS")

# Crear workflow
lda_wf_CI <- workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(lda_spec)

# Evaluar con validación cruzada
lda_res_CI <- fit_resamples(
  lda_wf_CI,
  resamples = cv_ci,
  metrics = metricas,
  control = ctrl
)

# Mostrar métricas
collect_metrics(lda_res_CI)

# Especificar modelo Naive Bayes
nb_spec <- naive_Bayes() %>%
  set_engine("klaR")

# Crear workflow
nb_wf_CI <- workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(nb_spec)

# Evaluar con validación cruzada
nb_res_CI <- fit_resamples(
  nb_wf_CI,
  resamples = cv_ci,
  metrics = metricas,
  control = ctrl
)

# Mostrar métricas
collect_metrics(nb_res_CI)

# Especificar modelo Random Forest
rf_spec <- rand_forest(
  trees = 1000,
  mtry = tune(),
  min_n = tune()
) %>%
  set_engine("ranger") %>%
  set_mode("classification")

# Crear workflow
rf_wf_CI <- workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(rf_spec)

# Definir grilla de hiperparámetros
rf_grid_CI <- grid_regular(
  mtry(range = c(1, 4)),
  min_n(range = c(2, 20)),
  levels = 6
)

# Optimizar hiperparámetros
set.seed(123)
rf_res_CI <- tune_grid(
  rf_wf_CI,
  resamples = cv_ci,
  grid = rf_grid_CI,
  metrics = metricas,
  control = ctrl
)

# Seleccionar la mejor configuración
best_rf_CI <- select_best(rf_res_CI, metric = "bal_accuracy")

# Especificar modelo Elastic Net
enet_spec <- logistic_reg(
  penalty = tune(),
  mixture = tune()
) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# Crear workflow
enet_wf_CI <- workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(enet_spec)

# Definir grilla de hiperparámetros
enet_grid <- grid_regular(
  penalty(range = c(-4, 0)),
  mixture(range = c(0, 1)),
  levels = 10
)

# Optimizar hiperparámetros
set.seed(123)
enet_res_CI <- tune_grid(
  enet_wf_CI,
  resamples = cv_ci,
  grid = enet_grid,
  control = ctrl,
  metrics = metricas
)

# Seleccionar la mejor configuración
best_enet_CI <- select_best(enet_res_CI, metric = "bal_accuracy")

# Especificar modelo SVM lineal
svm_lin <- svm_linear(cost = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# Crear workflow
svm_lin_wf_CI <- workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(svm_lin)

# Definir grilla de hiperparámetros
svm_lin_grid <- grid_regular(
  cost(range = c(-3, 2)),
  levels = 10
)

# Optimizar hiperparámetros
set.seed(123)
svm_lin_res_CI <- tune_grid(
  svm_lin_wf_CI,
  resamples = cv_ci,
  grid = svm_lin_grid,
  metrics = metricas,
  control = ctrl
)

# Seleccionar la mejor configuración
best_svm_lin_CI <- select_best(svm_lin_res_CI, metric = "bal_accuracy")

# Especificar modelo SVM con kernel radial
svm_rad <- svm_rbf(
  cost = tune(),
  rbf_sigma = tune()
) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# Crear workflow
svm_rad_wf_CI <- workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(svm_rad)

# Definir grilla de hiperparámetros
svm_rad_grid <- grid_regular(
  cost(range = c(-2, 2)),
  rbf_sigma(range = c(-3, -1)),
  levels = 6
)

# Optimizar hiperparámetros
set.seed(123)
svm_rad_res_CI <- tune_grid(
  svm_rad_wf_CI,
  resamples = cv_ci,
  grid = svm_rad_grid,
  metrics = metricas,
  control = ctrl
)

# Seleccionar la mejor configuración
best_svm_rad_CI <- select_best(svm_rad_res_CI, metric = "bal_accuracy")

# Especificar modelo XGBoost
xgb_spec <- boost_tree(
  trees = 800,
  learn_rate = tune(),
  mtry = tune(),
  tree_depth = tune(),
  min_n = tune(),
  loss_reduction = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# Crear workflow
xgb_wf_CI <- workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(xgb_spec)

# Definir grilla de hiperparámetros
xgb_grid_CI <- grid_space_filling(
  learn_rate(range = c(-3, -1)),
  mtry(range = c(1, 4)),
  tree_depth(range = c(2L, 6L)),
  min_n(range = c(2L, 15L)),
  loss_reduction(range = c(-5, -1)),
  size = 25
)

# Optimizar hiperparámetros
set.seed(123)
xgb_res_CI <- tune_grid(
  xgb_wf_CI,
  resamples = cv_ci,
  grid = xgb_grid_CI,
  metrics = metricas,
  control = ctrl
)

# Seleccionar la mejor configuración
best_xgb_CI <- select_best(xgb_res_CI, metric = "bal_accuracy")

# Recolección robusta de métricas para CI
comparacion_modelos_CI <- bind_rows(
  extraer_simple(log_res_CI, "Logistic"),
  extraer_tuneado(enet_res_CI, "Elastic Net"),
  extraer_simple(lda_res_CI, "LDA"),
  extraer_simple(nb_res_CI, "Naive Bayes"),
  extraer_tuneado(rf_res_CI, "Random Forest"),
  extraer_tuneado(svm_lin_res_CI, "SVM Linear"),
  extraer_tuneado(svm_rad_res_CI, "SVM Radial"),
  extraer_tuneado(xgb_res_CI, "XGBoost")
) %>%
  dplyr::arrange(desc(mean_BA))

# Finalizar workflow con los mejores hiperparámetros
best_svm_rad_CI <- select_best(svm_rad_res_CI, metric = "bal_accuracy")
svm_rad_final_wf_CI <- finalize_workflow(svm_rad_wf_CI, best_svm_rad_CI)
svm_rad_fit_CI <- fit(svm_rad_final_wf_CI, data = ci_train)

# Predecir con el modelo final
svm_pred_CI <- predict(svm_rad_fit_CI, ci_test, type = "prob") %>%
  bind_cols(predict(svm_rad_fit_CI, ci_test)) %>%
  bind_cols(dplyr::select(ci_test, Sexo))
svm_pred_CI$Sexo <- as.factor(svm_pred_CI$Sexo)

# Calcular métricas
metricas_clase(svm_pred_CI, truth = Sexo, estimate = .pred_class)
roc_auc(svm_pred_CI, truth = Sexo, .pred_F)

# Bootstrap para intervalos de confianza
boot_res_CI <- boot(
  data = svm_pred_CI,
  statistic = boot_metricas,
  R = 2000
)

# Preparar datos de prueba
ci_test_imp <- ci_test %>%
  dplyr::select(Sexo, MD, BL, MDCu, BLCu)

# Calcular importancia por permutación
set.seed(123)
imp_ci <- vip::vi(
  svm_rad_fit_CI,
  method = "permute",
  train = ci_test_imp,
  target = "Sexo",
  metric = "bal_accuracy",
  pred_wrapper = function(object, newdata) {
    predict(object, newdata)$.pred_class
  },
  nsim = 50
)

# Resumir y graficar importancia
imp_ci_summary <- imp_ci %>%
  group_by(Variable) %>%
  summarise(
    mean_imp = mean(Importance),
    sd_imp   = sd(Importance)
  )

ggplot(imp_ci_summary, aes(x = reorder(Variable, mean_imp), y = mean_imp)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_imp - sd_imp, ymax = mean_imp + sd_imp), width = 0.2) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = "",
    y = "Permutation Importance (Δ Balanced Accuracy)",
    title = "Variable Importance – SVM Radial (CI)"
  )

# Crear grilla para visualizar la superficie de decisión
grid_data <- expand.grid(
  MDCu = seq(min(ci_test$MDCu, na.rm = TRUE),
             max(ci_test$MDCu, na.rm = TRUE),
             length.out = 100),
  MD = seq(min(ci_test$MD, na.rm = TRUE),
           max(ci_test$MD, na.rm = TRUE),
           length.out = 100)
)

# Fijar otras variables en la media
grid_data$BLCu <- mean(ci_test$BLCu, na.rm = TRUE)
grid_data$BL <- mean(ci_test$BL, na.rm = TRUE)

# Predecir en la grilla
grid_pred <- predict(svm_rad_fit_CI, grid_data, type = "prob")
grid_plot <- cbind(grid_data, grid_pred)

# Filtrar datos limpios para puntos
ci_test_clean <- ci_test %>%
  filter(!is.na(MDCu), !is.na(MD))

# Graficar superficie de decisión
ggplot() +
  geom_raster(
    data = grid_plot,
    aes(x = MDCu, y = MD, fill = .pred_F)
  ) +
  scale_fill_viridis_c() +
  geom_point(
    data = ci_test_clean,
    aes(x = MDCu, y = MD, shape = Sexo),
    size = 3,
    color = "black"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Superficie de decisión – SVM Radial (CI)",
    x = "MDCu",
    y = "MD",
    fill = "P(F)"
  )

# Guardar modelos finales
saveRDS(xgb_fit, "xgb_fit_CS.rds")
saveRDS(svm_rad_fit_CI, "svm_rad_fit_CI.rds")

# Guardar conjuntos de prueba
write.table(cs_test, "Conjunto_testeo.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ci_test, "Conjunto_testeo_CI.txt", sep = "\t", row.names = FALSE, quote = FALSE)




# === INFORMACIÓN DE LA SESIÓN ===

session_info <- sessionInfo()
cat("\n\n", "========================================\n",
    "INFORMACIÓN DE SESIÓN - ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "========================================\n",
    "Versión de R: ", R.version.string, "\n",
    "Sistema operativo: ", session_info$running, "\n",
    "Plataforma: ", R.version$platform, "\n",
    "CPU: ", R.version$arch, "\n\n",
    "Paquetes cargados:\n", 
    paste(names(session_info$otherPkgs), sessionInfo()$otherPkgs$ver, 
          sep=" v", collapse="\n  - "),
    "\n\n", "========================================\n\n",
    file = "session_info.txt", append = FALSE)

#Tiempo de ejecucion: 1843.03 segundos



