# 1: Configuración Inicial
# Este chunk carga las bibliotecas necesarias para el análisis.
# Se utiliza tidymodels_prefer() para establecer preferencias de nombres de funciones comunes.

# Paquetes necesarios
required_packages <- c(
  "tidymodels", "ggplot2", "dplyr", "factoextra", "glmnet", "broom", "rsample",
  "discrim", "vip", "yardstick", "patchwork", "workflows", "shapviz", "tidyr",
  "klaR", "ranger", "kernlab", "ggbeeswarm", "boot"
)

# Cargar paquetes
invisible(lapply(required_packages, library, character.only = TRUE))
tidymodels_prefer()

# 2: Carga y Exploración Inicial de Datos
# Carga el archivo de datos 'Dientes.txt' y realiza una primera inspección.
# Se convierten las variables categóricas 'Ind' y 'Sexo' a factores.

# Verificar existencia del archivo
if (!file.exists("Dientes.txt")) {
  stop("Archivo 'Dientes.txt' no encontrado en el directorio de trabajo.")
}

# Cargar datos
Dent <- read.table("Dientes.txt", header = TRUE)

# Mostrar primeras filas y resumen estadístico
print(head(Dent, 20))
print(summary(Dent))

# Estructura del dataset
str(Dent)

# Convertir variables categóricas a factores
Dent$Ind<- as.factor(Dent$Ind)
Dent$Sexo <- as.factor(Dent$Sexo)

# 3: Filtrado y Preparación de Datos para Caninos Superiores (CS)
# Selecciona los registros correspondientes a Caninos Superiores (CS) y guarda un subconjunto.

# Filtrar para Caninos Superiores (CS)
CS <- Dent %>%
  dplyr::filter(Tipo_Diente == "CS")

# Inspección de CS
print(summary(CS[, c("MDCo", "BLCo", "MDCu", "BLCu")]))
print(table(CS$Sexo))
print(length(unique(CS$Ind)))


# Guardar archivo de CS para análisis posteriores
write.table(
  CS,
  file = "caninos_superiores.txt",
  sep = "\t",
  row.names = FALSE
)

# 4: Gráfico de Variables para CS y CI
# Crea un boxplot para visualizar las métricas dentales por sexo y tipo de diente (CS y CI).

# Preparar datos para el gráfico (CS y CI)
vars_plot <- c("MDCo", "BLCo", "MDCu", "BLCu")

# Verificar que las variables existan
missing_vars_plot <- setdiff(vars_plot, names(Dent))
if(length(missing_vars_plot) > 0) stop("Variables faltantes para el gráfico: ", paste(missing_vars_plot, collapse = ", "))

datos_long_plot <- Dent %>%
  dplyr::filter(Tipo_Diente %in% c("CS", "CI")) %>%
  tidyr::pivot_longer(
    cols = all_of(vars_plot),
    names_to = "Variable",
    values_to = "Valor",
    values_drop_na = TRUE
  )

# Crear y guardar el gráfico
p1 <- ggplot(datos_long_plot, aes(x = Sexo, y = Valor, fill = Sexo)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.3) +
  facet_grid(Tipo_Diente ~ Variable, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Sexo",
    y = "Tamaño (mm)",
    title = "Métricas dentales para CS y CI"
  ) +
  theme(legend.position = "none")

ggsave("boxplot_CS_CI.jpg", plot = p1, width = 20, height = 15, units = "cm", dpi = 300)

# 5: Cálculo del Índice de Garn y Bootstrap
# Calcula el índice de dimorfismo sexual (Garn) y sus intervalos de confianza mediante bootstrap.

# Calcular Garn index para CS y CI
vars_garn <- c("MDCo", "BLCo", "MDCu", "BLCu")

# Verificar que las variables existan
missing_vars_garn <- setdiff(vars_garn, names(Dent))
if(length(missing_vars_garn) > 0) stop("Variables faltantes para Garn: ", paste(missing_vars_garn, collapse = ", "))

datos_long_garn <- Dent %>%
  dplyr::filter(Tipo_Diente %in% c("CS", "CI")) %>%
  tidyr::pivot_longer(
    cols = all_of(vars_garn),
    names_to = "Variable",
    values_to = "Valor",
    values_drop_na = TRUE
  )

garn <- datos_long_garn %>%
  group_by(Tipo_Diente, Variable, Sexo) %>%
  summarise(media = mean(Valor, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sexo, values_from = media) %>%
  mutate(Garn = (M - F) / F * 100)

# Función de bootstrap para intervalos de confianza del índice Garn
# Corregido: la función ahora usa 'grupo_datos' en lugar de 'CS' global.
garn_boot <- function(grupo_datos, nboot = 5000) {
  # grupo_datos debe tener columnas 'Valor' y 'Sexo'
  F_vals <- grupo_datos$Valor[grupo_datos$Sexo == "F"]
  M_vals <- grupo_datos$Valor[grupo_datos$Sexo == "M"]

  if(length(F_vals) < 2 || length(M_vals) < 2) {
    warning("No hay suficientes datos para ambos sexos en este grupo.")
    return(c(SD = NA, CI_low = NA, CI_high = NA))
  }

  res <- replicate(nboot, {
    F_sample <- sample(F_vals, replace = TRUE)
    M_sample <- sample(M_vals, replace = TRUE)

    F_mean <- mean(F_sample, na.rm = TRUE)
    M_mean <- mean(M_sample, na.rm = TRUE)

    if(F_mean == 0) { # Evitar división por cero
      return(NA)
    }
    (M_mean - F_mean) / F_mean * 100
  })

  res <- res[!is.na(res)] # Remover NAs generados por división por cero
  if(length(res) == 0) return(c(SD = NA, CI_low = NA, CI_high = NA))

  c(
    SD = mean(res),
    CI_low = quantile(res, 0.025, names = FALSE),
    CI_high = quantile(res, 0.975, names = FALSE)
  )
}

# Aplicar bootstrap a cada combinación de Tipo_Diente y Variable
garn_results <- datos_long_garn %>%
  group_by(Tipo_Diente, Variable) %>%
  group_modify(~ {
    boot_result <- garn_boot(.x) # .x es el subconjunto del grupo actual
    data.frame(
      SD = boot_result[1],
      CI_low = boot_result[2],
      CI_high = boot_result[3]
    )
  }) %>%
  ungroup()

# Unir resultados Garn con intervalos
garn_final <- left_join(garn, garn_results, by = c("Tipo_Diente", "Variable"))

# Gráfico CS-CI con Garn index
p_base <- ggplot(datos_long_garn, aes(x = Sexo, y = Valor, fill = Sexo)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
  facet_grid(Tipo_Diente ~ Variable, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Sexo",
    y = "Tamaño (mm)"
  ) +
  theme(legend.position = "none")

# Agregar texto de Garn y CI
p_with_labels <- p_base +
  geom_text(
    data = garn_final,
    aes(
      x = 1.5, # Centrado entre las categorias
      y = Inf,
      label = paste0("SDI = ", round(SD, 1), "%\n", "95% CI [", round(CI_low, 1), ", ", round(CI_high, 1), "]")
    ),
    vjust = 1.3,
    hjust = 0.5, # Centrado horizontal
    inherit.aes = FALSE,
    size = 2.5
  )

print(p_with_labels)
ggsave("medidas_CS_CI+Garn.jpg", plot = p_with_labels, width = 20, height = 15, units = "cm", dpi = 300)

# 6: Gráfico de Intervalos de Confianza del Garn
# Visualiza los índices de Garn con sus intervalos de confianza bootstrap.

# Preparar datos para el gráfico de intervalos
garn_for_plot <- garn_final %>%
  mutate(
    Variable = factor(Variable, levels = c("MDCo", "BLCo", "MDCu", "BLCu"))
  )

#graficar Garn con intervalos
p_dim <- ggplot(garn_for_plot,
       aes(x = SD, y = Variable, color = Tipo_Diente)) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey60",
             linewidth = 0.5) +
  
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    width = 0.15,
    linewidth = 1,
    position = position_dodge(width = 0 )  +
  
  geom_point(
    size = 3,
    position = position_dodge(width = 0.5)
  ) +
  
  theme_classic(base_size = 13) +
  
  labs(
    x = "Sexual Dimorphism Index (%)",
    y = "Dental measurement",
    color = "Tooth",
    title = "Garn index with bootstrap 95% confidence intervals"
  ) +
  
  scale_color_manual(
    values = c(
      "CS" = "#0072B2",
      "CI" = "#D55E00"
    )
  )

print(p_dim)
ggsave("Intervalo_confianza_Garn.jpg", plot = p_dim, width = 20, height = 15, units = "cm", dpi = 300)

# 7: Función de División de Datos
# Define una función para dividir los datos manteniendo la proporción de sexo y evitando duplicación de individuos entre entrenamiento y testeo.

# Función para dividir datos estratificando por Sexo y agrupando por Individuo
split_dentales <- function(data, prop = 0.8) {
  required_cols <- c("Ind", "Sexo")
  missing <- setdiff(required_cols, names(data))
  if(length(missing) > 0){
    stop(paste("Faltan columnas requeridas:", 
               paste(missing, collapse = ", ")))
  }
  
  set.seed(123) # Semilla para reproducibilidad
  
  split_obj <- group_initial_split(
    data,
    prop = prop,
    group = Ind,     # Agrupa por individuo
    strata = Sexo    # Estratifica por sexo
  )
  
  list(
    train = training(split_obj),
    test  = testing(split_obj),
    split = split_obj
  )
}

# 8: División de Datos para CS
# Aplica la función de división al dataset de Caninos Superiores (CS).

# Dividir la muestra CS
cs_split <- split_dentales(CS)

cs_train <- cs_split$train # Entrenamiento
cs_test <- cs_split$test   # Testeo

# Verificar divisiones
print(summary(cs_train))
print(summary(cs_test))

# Chequeo de individuos repetidos
print(paste("Intersección de individuos entre train y test:", length(intersect(unique(cs_train$Ind), unique(cs_test$Ind)))))

# 9: Limpieza y Preparación de Receta para CS
# Imputa valores faltantes en el conjunto de entrenamiento y define la receta de preprocesamiento.

# Remover fila completamente NA en train (si existe)
n_rows_before <- nrow(cs_train)
cs_train <- cs_train %>%
  dplyr::filter(!(is.na(MDCo) & is.na(BLCo) & is.na(MDCu) & is.na(BLCu)))
n_rows_after <- nrow(cs_train)
if(n_rows_before != n_rows_after) cat("Eliminada 1 fila con todos los predictores NA en cs_train.\n")

# Chequeo post-limpieza
check_na_all_na <- cs_train %>%
  dplyr::filter(is.na(MDCo) & is.na(BLCo) & is.na(MDCu) & is.na(BLCu))
if(nrow(check_na_all_na) > 0) cat("Advertencia: Aún hay filas con todos los predictores NA en cs_train.\n")

# Definir receta de preprocesamiento para CS
# La receta se ajusta solo al conjunto de entrenamiento
cs_recipe <- recipe(Sexo ~ MDCo + BLCo + MDCu + BLCu, data = cs_train) %>%
  step_impute_knn(all_predictors()) %>% # Imputación KNN
  step_normalize(all_predictors())      # Estandarización

# Preparar la receta (esto ajusta los parámetros de imputación y normalización)
cs_recipe_prep <- prep(cs_recipe)
print(juice(cs_recipe_prep)) # Ver datos procesados

# 10: Validación Cruzada para CS
# Configura la validación cruzada estratificada y agrupada para evitar pseudorreplicación.

# Configurar validación cruzada agrupada por Ind y estratificada por Sexo
set.seed(123)
cv_cs <- group_vfold_cv(
  cs_train,
  v = 10,       # 10 folds
  repeats = 5,  # 5 repeticiones
  group = Ind,  # Agrupar por individuo
  strata = Sexo # Estratificar por sexo
)

# Verificar estructura de los folds
print(cv_cs)

# 11: Definición de Métricas y Control
# Define las métricas de evaluación y el control para los remuestreos.

# Métricas de evaluación
metricas <- metric_set(
  accuracy,
  sens,
  spec,
  bal_accuracy,
  roc_auc
)

# Control de remuestreo
ctrl <- control_resamples(save_pred = TRUE)

# Workflow base para CS
cs_wf <- workflow() %>%
  add_recipe(cs_recipe)

# 12: Definición y Evaluación de Modelos para CS
# Define y evalúa múltiples modelos de clasificación en el conjunto de validación cruzada de CS.

# --- Modelo 1: Regresión Logística ---
log_spec <- logistic_reg(mode = "classification") %>%
  set_engine("glm")

log_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(log_spec)

log_res <- fit_resamples(
  log_wf,
  resamples = cv_cs,
  control = ctrl,
  metrics = metricas
)
print(collect_metrics(log_res))

# --- Modelo 2: Elastic Net ---
enet_spec <- logistic_reg(
  penalty = tune(), # Parámetro lambda para penalización
  mixture = tune()  # Parámetro alpha (0=Ridge, 1=Lasso)
) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

enet_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(enet_spec)

enet_grid <- grid_regular(
  penalty(range = c(-4, 0)),  # Rango en escala logarítmica
  mixture(range = c(0, 1)),   # Valores entre Ridge y Lasso
  levels = 10                 # 10x10 = 100 combinaciones
)

set.seed(123)
enet_res <- tune_grid(
  enet_wf,
  resamples = cv_cs,
  grid = enet_grid,
  control = ctrl,
  metrics = metricas
)
print(collect_metrics(enet_res))

# --- Modelo 3: LDA ---
lda_spec <- discrim_linear() %>%
  set_engine("MASS")

lda_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(lda_spec)

lda_res <- fit_resamples(
  lda_wf,
  resamples = cv_cs,
  metrics = metricas,
  control = ctrl
)
print(collect_metrics(lda_res))

# --- Modelo 4: Naive Bayes ---
nb_spec <- naive_Bayes() %>%
  set_engine("klaR")

nb_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(nb_spec)

nb_res <- fit_resamples(
  nb_wf,
  resamples = cv_cs,
  metrics = metricas,
  control = ctrl
)
print(collect_metrics(nb_res))

# --- Modelo 5: Random Forest ---
rf_spec <- rand_forest(
  trees = 1000,      # Número de árboles
  mtry = tune(),     # Número de variables para splits
  min_n = tune()     # Mínimo número de observaciones en nodos hoja
) %>%
  set_engine("ranger") %>%
  set_mode("classification")

rf_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(rf_spec)

rf_grid <- grid_regular(
  mtry(range = c(1, 4)),
  min_n(range = c(2, 20)),
  levels = 6 # 4x6 = 24 combinaciones
)

set.seed(123)
rf_res <- tune_grid(
  rf_wf,
  resamples = cv_cs,
  grid = rf_grid,
  metrics = metricas,
  control = ctrl
)
print(collect_metrics(rf_res))

# --- Modelo 6: SVM Lineal ---
svm_lin_spec <- svm_linear(cost = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

svm_lin_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(svm_lin_spec)

svm_lin_grid <- grid_regular(
  cost(range = c(-3, 2)), # Rango del parámetro de penalización C
  levels = 10
)

set.seed(123)
svm_lin_res <- tune_grid(
  svm_lin_wf,
  resamples = cv_cs,
  grid = svm_lin_grid,
  metrics = metricas,
  control = ctrl
)
print(collect_metrics(svm_lin_res))

# --- Modelo 7: SVM Radial ---
svm_rad_spec <- svm_rbf(
  cost = tune(),
  rbf_sigma = tune()
) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

svm_rad_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(svm_rad_spec)

svm_rad_grid <- grid_regular(
  cost(range = c(-2, 2)),
  rbf_sigma(range = c(-3, -1)),
  levels = 6
)

set.seed(123)
svm_rad_res <- tune_grid(
  svm_rad_wf,
  resamples = cv_cs,
  grid = svm_rad_grid,
  metrics = metricas,
  control = ctrl
)
print(collect_metrics(svm_rad_res))

# --- Modelo 8: XGBoost ---
xgb_spec <- boost_tree(
  trees = 800,              # Número de árboles
  learn_rate = tune(),      # Tasa de aprendizaje
  mtry = tune(),            # Variables por split
  tree_depth = tune(),      # Profundidad máxima de los árboles
  min_n = tune(),           # Mínimo nodos hoja
  loss_reduction = tune()   # Reducción mínima de pérdida
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

xgb_wf <- workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(xgb_spec)

xgb_grid <- grid_space_filling(
  learn_rate(range = c(-3, -1)),      # 0.001 a 0.1
  mtry(range = c(1, 4)),
  tree_depth(range = c(2L, 6L)),
  min_n(range = c(2L, 15L)),
  loss_reduction(range = c(-5, -1)),  # 1e-5 a 0.1
  size = 25                           # 25 combinaciones
)

set.seed(123)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = cv_cs,
  grid = xgb_grid,
  metrics = metricas,
  control = ctrl
)
print(collect_metrics(xgb_res))

# 13: Comparación de Modelos para CS
# Resume las métricas de desempeño de todos los modelos de CS y determina el mejor.

# Funciones auxiliares para extraer métricas de los resultados de los modelos
extraer_simple <- function(res_obj, nombre_modelo) {
  collect_metrics(res_obj, summarize = FALSE) %>%
    dplyr::filter(.metric == "bal_accuracy") %>%
    dplyr::summarise(
      mean_BA = mean(.estimate),
      sd_BA   = sd(.estimate)
    ) %>%
    dplyr::mutate(modelo = nombre_modelo)
}

extraer_tuneado <- function(res_obj, nombre_modelo) {
  best_id <- select_best(res_obj, metric = "bal_accuracy")$.config
  collect_metrics(res_obj, summarize = FALSE) %>%
    dplyr::filter(.metric == "bal_accuracy", .config == best_id) %>%
    dplyr::summarise(
      mean_BA = mean(.estimate),
      sd_BA   = sd(.estimate)
    ) %>%
    dplyr::mutate(modelo = nombre_modelo)
}

# Recopilar comparación
comparacion_modelos_CS <-
  bind_rows(
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

print(comparacion_modelos_CS)

# Seleccionar el mejor modelo basado en mean_BA (XGBoost suele ganar según el script original)
best_model_name_CS <- comparacion_modelos_CS$modelo[1]
cat("El mejor modelo para CS es:", best_model_name_CS, "\n")

# 14: Entrenamiento Final y Predicción del Mejor Modelo para CS
# Entrena el modelo final de CS con todo el conjunto de entrenamiento y lo evalúa en el de testeo.

# Basado en el script original, se asume que XGBoost fue el mejor (índice 8)
# Si el orden cambia, este paso debe ajustarse.
if(best_model_name_CS == "XGBoost") {
  best_config_xgb <- select_best(xgb_res, metric = "bal_accuracy")
  xgb_final_wf <- finalize_workflow(xgb_wf, best_config_xgb)
  xgb_fit <- fit(xgb_final_wf, data = cs_train)
  
  # Predicciones en test
  xgb_pred <- predict(xgb_fit, cs_test, type = "prob") %>%
    bind_cols(predict(xgb_fit, cs_test)) %>%
    bind_cols(dplyr::select(cs_test, Sexo))
  
  xgb_pred$Sexo <- as.factor(xgb_pred$Sexo)
  
  # Métricas finales en test
  final_metrics_cs <- metricas(xgb_pred, truth = Sexo, estimate = .pred_class)
  print(final_metrics_cs)
  
  # AUC
  auc_value_cs <- roc_auc(xgb_pred, truth = Sexo, .pred_F)
  print(auc_value_cs)
  
  # Curva ROC
  roc_plot_cs <- roc_curve(xgb_pred, truth = Sexo, .pred_F) %>%
    autoplot()
  print(roc_plot_cs)
  ggsave("ROC_CS.jpg", plot = roc_plot_cs, width = 20, height = 15, units = "cm", dpi = 300)
  
} else {
  stop("El mejor modelo para CS no es XGBoost. Este script está configurado para continuar con XGBoost. Ajuste manual requerido.")
}

# 15: Bootstrap de Intervalos de Confianza en Test para CS
# Calcula intervalos de confianza bootstrap para las métricas del modelo final de CS en el conjunto de testeo.

library(boot)

boot_metricas <- function(data, indices) {
  d <- data[indices, ]
  # Asegurar que .pred_class sea factor con los mismos niveles que Sexo
  d$.pred_class <- factor(d$.pred_class, levels = levels(d$Sexo))
  c(
    bal_acc = bal_accuracy(d, truth = Sexo, estimate = .pred_class)$.estimate,
    sens    = sens(d, truth = Sexo, estimate = .pred_class)$.estimate,
    spec    = spec(d, truth = Sexo, estimate = .pred_class)$.estimate,
    acc     = accuracy(d, truth = Sexo, estimate = .pred_class)$.estimate
  )
}

set.seed(123)
boot_res_cs <- boot(
  data = xgb_pred,
  statistic = boot_metricas,
  R = 2000
)

# Intervalos de confianza percentiles 95%
cat("\n--- Intervalos de Confianza Bootstrap (95%) para CS ---\n")
print(boot.ci(boot_res_cs, type = "perc", index = 1)) # Balanced Acc
print(boot.ci(boot_res_cs, type = "perc", index = 2)) # Sensibilidad
print(boot.ci(boot_res_cs, type = "perc", index = 3)) # Especificidad
print(boot.ci(boot_res_cs, type = "perc", index = 4)) # Accuracy

# 16: Análisis de Umbral y Matriz de Confusión para CS
# Analiza el umbral óptimo y genera la matriz de confusión para el modelo final de CS.

# Gráfico de Sensibilidad y Especificidad vs Umbral
threshold_df_cs <- roc_curve(xgb_pred, truth = Sexo, .pred_F) %>%
  dplyr::select(.threshold, sensitivity, specificity)

p_thresh_cs <- ggplot(threshold_df_cs, aes(x = .threshold)) +
  geom_line(aes(y = sensitivity, color = "Sensibilidad"), size = 1) +
  geom_line(aes(y = specificity, color = "Especificidad"), size = 1) +
  scale_color_manual(values = c("Sensibilidad" = "red", "Especificidad" = "blue")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Umbral",
    y = "Valor",
    color = " ",
    title = "Sensibilidad y Especificidad según Umbral (XGBoost - CS)"
  )
print(p_thresh_cs)

# Umbral más balanceado (dif mínima entre sens y spec)
best_threshold_cs <- threshold_df_cs %>%
  dplyr::mutate(diff = abs(sensitivity - specificity)) %>%
  dplyr::arrange(diff) %>%
  dplyr::slice(1)
print(best_threshold_cs)

# Matriz de confusión
conf_mat_cs <- conf_mat(xgb_pred, truth = Sexo, estimate = .pred_class)
print(conf_mat_cs)

# 17: Importancia de Variables para CS
# Calcula y grafica la importancia de las variables para el modelo final de CS.

# Importancia de variables (VIP)
xgb_importance <- xgb_fit %>%
  extract_fit_parsnip() %>%
  vip::vi()

# Gráfico de importancia
p_vip_cs <- vip(xgb_fit, num_features = 10, geom = "col", aesthetics = list(fill = "steelblue"))
print(p_vip_cs)
ggsave("Var_Imp_CS.jpg", plot = p_vip_cs, width = 20, height = 15, units = "cm", dpi = 300)

# 18: Análisis SHAP para CS
# Calcula y visualiza las contribuciones SHAP de las variables para el modelo final de CS.

# Extraer modelo core
xgb_core <- extract_fit_parsnip(xgb_fit)$fit

# Obtener datos procesados por la receta
rec <- extract_recipe(xgb_fit)
X_df <- bake(rec, new_data = cs_train) %>%
  dplyr::select(-Sexo)

X_matrix <- as.matrix(X_df)

# Calcular SHAP
sv <- shapviz(
  xgb_core,
  X = X_matrix,
  X_pred = X_matrix,
  baseline = "auto"
)

# Gráfico de importancia SHAP
p_shap_imp_cs <- sv_importance(sv, kind = "beeswarm")
print(p_shap_imp_cs)
ggsave("Shap_CS.jpg", plot = p_shap_imp_cs, width = 20, height = 15, units = "cm", dpi = 300)

# Gráficos de dependencia SHAP
vars_shap <- c("MDCu", "BLCu", "BLCo", "MDCo")
plots_shap <- lapply(vars_shap, function(v) {
  sv_dependence(sv, v) +
    ggtitle(paste("Dependencia SHAP -", v)) +
    theme_minimal(base_size = 12)
})
p_shap_dep_cs <- wrap_plots(plots_shap, ncol = 2)
print(p_shap_dep_cs)

# 19: Preparación de Datos para Caninos Inferiores (CI)
# Repite pasos similares para el dataset de Caninos Inferiores (CI).

CI <- Dent %>%
  dplyr::filter(Tipo_Diente == "CI")

# Inspección de CI
print(summary(CI[, c("MDCo", "BLCo", "MDCu", "BLCu")]))
print(table(CI$Sexo))
print(length(unique(CI$Ind)))
cat("Dataset CI balanceado y con pocos datos faltantes\n")

# Guardar archivo de CI para análisis posteriores
write.table(
  CI,
  file = "caninos_inferiores.txt",
  sep = "\t",
  row.names = FALSE
)

# 20: División y Preparación de Datos para CI
# Divide los datos de CI y define la receta de preprocesamiento.

# Dividir la muestra CI
ci_split <- split_dentales(CI)
ci_train <- ci_split$train
ci_test <- ci_split$test

# Verificar divisiones
print(summary(ci_train))
print(summary(ci_test))
print(paste("Intersección de individuos entre train y test (CI):", length(intersect(unique(ci_train$Ind), unique(ci_test$Ind)))))

# Definir receta de preprocesamiento para CI
ci_recipe <- recipe(Sexo ~ MDCo + BLCo + MDCu + BLCu, data = ci_train) %>%
  step_impute_knn(all_predictors()) %>%
  step_normalize(all_predictors())

# Preparar la receta
ci_recipe_prep <- prep(ci_recipe)

# Configurar validación cruzada para CI
set.seed(123)
cv_ci <- group_vfold_cv(
  ci_train,
  v = 10,
  repeats = 5,
  group = Ind,
  strata = Sexo
)

# 21: Evaluación de Modelos para CI
# Evalúa los mismos modelos utilizados para CS en el conjunto de validación cruzada de CI.

# Reutilizamos las mismas especificaciones de modelos, cambiando solo el workflow y los datos.
# Regresión Logística CI
log_res_CI <- fit_resamples(
  log_wf, # Reutiliza el workflow de CS con receta de CI
  resamples = cv_ci,
  control = ctrl,
  metrics = metricas
)
print(collect_metrics(log_res_CI))

# LDA CI
lda_res_CI <- fit_resamples(
  lda_wf, # Reutiliza el workflow de CS con receta de CI
  resamples = cv_ci,
  control = ctrl,
  metrics = metricas
)
print(collect_metrics(lda_res_CI))

# Naive Bayes CI
nb_res_CI <- fit_resamples(
  nb_wf, # Reutiliza el workflow de CS con receta de CI
  resamples = cv_ci,
  control = ctrl,
  metrics = metricas
)
print(collect_metrics(nb_res_CI))

# Random Forest CI
rf_res_CI <- tune_grid(
  rf_wf, # Reutiliza el workflow de CS con receta de CI
  resamples = cv_ci,
  grid = rf_grid,
  metrics = metricas,
  control = ctrl
)
print(collect_metrics(rf_res_CI))

# Elastic Net CI
enet_res_CI <- tune_grid(
  enet_wf, # Reutiliza el workflow de CS con receta de CI
  resamples = cv_ci,
  grid = enet_grid,
  control = ctrl,
  metrics = metricas
)
print(collect_metrics(enet_res_CI))

# SVM Lineal CI
svm_lin_res_CI <- tune_grid(
  svm_lin_wf, # Reutiliza el workflow de CS con receta de CI
  resamples = cv_ci,
  grid = svm_lin_grid,
  metrics = metricas,
  control = ctrl
)
print(collect_metrics(svm_lin_res_CI))

# SVM Radial CI
 <- tune_grid(
  svm_rad_wf, # Reutiliza el workflow de CS con receta de CI
  resamples = cv_ci,
  grid = svm_rad_grid,
  metrics = metricas,
  control = ctrl
)
print(collect_metrics(svm_rad_res_CI))

# XGBoost CI
xgb_res_CI <- tune_grid(
  xgb_wf, # Reutiliza el workflow de CS con receta de CI
  resamples = cv_ci,
  grid = xgb_grid,
  metrics = metricas,
  control = ctrl
)
print(collect_metrics(xgb_res_CI))

# 22: Comparación y Selección del Mejor Modelo para CI
# Resume las métricas de los modelos para CI y selecciona el mejor.

# Recopilar comparación para CI
comparacion_modelos_CI <-
  bind_rows(
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

print(comparacion_modelos_CI)

# Seleccionar el mejor modelo basado en mean_BA (SVM Radial suele ganar según el script original)
best_model_name_CI <- comparacion_modelos_CI$modelo[1]
cat("El mejor modelo para CI es:", best_model_name_CI, "\n")

# 23: Entrenamiento Final y Predicción del Mejor Modelo para CI
# Entrena el modelo final de CI y lo evalúa en el de testeo.

# Basado en el script original, se asume que SVM Radial fue el mejor (índice 7)
if(best_model_name_CI == "SVM Radial") {
  best_config_svm_rad_CI <- select_best(svm_rad_res_CI, metric = "bal_accuracy")
  svm_rad_final_wf_CI <- finalize_workflow(svm_rad_wf, best_config_svm_rad_CI)
  svm_rad_fit_CI <- fit(svm_rad_final_wf_CI, data = ci_train)
  
  # Predicciones en test
  svm_pred_CI <- predict(svm_rad_fit_CI, ci_test, type = "prob") %>%
    bind_cols(predict(svm_rad_fit_CI, ci_test)) %>%
    bind_cols(dplyr::select(ci_test, Sexo))
  
  svm_pred_CI$Sexo <- as.factor(svm_pred_CI$Sexo)
  
  # Métricas finales en test
  final_metrics_ci <- metricas(svm_pred_CI, truth = Sexo, estimate = .pred_class)
  print(final_metrics_ci)
  
  # AUC
  auc_value_ci <- roc_auc(svm_pred_CI, truth = Sexo, .pred_F)
  print(auc_value_ci)
  
  # Curva ROC
  roc_plot_ci <- roc_curve(svm_pred_CI, truth = Sexo, .pred_F) %>%
    autoplot()
  print(roc_plot_ci)
  ggsave("ROC_CI.jpg", plot = roc_plot_ci, width = 20, height = 15, units = "cm", dpi = 300)
  
} else {
  stop("El mejor modelo para CI no es SVM Radial. Este script está configurado para continuar con SVM Radial. Ajuste manual requerido.")
}

# 24: Análisis de Umbral, Distribución de Probabilidades y Matriz de Confusión para CI
# Analiza el umbral óptimo, distribución de probabilidades y genera la matriz de confusión para el modelo final de CI.

# Gráfico de Sensibilidad y Especificidad vs Umbral (CI)
threshold_df_CI <- roc_curve(svm_pred_CI, truth = Sexo, .pred_F) %>%
  dplyr::select(.threshold, sensitivity, specificity)

p_thresh_ci <- ggplot(threshold_df_CI, aes(x = .threshold)) +
  geom_line(aes(y = sensitivity, color = "Sensibilidad"), size = 1) +
  geom_line(aes(y = specificity, color = "Especificidad"), size = 1) +
  scale_color_manual(values = c("Sensibilidad" = "red", "Especificidad" = "blue")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Umbral",
    y = "Valor",
    color = " ",
    title = "Sensibilidad y Especificidad según Umbral (SVM Radial - CI)"
  )
print(p_thresh_ci)

# Distribución de probabilidades
p_dist_prob_ci <- ggplot(svm_pred_CI, aes(x = .pred_F, fill = Sexo)) +
  geom_density(alpha = 0.4) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Probabilidad predicha de F",
    y = "Densidad",
    title = "Distribución de Probabilidades (SVM Radial - CI)"
  )
print(p_dist_prob_ci)
ggsave("Dist_prob.jpg", plot = p_dist_prob_ci, width = 20, height = 15, units = "cm", dpi = 300)

# Matriz de confusión
conf_mat_ci <- conf_mat(svm_pred_CI, truth = Sexo, estimate = .pred_class)
print(conf_mat_ci)

# 25: Importancia de Variables para CI
# Calcula la importancia de las variables para el modelo final de CI usando permutación.

# Preparar datos para cálculo de importancia
ci_test_imp <- ci_test %>%
  dplyr::select(Sexo, MDCo, BLCo, MDCu, BLCu)

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

imp_ci_summary <- imp_ci %>%
  group_by(Variable) %>%
  summarise(
    mean_imp = mean(Importance),
    sd_imp   = sd(Importance)
  )

print(imp_ci_summary)

# Gráfico de importancia
p_imp_ci <- ggplot(imp_ci_summary,
       aes(x = reorder(Variable, mean_imp), y = mean_imp)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_imp - sd_imp, ymax = mean_imp + sd_imp), width = 0.2) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = "",
    y = "Permutation Importance (Δ Balanced Accuracy)",
    title = "Variable Importance – SVM Radial (CI)"
  )
print(p_imp_ci)
ggsave("Variable_imp_CI.jpg", plot = p_imp_ci, width = 20, height = 15, units = "cm", dpi = 300)

# 26: Gráfico de Superficie de Decisión para CI
# Visualiza la superficie de decisión del modelo final de CI en un subespacio de variables.

# Crear grilla para la superficie de decisión (ej: MDCu vs MDCo)
grid_data <- expand.grid(
  MDCu = seq(min(ci_test$MDCu, na.rm = TRUE), max(ci_test$MDCu, na.rm = TRUE), length.out = 100),
  MDCo = seq(min(ci_test$MDCo, na.rm = TRUE), max(ci_test$MDCo, na.rm = TRUE), length.out = 100)
)
grid_data$BLCu <- mean(ci_test$BLCu, na.rm = TRUE)
grid_data$BLCo <- mean(ci_test$BLCo, na.rm = TRUE)

# Predecir en la grilla
grid_pred <- predict(svm_rad_fit_CI, grid_data, type = "prob")

grid_plot <- cbind(grid_data, grid_pred)

# Filtrar datos de test para plotear puntos
ci_test_clean <- ci_test %>%
  dplyr::filter(!is.na(MDCu), !is.na(MDCo))

# Plotear superficie y puntos
p_surf_ci <- ggplot() +
  geom_raster(
    data = grid_plot,
    aes(x = MDCu, y = MDCo, fill = .pred_F)
  ) +
  scale_fill_viridis_c() +
  geom_point(
    data = ci_test_clean,
    aes(x = MDCu, y = MDCo, shape = Sexo),
    size = 3,
    color = "black"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Superficie de decisión – SVM Radial (CI)",
    x = "MDCu",
    y = "MDCo",
    fill = "P(F)"
  )
print(p_surf_ci)
ggsave("SVM_sup.jpg", plot = p_surf_ci, width = 20, height = 15, units = "cm", dpi = 300)

# 27: Resumen Final y Guardado de Resultados
# Presenta un resumen comparativo final y guarda los modelos entrenados y los conjuntos de testeo.

# Tabla resumen de resultados finales
results_summary <- data.frame(
  Diente = c("CS", "CI"),
  Modelo_Final = c("XGBoost", "SVM Radial"),
  BA = c(final_metrics_cs %>% filter(.metric == "bal_accuracy") %>% pull(.estimate),
         final_metrics_ci %>% filter(.metric == "bal_accuracy") %>% pull(.estimate)),
  IC95_BA = c(
    paste0(round(boot.ci(boot_res_cs, type = "perc", index = 1)$percent[4], 2), " – ", round(boot.ci(boot_res_cs, type = "perc", index = 1)$percent[5], 2)),
    paste0(round(boot.ci(boot_res_CI, type = "perc", index = 1)$percent[4], 2), " – ", round(boot.ci(boot_res_CI, type = "perc", index = 1)$percent[5], 2))
  ),
  Sens = c(final_metrics_cs %>% filter(.metric == "sens") %>% pull(.estimate),
           final_metrics_ci %>% filter(.metric == "sens") %>% pull(.estimate)),
  IC95_Sens = c(
    paste0(round(boot.ci(boot_res_cs, type = "perc", index = 2)$percent[4], 2), " – ", round(boot.ci(boot_res_cs, type = "perc", index = 2)$percent[5], 2)),
    paste0(round(boot.ci(boot_res_CI, type = "perc", index = 2)$percent[4], 2), " – ", round(boot.ci(boot_res_CI, type = "perc", index = 2)$percent[5], 2))
  ),
  Spec = c(final_metrics_cs %>% filter(.metric == "spec") %>% pull(.estimate),
           final_metrics_ci %>% filter(.metric == "spec") %>% pull(.estimate)),
  IC95_Spec = c(
    paste0(round(boot.ci(boot_res_cs, type = "perc", index = 3)$percent[4], 2), " – ", round(boot.ci(boot_res_cs, type = "perc", index = 3)$percent[5], 2)),
    paste0(round(boot.ci(boot_res_CI, type = "perc", index = 3)$percent[4], 2), " – ", round(boot.ci(boot_resperc", index = 3)$percent[5], 2))
  ),
  AUC = c(auc_value_cs$.estimate, auc_value_ci$.estimate)
)

print(results_summary)

# Guardar modelos entrenados
saveRDS(xgb_fit, "xgb_fit_CS.rds")
saveRDS(svm_rad_fit_CI, "svm_rad_fit_CI.rds")

# Guardar conjuntos de testeo
write.table(cs_test, "Conjunto_testeo_CS.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ci_test, "Conjunto_testeo_CI.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n--- Proceso completado ---\n")
cat("- Modelos guardados: 'xgb_fit_CS.rds', 'svm_rad_fit_CI.rds'\n")
cat("- Conjuntos de testeo guardados: 'Conjunto_testeo_CS.txt', 'Conjunto_testeo_CI.txt'\n")
cat("- Gráficos guardados en el directorio actual.\n")


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
    file = "session_info_ajuste_modelos.txt", append = FALSE)


# Tiempo de ejecución aproximado: ~1800-2000 segundos (30-35 minutos)
