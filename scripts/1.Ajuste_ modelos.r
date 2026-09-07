# ==============================================================================
# SCRIPT INTEGRADO: EVALUACIÓN COMPARATIVA DE MODELOS DE APRENDIZAJE AUTOMÁTICO
#APLICADOS A LA ESTIMACIÓN SEXUAL MEDIANTE EL ANÁLISIS MÉTRICO DE CANINOS
#PERMANENTES.

#Leandro Luna, Claudia Aranda y Marcelo Cardillo

# Descripción: Codigo de entrenamiento y testeo mediante Machine Learning 
#para Canino Superior (CS) y Canino inferior (CI) 
# ==============================================================================

# --- BLOQUE R 001 ----------------------------------------------------------
#paquetes
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
tidymodels_prefer()

# --- BLOQUE R 002 ----------------------------------------------------------
Dent=read.table("Dientes.txt", header=TRUE)

head(Dent,20)

# --- BLOQUE R 003 ----------------------------------------------------------
summary(Dent)

str(Dent)

# --- BLOQUE R 004 ----------------------------------------------------------
Dent$Ind=as.factor(Dent$Ind)
Dent$Sexo=as.factor(Dent$Sexo)

# --- BLOQUE R 005 ----------------------------------------------------------
###arrancamos con caninos superiores (CS)


CS <- Dent %>%
  filter(Tipo_Diente == "CS")

#chequeando
#cuántos dientes,
#distribución por sexo,
#n individuos

summary(CS[,c("MDCo","BLCo","MDCu","BLCu")])
table(CS$Sexo)
length(unique(CS$Ind))



# --- BLOQUE R 006 ----------------------------------------------------------
##archivo de CS para análisis posteriores

write.table(
  CS,
  file = "caninos_superiores.txt",
  sep = "\t",
  row.names = FALSE)


# --- BLOQUE R 007 ----------------------------------------------------------
#grafico variables para CS-CI (todos menos los faltantes)

vars <- c("MDCo","BLCo","MDCu","BLCu")

datos_long <- Dent%>% 
  filter(Tipo_Diente %in% c("CS","CI")) %>% 
  pivot_longer(
    cols = all_of(vars),
    names_to = "Variable",
    values_to = "Valor",
    values_drop_na = TRUE
  )

ggplot(datos_long, aes(x = Sexo, y = Valor, fill = Sexo)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.3) +
  facet_grid(Tipo_Diente ~ Variable, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Sexo",
    y = "Tamaño (mm)",
    title = "Métricas dentales"
  ) +
  theme(legend.position = "none")

# --- BLOQUE R 008 ----------------------------------------------------------
ggsave("boxplot_CS_CI.jpg", width = 20, height = 15, units = "cm", dpi = 300)

# --- BLOQUE R 009 ----------------------------------------------------------
# Comparación entre sexos: Mann-Whitney U con ajuste de Holm-Bonferroni


comparacion_sexos <- datos_long %>%
  group_by(Tipo_Diente, Variable) %>%
  summarise(
    n_M = sum(Sexo == "M"),
    n_F = sum(Sexo == "F"),
    mediana_M = median(Valor[Sexo == "M"], na.rm = TRUE),
    mediana_F = median(Valor[Sexo == "F"], na.rm = TRUE),
    p = wilcox.test(Valor ~ Sexo, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_ajustado = p.adjust(p, method = "holm")
  )

comparacion_sexos

# --- BLOQUE R 010 ----------------------------------------------------------
## computar Garn index
vars <- c("MDCo","BLCo","MDCu","BLCu")

garn <- Dent %>% 
  filter(Tipo_Diente %in% c("CS","CI")) %>% 
  pivot_longer(
    cols = all_of(vars),
    names_to = "Variable",
    values_to = "Valor",
    values_drop_na = TRUE
  ) %>% 
  group_by(Tipo_Diente, Variable, Sexo) %>% 
  summarise(media = mean(Valor), .groups = "drop") %>% 
  pivot_wider(names_from = Sexo, values_from = media) %>% 
  mutate(Garn = (M - F) / F * 100)

# --- BLOQUE R 011 ----------------------------------------------------------
##bootstrap para CI
garn_boot <- function(Dent, nboot = 5000){
  
  res <- replicate(nboot, {
    
    F_sample <- sample(Dent$Valor[Dent$Sexo == "F"], replace = TRUE)
    M_sample <- sample(Dent$Valor[Dent$Sexo == "M"], replace = TRUE)
    
    F_mean <- mean(F_sample, na.rm = TRUE)
    M_mean <- mean(M_sample, na.rm = TRUE)
    
    (M_mean - F_mean) / F_mean * 100
    
  })
  
  c(
    SD = mean(res),
    CI_low = quantile(res, 0.025),
    CI_high = quantile(res, 0.975)
  )
}

###aplicarlo a cada diente

vars <- c("MDCo","BLCo","MDCu","BLCu")

datos_long <- Dent %>%
  filter(Tipo_Diente %in% c("CS","CI")) %>%
  pivot_longer(
    cols = all_of(vars),
    names_to = "Variable",
    values_to = "Valor",
    values_drop_na = TRUE
  )
##boot x diente
garn_results <- datos_long %>%
  group_by(Tipo_Diente, Variable) %>%
  group_modify(~{
    
    boot <- garn_boot(.x)
    
    data.frame(
      SD = boot[1],
      CI_low = boot[2],
      CI_high = boot[3]
    )
    
  }) %>%
  ungroup()

# --- BLOQUE R 012 ----------------------------------------------------------
##Grafico CS-CI con Garn index

p <- ggplot(datos_long, aes(x = Sexo, y = Valor, fill = Sexo)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
  facet_grid(Tipo_Diente ~ Variable, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Sexo",
    y = "Tamaño (mm)"
  ) +
  theme(legend.position = "none")


p +
geom_text(
  data = garn,
  aes(
    x = 1,
    y = Inf,
    label = paste0("SD = ", round(Garn,1), "%")
  ),
  vjust = 1.5,
  inherit.aes = FALSE,
  size = 3
)

##integrar al grafico 
garn_results <- garn_results %>%
  mutate(
    label = paste0(
      "SDI = ", round(SD,1), "%\n",
      "95% CI [", round(CI_low,1), ", ", round(CI_high,1), "]"
    )
  )

p +
geom_text(
  data = garn_results,
  aes(
    x = 1,
    y = Inf,
    label = label
  ),
  vjust = 1.3,
  inherit.aes = FALSE,
  size = 2
)

# --- BLOQUE R 013 ----------------------------------------------------------
ggsave("medidas_CS_CI+Garn.jpg", width = 20, height = 15, units = "cm", dpi = 300)

# --- BLOQUE R 014 ----------------------------------------------------------
#Confidence intervals+media
garn_results <- garn_results %>%
  mutate(Variable = factor(
    Variable,
    levels = c("MDCo","BLCo","MDCu","BLCu")
  ))


p_dim <- ggplot(garn_results,
       aes(x = SD, y = Variable, color = Tipo_Diente)) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey60",
             linewidth = 0.5) +
  
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    width = 0.15,
    linewidth = 1,
    position = position_dodge(width = 0.5)
  ) +
  
  geom_point(
    size = 3,
    position = position_dodge(width = 0.5)
  ) +
  
  theme_classic(base_size = 13) +
  
  labs(
    x = "Dimorfismo sexual (%)",
    y = "Medida dental",
    color = "Diente") +
  
  scale_color_manual(
    values = c(
      "CS" = "#0072B2",
      "CI" = "#D55E00"
    )
  )

p_dim 

# --- BLOQUE R 015 ----------------------------------------------------------
ggsave("Intervalo_confianza_Garn.jpg", width = 20, height = 15, units = "cm", dpi = 300)

# --- BLOQUE R 016 ----------------------------------------------------------
##ANALISIS##
##funcion de division por tipo de diente 80 20 tomando en cuenta por sexo e 
#individuo y controlando por completitud de datos


split_dentales <- function(data, prop = 0.8) {
  
  required_cols <- c("Ind", "Sexo")
  missing <- setdiff(required_cols, names(data))
  if(length(missing) > 0){
    stop(paste("Faltan columnas:", 
               paste(missing, collapse = ", ")))
  }
  
  set.seed(123)
  
  split_obj <- group_initial_split(
    data,
    prop = prop,
    group = Ind,     # x individuo
    strata = Sexo)   # por sexo
  
  list(
    train = training(split_obj),
    test  = testing(split_obj),
    split = split_obj
  )
}



# --- BLOQUE R 017 ----------------------------------------------------------
#division de la muestra con la funcion anterior

cs_split =split_dentales(CS)

cs_train = cs_split$train#entrenamiento
cs_test = cs_split$test#testeo


# --- BLOQUE R 018 ----------------------------------------------------------
summary(cs_train)
summary(cs_test)

# --- BLOQUE R 019 ----------------------------------------------------------
###chequeo de que no estemos repitiendo individuos en el muestreo
length(intersect(unique(cs_train$Ind), unique(cs_test$Ind)))


# --- BLOQUE R 020 ----------------------------------------------------------

###chequeo de los datos a imputar
cs_train %>%
  filter(is.na(MDCo) & is.na(BLCo) & is.na(MDCu) & is.na(BLCu))


# --- BLOQUE R 021 ----------------------------------------------------------
##remover NA
cs_train <- cs_train %>%
  filter(!(is.na(MDCo) & is.na(BLCo) & is.na(MDCu) & is.na(BLCu)))

##chequeo
cs_train %>%
  filter(is.na(MDCo) & is.na(BLCo) & is.na(MDCu) & is.na(BLCu))



# --- BLOQUE R 022 ----------------------------------------------------------
##imputacion y estandarizacion. Aca usamos el knn como se utilizo antes apara 
#armar la formula solo en el training data, para evitar fuga de informacion 
#del test set. Despues vamos a aplicar la misma receta al test set


cs_recipe <- recipe(Sexo ~ MDCo + BLCo + MDCu + BLCu, data = cs_train) %>%
  step_impute_knn(all_predictors()) %>%
  step_normalize(all_predictors())


# --- BLOQUE R 023 ----------------------------------------------------------
###chequeo de los datos
prep(cs_recipe)


# --- BLOQUE R 024 ----------------------------------------------------------
##base de control con validacion cruzada por individuo y sexo para evitar 
#pseudorreplicacion 

set.seed(123)

cv_cs =
  group_vfold_cv(
    cs_train,
    v = 10,
    repeats=5,
    group = Ind,
    strata = Sexo)


# --- BLOQUE R 025 ----------------------------------------------------------
##metricas para extraer

metricas =
  metric_set(
    accuracy,
    sens,
    spec,
    bal_accuracy,
    roc_auc
  )

ctrl =
  control_resamples(save_pred = TRUE)


# --- BLOQUE R 026 ----------------------------------------------------------
##arma el workflow para el modelo
cs_wf =
  workflow() %>%
  add_recipe(cs_recipe)


# --- BLOQUE R 027 ----------------------------------------------------------

####Modelos###############

# Regresion logistica


log_spec =
  logistic_reg(mode = "classification") %>%
  set_engine("glm")

log_wf =
  workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(log_spec)

log_res =
  fit_resamples(
    log_wf,
    resamples = cv_cs,
    control = ctrl,
    metrics = metricas
  )

collect_metrics(log_res)   ##metricas


# --- BLOQUE R 028 ----------------------------------------------------------
# ELASTIC NET

##regresion logistica con regularizacion
enet_spec =
  logistic_reg(
    penalty = tune(),#valor de lamda que va a penalizar los coeficientes para
    #obtener los mas pequeños posibles, en este caso esta en automatico
    mixture = tune()#tipo de penalizacion o alfa, 0 ridge, 1 lasso, entre 0 y 1.
    #automatica tambien
  ) %>%
  set_engine("glmnet") %>% #motor de ajuste del modelo
  set_mode("classification")

enet_wf =
  workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(enet_spec)

#recipe = cómo se procesan datos (imputación + estandarización, etc.)
#model = qué modelo entrenar

enet_grid =
  grid_regular(penalty(range = c(-4, 0)), # grilla de combinaciones de las variable 
  # testeo y  penalizacion (en escala log)
    mixture(range = c(0, 1)), # 0=ridge, 1=lasso. valores intermedios = elastic net
    levels = 10 #se evaluan 10x10 combinaciones alpha y lambda (100 modelos)
  )

set.seed(123)
#todos los modelos se prueban con validacion cruzada y se estima balanced accuracy
#como metrica principal
enet_res =
  tune_grid(
    enet_wf,
    resamples = cv_cs,
    grid = enet_grid,
    control = ctrl,
    metrics = metricas
  )

collect_metrics(enet_res)

best_enet =
  select_best(enet_res, metric = "bal_accuracy") #seleccion y almacenamiento del mejor
  #modelo

best_enet


# --- BLOQUE R 029 ----------------------------------------------------------

# LDA

lda_spec =
  discrim_linear() %>%
  set_engine("MASS")

lda_wf =
  workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(lda_spec)

##prueba de los modelos con validacion cruzada
lda_res =
  fit_resamples(
    lda_wf,
    resamples = cv_cs,
    metrics = metricas,
    control = ctrl
  )

collect_metrics(lda_res)


# --- BLOQUE R 030 ----------------------------------------------------------

# NAIVE BAYES

#especificacion del modelo
nb_spec =
  naive_Bayes() %>%
  set_engine("klaR")##usa paquete klaR para ajuste

#receta del modelo
nb_wf =
  workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(nb_spec)

##prueba de los modelos con validacion cruzada
nb_res =
  fit_resamples(
    nb_wf,
    resamples = cv_cs,
    metrics = metricas,
    control = ctrl)

collect_metrics(nb_res)


# --- BLOQUE R 031 ----------------------------------------------------------

# RANDOM FOREST (con tuning )

rf_spec =
  rand_forest(
    trees = 1000,#1000 arboles estabilizar busqueda y parametros
    mtry = tune(),
    min_n = tune()
  ) %>%
  set_engine("ranger") %>%
  set_mode("classification")

rf_wf =
  workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(rf_spec)

rf_grid =
  grid_regular(
    mtry(range = c(1, 4)),##entrenamos con todos los split para el nodo con 1,
    #2, 3 y hasta 4 predictores, va a dar arboles mas y menos aletorios
    min_n(range = c(2, 20)),
    levels = 6 #da unos 36 modelos o combinaciones x prueba
  )


set.seed(123)

rf_res =
  tune_grid(
    rf_wf,
    resamples = cv_cs,
    grid = rf_grid,
    metrics = metricas,
    control = ctrl)

collect_metrics(rf_res)

best_rf =
  select_best(rf_res, metric = "bal_accuracy")

best_rf


# --- BLOQUE R 032 ----------------------------------------------------------

# SVM lineal


svm_lin =
  svm_linear(cost = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

svm_lin_wf =
  workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(svm_lin)

svm_lin_grid =
  grid_regular(cost(range = c(-3, 2)), levels = 10)

set.seed(123)

svm_lin_res =
  tune_grid(
    svm_lin_wf,
    resamples = cv_cs,
    grid = svm_lin_grid,
    metrics = metricas,
    control = ctrl
  )

collect_metrics(svm_lin_res)

best_svm_lin =
  select_best(svm_lin_res, metric = "bal_accuracy")
best_svm_lin


# --- BLOQUE R 033 ----------------------------------------------------------

# SVM RADIAL CS

svm_rad =
  svm_rbf(
    cost = tune(),
    rbf_sigma = tune()
  ) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

svm_rad_wf =
  workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(svm_rad)

svm_rad_grid =
  grid_regular(
    cost(range = c(-2, 2)),        # 0.01 a 100
    rbf_sigma(range = c(-3, -1)),  # rango más razonable
    levels = 6
  )

set.seed(123)

svm_rad_res =
  tune_grid(
    svm_rad_wf,
    resamples = cv_cs,
    grid = svm_rad_grid,
    metrics = metricas,
    control = ctrl
  )

collect_metrics(svm_rad_res)

best_svm_rad =
  select_best(svm_rad_res, metric = "bal_accuracy")
best_svm_rad


# --- BLOQUE R 034 ----------------------------------------------------------

# XGBOOST CS


xgb_spec =
  boost_tree(
    trees = 800,#800 arboles
    learn_rate = tune(),
    mtry = tune(),
    tree_depth = tune(),
    min_n = tune(),
    loss_reduction = tune()
  ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

xgb_wf =
  workflow() %>%
  add_recipe(cs_recipe) %>%
  add_model(xgb_spec)

xgb_grid =
  grid_space_filling(
    learn_rate(range = c(-3, -1)),      # 0.001 a 0.1
    mtry(range = c(1, 4)),
    tree_depth(range = c(2L, 6L)),      # profundidad moderada
    min_n(range = c(2L, 15L)),
    loss_reduction(range = c(-5, -1)),  # 1e-5 a 0.1
    size = 25
  )

set.seed(123)

xgb_res =
  tune_grid(
    xgb_wf,
    resamples = cv_cs,
    grid = xgb_grid,
    metrics = metricas,
    control = ctrl
  )

collect_metrics(xgb_res)

best_xgb =
  select_best(xgb_res, metric = "bal_accuracy")

xgb_final_wf <- finalize_workflow(xgb_wf, best_xgb)

xgb_fit <- fit(
  xgb_final_wf,
  data = cs_train
)

best_xgb


# --- BLOQUE R 035 ----------------------------------------------------------
##recoleccion de la metrica de clasificacion

## ---- funcion auxiliar para modelos sin tuning ----
extraer_simple <- function(res_obj, nombre_modelo) {
  collect_metrics(res_obj, summarize = FALSE) %>%
    dplyr::filter(.metric == "bal_accuracy") %>%
    dplyr::summarise(
      mean_BA = mean(.estimate),
      sd_BA   = sd(.estimate)
    ) %>%
    dplyr::mutate(modelo = nombre_modelo)
}

## ---- funcion auxiliar para modelos con tuning ----
extraer_tuneado <- function(res_obj, nombre_modelo) {

  # identificar mejor configuración
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

## ---- recoleccion robusta ----
comparacion_modelos <-
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

comparacion_modelos



# --- BLOQUE R 036 ----------------------------------------------------------
##Compracion estadistica entre modelos
#extraccion de parametros estimados

##1- funcion

extraer_resamples_simple <- function(res_obj, nombre_modelo) {
  
  collect_metrics(res_obj, summarize = FALSE) %>%
    filter(.metric == "bal_accuracy") %>%
    select(id, id2, .estimate) %>%
    mutate(modelo = nombre_modelo)
}


extraer_resamples_tuneado <- function(res_obj, nombre_modelo) {
  
  best_id <- select_best(res_obj, metric = "bal_accuracy")$.config
  
  collect_metrics(res_obj, summarize = FALSE) %>%
    filter(.metric == "bal_accuracy",
           .config == best_id) %>%
    select(id, id2, .estimate) %>%
    mutate(modelo = nombre_modelo)
}


##2-objeto
resamples_CS <-
  bind_rows(
    
    extraer_resamples_simple(log_res, "Logistic"),
    extraer_resamples_tuneado(enet_res, "Elastic Net"),
    extraer_resamples_simple(lda_res, "LDA"),
    extraer_resamples_simple(nb_res, "Naive Bayes"),
    extraer_resamples_tuneado(rf_res, "Random Forest"),
    extraer_resamples_tuneado(svm_lin_res, "SVM Linear"),
    extraer_resamples_tuneado(svm_rad_res, "SVM Radial"),
    extraer_resamples_tuneado(xgb_res, "XGBoost")
    
  )


# --- BLOQUE R 037 ----------------------------------------------------------
#check
resamples_CS

# --- BLOQUE R 038 ----------------------------------------------------------
##matriz
resamples_CS_wide <-
  resamples_CS %>%
  select(id, id2, modelo, .estimate) %>%
  tidyr::pivot_wider(
    names_from = modelo,
    values_from = .estimate
  )

# --- BLOQUE R 039 ----------------------------------------------------------
##check (dimensiones y composicion de la matriz)
dim(resamples_CS_wide)
head(resamples_CS_wide)

# --- BLOQUE R 040 ----------------------------------------------------------
##test pareadeo de Friedman (global) con comparaciones de a pares ajustadas x Holm-Bonferroni

friedman.test(
  y = as.matrix(
    resamples_CS_wide %>%
      select(
        Logistic,
        `Elastic Net`,
        LDA,
        `Naive Bayes`,
        `Random Forest`,
        `SVM Linear`,
        `SVM Radial`,
        XGBoost
      )))


# --- BLOQUE R 041 ----------------------------------------------------------
datos_CS <- resamples_CS_wide %>%
  select(
    id, id2,
    Logistic,
    `Elastic Net`,
    LDA,
    `Naive Bayes`,
    `Random Forest`,
    `SVM Linear`,
    `SVM Radial`,
    XGBoost
  )

modelos <- names(datos_CS)[3:10]

comparaciones_CS <- combn(modelos, 2, simplify = FALSE) %>%
  map_dfr(function(par) {
    
    test <- wilcox.test(
      datos_CS[[par[1]]],
      datos_CS[[par[2]]],
      paired = TRUE,
      exact = FALSE
    )
    
    tibble(
      modelo_1 = par[1],
      modelo_2 = par[2],
      p = test$p.value
    )
    
  }) %>%
  mutate(
    p_ajustado = p.adjust(p, method = "holm")
  ) %>%
  arrange(p_ajustado)

##ver
print(comparaciones_CS, n=28)


# --- BLOQUE R 042 ----------------------------------------------------------
##diferencias medias de BA entre modelos
resamples_CS_wide %>%
  summarise(
    XGB_vs_RF = mean(XGBoost - `Random Forest`),
    XGB_vs_SVMrad = mean(XGBoost - `SVM Radial`),
    XGB_vs_NB = mean(XGBoost - `Naive Bayes`),
    XGB_vs_EN = mean(XGBoost - `Elastic Net`),
    XGB_vs_Logistic = mean(XGBoost - Logistic),
    XGB_vs_SVMlin = mean(XGBoost - `SVM Linear`),
    XGB_vs_LDA = mean(XGBoost - LDA)
  )


# --- BLOQUE R 043 ----------------------------------------------------------
##XGBOOSt el mejor modelo, asi que predecimos con el


#  Generar predicciones (clase + probabilidades)

xgb_pred =
  predict(xgb_fit, cs_test, type = "prob") %>%
  bind_cols(predict(xgb_fit, cs_test)) %>%
  bind_cols(dplyr::select(cs_test, Sexo))

# asegurar que Sexo sea factor x las dudas
xgb_pred$Sexo = as.factor(xgb_pred$Sexo)



# Métricas principales de clasificación (Balanced Accuracy = métrica primaria)


metricas_clase =
  metric_set(
    accuracy,
    sens,
    spec,
    bal_accuracy
  )

metricas_clase(
  xgb_pred,
  truth = Sexo,
  estimate = .pred_class
)


# ROC – AUC (probabilidades) Seteo F como la clase positiva

roc_auc(
  xgb_pred,
  truth = Sexo,
  .pred_F)


# --- BLOQUE R 044 ----------------------------------------------------------
#curva roc para explorar los umbrales
roc_curve(
  xgb_pred,
  truth = Sexo,
  .pred_F
) %>%
  autoplot()

# --- BLOQUE R 045 ----------------------------------------------------------
ggsave("ROC_CS.jpg", width = 20, height = 15, units = "cm", dpi = 300)

# --- BLOQUE R 046 ----------------------------------------------------------
#Check sobre los resultados antes del remuestreo
table(xgb_pred$Sexo)

# --- BLOQUE R 047 ----------------------------------------------------------
############################################################
## BOOTSTRAP DEL TEST SET – Intervalos de confianza
############################################################

library(boot)

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

set.seed(123)

boot_res <- boot(
  data = xgb_pred,
  statistic = boot_metricas,
  R = 2000
)

############################################################
## Intervalos percentiles 95%
############################################################

# Balanced Accuracy
boot.ci(boot_res, type = "perc", index = 1)

# Sensibilidad
boot.ci(boot_res, type = "perc", index = 2)

# Especificidad
boot.ci(boot_res, type = "perc", index = 3)

# Accuracy
boot.ci(boot_res, type = "perc", index = 4)


# --- BLOQUE R 048 ----------------------------------------------------------

# Tabla Sensibilidad – Especificidad para grafico

threshold_df =
  roc_curve(xgb_pred, truth = Sexo, .pred_F) %>%
  dplyr::select(.threshold, sensitivity, specificity)

# Gráfico Sensibilidad y Especificidad

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



# --- BLOQUE R 049 ----------------------------------------------------------

#  Umbral más Balanceado (Sens y Spec)

best_threshold =
  threshold_df %>%
  dplyr::mutate(diff = abs(sensitivity - specificity)) %>%
  dplyr::arrange(diff) %>%
  dplyr::slice(1)

best_threshold

# --- BLOQUE R 050 ----------------------------------------------------------

#Umbral de Youden J busca el mas explicativo, no necesariamente balanceado

best_youden =
  threshold_df %>%
  mutate(J = sensitivity + specificity - 1) %>%
  arrange(desc(J)) %>%
 dplyr::slice(1)

best_youden


# --- BLOQUE R 051 ----------------------------------------------------------
##matriz de confusion para predicciones XGBOOST
conf_mat(
  xgb_pred,
  truth = Sexo,
  estimate = .pred_class
)


# --- BLOQUE R 052 ----------------------------------------------------------
##importancia de la variables con VIP
xgb_importance =
  xgb_fit %>%
  extract_fit_parsnip() %>%
  vip::vi()

xgb_importance

# --- BLOQUE R 053 ----------------------------------------------------------
#grafico
vip(
  xgb_fit,
  num_features = 10,
  geom = "col",
  aesthetics = list(fill = "steelblue")
) +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Importancia")


# --- BLOQUE R 054 ----------------------------------------------------------
ggsave("Var_Imp_CS.jpg", width = 20, height = 15, units = "cm", dpi = 300)


# --- BLOQUE R 055 ----------------------------------------------------------

###valores del modelo y contribucion de los casos con SHAP

# Extraer el modelo XGBoost

xgb_core = extract_fit_parsnip(xgb_fit)$fit

# Obtener los datos procesados por la formula (recipe)

rec = extract_recipe(xgb_fit)

X_df = bake(rec, new_data = cs_train) %>%
  dplyr::select(-Sexo)

# Convertir a matriz numérica (XGBoost requiere tipos numéricos) sino da error
X_matrix = as.matrix(X_df)

## Calcular SHAP 

sv = shapviz(
  xgb_core,
  X = X_matrix,      # Datos para los gráficos
  X_pred = X_matrix, # Datos para la predicción (debe ser matriz)
  baseline = "auto"
)


# --- BLOQUE R 056 ----------------------------------------------------------
# Gráfico de importancia de variables SHAP
sv_importance(sv, kind = "beeswarm")

# --- BLOQUE R 057 ----------------------------------------------------------
ggsave("Shap_CS.jpg", width = 20, height = 15, units = "cm", dpi = 300)

# --- BLOQUE R 058 ----------------------------------------------------------
# Gráficos de dependencia SHAP para cada variable

vars = c("MDCu", "BLCu", "BLCo", "MDCo")

plots =
  lapply(vars, function(v){
    sv_dependence(sv, v) +
      ggtitle(paste("Dependencia SHAP -", v)) +
      theme_minimal(base_size = 12)
  })

wrap_plots(plots, ncol = 2)   # equivalente a mfrow=c(2,2)

# --- BLOQUE R 059 ----------------------------------------------------------
################Ajuste para CI#################
CI <- Dent %>%
  filter(Tipo_Diente == "CI")

#chequeando
#cuántos dientes,
#distribución por sexo,
#n individuos

summary(CI[,c("MDCo","BLCo","MDCu","BLCu")])
table(CI$Sexo)
length(unique(CI$Ind))



# --- BLOQUE R 060 ----------------------------------------------------------
##archivo de CI para análisis posteriores

write.table(
  CI,
  file = "caninos_inferiores.txt",
  sep = "\t",
  row.names = FALSE)

# --- BLOQUE R 061 ----------------------------------------------------------
##division muestra con la funcion split_dentales

ci_split =split_dentales(CI)

ci_train = ci_split$train
ci_test = ci_split$test


# --- BLOQUE R 062 ----------------------------------------------------------
summary(ci_train) 
summary(ci_test)

# --- BLOQUE R 063 ----------------------------------------------------------
###chequeo de que no estemos repitiendo individuos en el muestreo
length(intersect(unique(ci_train$Ind), unique(ci_test$Ind)))


# --- BLOQUE R 064 ----------------------------------------------------------
##imputacion y estandarizacion para CI (train) . 

ci_recipe <- recipe(Sexo ~ MDCo + BLCo + MDCu + BLCu, data = ci_train) %>%
  step_impute_knn(all_predictors()) %>%
  step_normalize(all_predictors())


# --- BLOQUE R 065 ----------------------------------------------------------
###chequeo de los datos
prep(ci_recipe)


# --- BLOQUE R 066 ----------------------------------------------------------
##base de control con validacion cruzada por individuo y sexo para evitar pseudorreplicacion para CI igual que CS

set.seed(123)

cv_ci =
  group_vfold_cv(
    ci_train,
    v = 10,
    repeats=5,
    group = Ind,
    strata = Sexo)


# --- BLOQUE R 067 ----------------------------------------------------------
##metricas para extraer de CI la formula es la misma

metricas =
  metric_set(
    accuracy,
    sens,
    spec,
    bal_accuracy,
    roc_auc
  )

ctrl =
  control_resamples(save_pred = TRUE)


# --- BLOQUE R 068 ----------------------------------------------------------
##arma el workflow para el modelo de CI
ci_wf =
  workflow() %>%
  add_recipe(ci_recipe)


# --- BLOQUE R 069 ----------------------------------------------------------

####Modelos para CI###############

# Regresion logistica


log_spec_CI =
  logistic_reg(mode = "classification") %>%
  set_engine("glm")

log_wf_CI =
  workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(log_spec_CI)

log_res_CI =
  fit_resamples(
    log_wf_CI,
    resamples = cv_ci,
    control = ctrl,
    metrics = metricas
  )

collect_metrics(log_res_CI)   ##metricas CI


# --- BLOQUE R 070 ----------------------------------------------------------

# LDA para CI

lda_spec =#igual
  discrim_linear() %>%
  set_engine("MASS")

lda_wf_CI =
  workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(lda_spec)

##prueba de los modelos con validacion cruzada igual
lda_res_CI =
  fit_resamples(
    lda_wf_CI,
    resamples = cv_ci,
    metrics = metricas,
    control = ctrl
  )

collect_metrics(lda_res_CI)


# --- BLOQUE R 071 ----------------------------------------------------------

# NAIVE BAYES para CI

#especificacion del modelo
nb_spec =#mismo motor
  naive_Bayes() %>%
  set_engine("klaR")##usa paquete klaR para ajuste

#receta del modelo
nb_wf_CI =
  workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(nb_spec)

##prueba de los modelos con validacion cruzada
nb_res_CI =
  fit_resamples(
    nb_wf_CI,
    resamples = cv_ci,
    metrics = metricas,
    control = ctrl)

nb_fit_CI =
  fit(
    nb_wf_CI,
    data = ci_train
  )


collect_metrics(nb_res_CI)


# --- BLOQUE R 072 ----------------------------------------------------------

# RANDOM FOREST (con tuning ) para CI

rf_spec =#idem motor de busqueda
  rand_forest(
    trees = 1000,
    mtry = tune(),
    min_n = tune()
  ) %>%
  set_engine("ranger") %>%
  set_mode("classification")

rf_wf_CI =
  workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(rf_spec)

rf_grid_CI =
  grid_regular(
    mtry(range = c(1, 4)),##entrenamos con todos los split para el nodo con 1,
    #2, 3 y hasta 4 predictores, va a dar arboles mas y menos aletorios
    min_n(range = c(2, 20)),
    levels = 6 #da unos 36 modelos o combinaciones x prueba
  )


set.seed(123)

rf_res_CI =
  tune_grid(
    rf_wf_CI,
    resamples = cv_ci,
    grid = rf_grid_CI,
    metrics = metricas,
    control = ctrl)

collect_metrics(rf_res_CI)

best_rf_CI =
  select_best(rf_res_CI, metric = "bal_accuracy")

best_rf_CI


# --- BLOQUE R 073 ----------------------------------------------------------
# ELASTIC NET para CI

##regresion logistica con regularizacion, esto todo igual porque aplico el mismo modelo
enet_spec =
  logistic_reg(
    penalty = tune(),#valor de lamda que va a penalizar los coeficientes para obtener los mas pequeños posibles, en este caso esta en automatico
    mixture = tune()#tipo de penalizacion o alfa, 0 ridge, 1 lasso, entre 0 y 1. automatica tambien
  ) %>%
  set_engine("glmnet") %>% #motor de ajuste del modelo
  set_mode("classification")

enet_wf_CI =
  workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(enet_spec)

#recipe = cómo se procesan datos (imputación + estandarización, etc.)
#model = qué modelo entrenar

enet_grid =#la grilla es igual ya que usamos siempre los mismos hiperparametros
  grid_regular(penalty(range = c(-4, 0)), # grilla de combinaciones de las variable a testea y  penalizacion (en escala log)
    mixture(range = c(0, 1)), # 0=ridge, 1=lasso. valores intermedios = elastic net
    levels = 10 #se evaluan 10x10 combinaciones alpha y lambda (100 modelos)
  )

set.seed(123)
##todos los modelos se prueban con validacion cruzada y se estima balanced accuracy como metrica principal
enet_res_CI =
  tune_grid(
    enet_wf_CI,
    resamples = cv_ci,
    grid = enet_grid,
    control = ctrl,
    metrics = metricas
  )

collect_metrics(enet_res_CI)

best_enet_CI =
  select_best(enet_res_CI, metric = "bal_accuracy") #seleccion y almacenamiento del mejor modelo

best_enet_CI


# --- BLOQUE R 074 ----------------------------------------------------------

# SVM lineal CI


svm_lin = #tuning igual
  svm_linear(cost = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

svm_lin_wf_CI = 
  workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(svm_lin)#tuning igual

svm_lin_grid =#idem de parametros
  grid_regular(cost(range = c(-3, 2)), levels = 10)


set.seed(123)

svm_lin_res_CI =
  tune_grid(
    svm_lin_wf_CI,
    resamples = cv_ci,
    grid = svm_lin_grid,#idem
    metrics = metricas,
    control = ctrl
  )

collect_metrics(svm_lin_res_CI)

best_svm_lin_CI =
  select_best(svm_lin_res_CI, metric = "bal_accuracy")

best_svm_lin_CI


# --- BLOQUE R 075 ----------------------------------------------------------

# SVM RADIAL CI

svm_rad =#motor y modo igual
  svm_rbf(
    cost = tune(),
    rbf_sigma = tune()
  ) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

svm_rad_wf_CI =
  workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(svm_rad)

svm_rad_grid = #iguales hiperparametros
  grid_regular(
    cost(range = c(-2, 2)),       
    rbf_sigma(range = c(-3, -1)),  
    levels = 6
  )

set.seed(123)

svm_rad_res_CI =
  tune_grid(
    svm_rad_wf_CI,
    resamples = cv_ci,
    grid = svm_rad_grid,
    metrics = metricas,
    control = ctrl
  )

collect_metrics(svm_rad_res_CI)

best_svm_rad_CI =
  select_best(svm_rad_res_CI, metric = "bal_accuracy")

best_svm_rad_CI


# --- BLOQUE R 076 ----------------------------------------------------------

# XGBOOST CI


xgb_spec =#modelo igual
  boost_tree(
    trees = 800,
    learn_rate = tune(),
    mtry = tune(),
    tree_depth = tune(),
    min_n = tune(),
    loss_reduction = tune()
  ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

xgb_wf_CI =
  workflow() %>%
  add_recipe(ci_recipe) %>%
  add_model(xgb_spec)

xgb_grid_CI =
  grid_space_filling(
    learn_rate(range = c(-3, -1)),      # 0.001 a 0.1
    mtry(range = c(1, 4)),
    tree_depth(range = c(2L, 6L)),      # profundidad moderada
    min_n(range = c(2L, 15L)),
    loss_reduction(range = c(-5, -1)),  # 1e-5 a 0.1
    size = 25
  )

set.seed(123)

xgb_res_CI =
  tune_grid(
    xgb_wf_CI,
    resamples = cv_ci,
    grid = xgb_grid_CI,
    metrics = metricas,
    control = ctrl
  )

collect_metrics(xgb_res_CI)

best_xgb_CI =
  select_best(xgb_res_CI, metric = "bal_accuracy")

best_xgb_CI


# --- BLOQUE R 077 ----------------------------------------------------------
##recoleccion de la metrica de clasificacion para CI

############################################################
## FUNCIONES AUXILIARES (si ya las definiste en CS, no repetir)
############################################################

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
    dplyr::filter(.metric == "bal_accuracy",
                  .config == best_id) %>%
    dplyr::summarise(
      mean_BA = mean(.estimate),
      sd_BA   = sd(.estimate)
    ) %>%
    dplyr::mutate(modelo = nombre_modelo)
}

############################################################
## RECOLECCIÓN ROBUSTA DE MÉTRICAS – CI
############################################################

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

comparacion_modelos_CI




# --- BLOQUE R 078 ----------------------------------------------------------
## Comparación estadística entre modelos - CI

resamples_CI <-
  bind_rows(
    
    extraer_resamples_simple(log_res_CI, "Logistic"),
    extraer_resamples_tuneado(enet_res_CI, "Elastic Net"),
    extraer_resamples_simple(lda_res_CI, "LDA"),
    extraer_resamples_simple(nb_res_CI, "Naive Bayes"),
    extraer_resamples_tuneado(rf_res_CI, "Random Forest"),
    extraer_resamples_tuneado(svm_lin_res_CI, "SVM Linear"),
    extraer_resamples_tuneado(svm_rad_res_CI, "SVM Radial"),
    extraer_resamples_tuneado(xgb_res_CI, "XGBoost")
    
  )

## check
resamples_CI

# --- BLOQUE R 079 ----------------------------------------------------------
## matriz

resamples_CI_wide <-
  resamples_CI %>%
  select(id, id2, modelo, .estimate) %>%
  tidyr::pivot_wider(
    names_from = modelo,
    values_from = .estimate
  )

## check
dim(resamples_CI_wide)
head(resamples_CI_wide)

# --- BLOQUE R 080 ----------------------------------------------------------
##Friedman (global y comparacion de a pares)
friedman.test(
  y = as.matrix(
    resamples_CI_wide %>%
      select(
        Logistic,
        `Elastic Net`,
        LDA,
        `Naive Bayes`,
        `Random Forest`,
        `SVM Linear`,
        `SVM Radial`,
        XGBoost
      ))
)

# --- BLOQUE R 081 ----------------------------------------------------------

##Wilcoxon pareado + Holm entre pares

datos_CI <- resamples_CI_wide %>%
  select(
    id, id2,
    Logistic,
    `Elastic Net`,
    LDA,
    `Naive Bayes`,
    `Random Forest`,
    `SVM Linear`,
    `SVM Radial`,
    XGBoost
  )

modelos <- names(datos_CI)[3:10]

comparaciones_CI <- combn(modelos, 2, simplify = FALSE) %>%
  map_dfr(function(par) {
    
    test <- wilcox.test(
      datos_CI[[par[1]]],
      datos_CI[[par[2]]],
      paired = TRUE,
      exact = FALSE
    )
    
    tibble(
      modelo_1 = par[1],
      modelo_2 = par[2],
      p = test$p.value
    )
    
  }) %>%
  mutate(
    p_ajustado = p.adjust(p, method = "holm")
  ) %>%
  arrange(p_ajustado)

##ver
print(comparaciones_CI, n=28)


# --- BLOQUE R 082 ----------------------------------------------------------
##diferencias medias de BA entre modelos
resamples_CI_wide %>%
  summarise(
    XGB_vs_RF = mean(XGBoost - `Random Forest`),
    XGB_vs_SVMrad = mean(XGBoost - `SVM Radial`),
    XGB_vs_NB = mean(XGBoost - `Naive Bayes`),
    XGB_vs_EN = mean(XGBoost - `Elastic Net`),
    XGB_vs_Logistic = mean(XGBoost - Logistic),
    XGB_vs_SVMlin = mean(XGBoost - `SVM Linear`),
    XGB_vs_LDA = mean(XGBoost - LDA)
  )

# --- BLOQUE R 083 ----------------------------------------------------------
##Extraer mejor SVM radial del modelo entrenado
best_svm_rad_CI <- select_best(
  svm_rad_res_CI,
  metric = "bal_accuracy"
)

#cerrar workflow
svm_rad_final_wf_CI <-
  finalize_workflow(
    svm_rad_wf_CI,
    best_svm_rad_CI
  )

#volver a correr con training
svm_rad_fit_CI <-
  fit(
    svm_rad_final_wf_CI,
    data = ci_train
  )


# --- BLOQUE R 084 ----------------------------------------------------------
##Predecir con SVM radial sobre test set

svm_pred_CI =
  predict(svm_rad_fit_CI, ci_test, type = "prob") %>%
  bind_cols(predict(svm_rad_fit_CI, ci_test)) %>%
  bind_cols(dplyr::select(ci_test, Sexo))

svm_pred_CI$Sexo = as.factor(svm_pred_CI$Sexo)


# --- BLOQUE R 085 ----------------------------------------------------------
##Metricas para SVM radial

metricas_clase(
  svm_pred_CI,
  truth = Sexo,
  estimate = .pred_class
)

roc_auc(
  svm_pred_CI,
  truth = Sexo,
  .pred_F
)


# --- BLOQUE R 086 ----------------------------------------------------------
############################################################
## BOOTSTRAP DEL TEST SET – CI (SVM Radial)
############################################################

boot_res_CI <- boot(
  data = svm_pred_CI,
  statistic = boot_metricas,
  R = 2000
)

# Balanced Accuracy
boot.ci(boot_res_CI, type = "perc", index = 1)

# Sensibilidad
boot.ci(boot_res_CI, type = "perc", index = 2)

# Especificidad
boot.ci(boot_res_CI, type = "perc", index = 3)

# Accuracy
boot.ci(boot_res_CI, type = "perc", index = 4)


# --- BLOQUE R 087 ----------------------------------------------------------
roc_curve(
  svm_pred_CI,
  truth = Sexo,
  .pred_F
) %>%
  autoplot()


# --- BLOQUE R 088 ----------------------------------------------------------
ggsave("ROC_CI.jpg", width = 20, height = 15, units = "cm", dpi = 300)

# --- BLOQUE R 089 ----------------------------------------------------------

# Tabla Sensibilidad – Especificidad

threshold_df_CI =
  roc_curve(svm_pred_CI, truth = Sexo, .pred_F) %>%
  dplyr::select(.threshold, sensitivity, specificity)


# Gráfico

ggplot(threshold_df_CI, aes(x = .threshold)) +
  geom_line(aes(y = sensitivity, color = "Sensibilidad"), size = 1) +
  geom_line(aes(y = specificity, color = "Especificidad"), size = 1) +
  scale_color_manual(values = c("Sensibilidad" = "red",
                                "Especificidad" = "blue")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Umbral",
    y = "Valor",
    color = "",
    title = "Sensibilidad y Especificidad según Umbral (SVM Radial)"
  )


# --- BLOQUE R 090 ----------------------------------------------------------
##DEistribucion de probabilidades por clase
ggplot(svm_pred_CI,
       aes(x = .pred_F,
           fill = Sexo)) +
  geom_density(alpha = 0.4) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Probabilidad predicha de F",
    y = "Densidad"
  )


# --- BLOQUE R 091 ----------------------------------------------------------
ggsave("Dist_prob.jpg", width = 20, height = 15, units = "cm", dpi = 300)

# --- BLOQUE R 092 ----------------------------------------------------------
##usar un esquema de permutaciones para estimar la importancia de las variables SVM
##library(vip)#Uso vip pero hay conflicto entre paquetes recordar vip::vi


ci_test_imp <- ci_test %>%
  dplyr::select(Sexo, MDCo, BLCo, MDCu, BLCu)

set.seed(123)
imp_ci <- vip::vi(
  svm_rad_fit_CI,
  method = "permute",
  train = ci_test_imp ,
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

#check
print(imp_ci_summary) 

# --- BLOQUE R 093 ----------------------------------------------------------
ggplot(imp_ci_summary,
       aes(x = reorder(Variable, mean_imp),
           y = mean_imp)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_imp - sd_imp,
                    ymax = mean_imp + sd_imp),
                width = 0.2) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = "",
    y = "Importancia (Δ Balanced Accuracy)"
  )


# --- BLOQUE R 094 ----------------------------------------------------------
ggsave("Variable_imp_CI.jpg", width = 20, height = 15, units = "cm", dpi = 300)

# --- BLOQUE R 095 ----------------------------------------------------------
#Grafico del espacio clasificatorio para MDCu y MDCo

grid_data <- expand.grid(
  MDCu = seq(min(ci_test$MDCu, na.rm = TRUE),
             max(ci_test$MDCu, na.rm = TRUE),
             length.out = 100),
  MDCo = seq(min(ci_test$MDCo, na.rm = TRUE),
             max(ci_test$MDCo, na.rm = TRUE),
             length.out = 100)
)

# fijamos las otras variables en la media
grid_data$BLCu <- mean(ci_test$BLCu, na.rm = TRUE)
grid_data$BLCo <- mean(ci_test$BLCo, na.rm = TRUE)


# --- BLOQUE R 096 ----------------------------------------------------------
#predecir grilla
grid_pred <- predict(
  svm_rad_fit_CI,
  grid_data,
  type = "prob"
)


grid_plot <- cbind(grid_data, grid_pred)


# --- BLOQUE R 097 ----------------------------------------------------------
#Filtrar valores NA

ci_test_clean <- ci_test %>%
  filter(!is.na(MDCu), !is.na(MDCo))

#plot
ggplot() +
  geom_raster(
    data = grid_plot,
    aes(x = MDCu,
        y = MDCo,
        fill = .pred_F)
  ) +
  scale_fill_viridis_c() +
  geom_point(
    data = ci_test_clean,
    aes(x = MDCu,
        y = MDCo,
        shape = Sexo),
    size = 3,
    color = "black"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "MDCu",
    y = "MDCo",
    fill = "P(F)"
  )



# --- BLOQUE R 098 ----------------------------------------------------------
ggsave("SVM_sup.jpg", width = 20, height = 15, units = "cm", dpi = 300)

# --- BLOQUE R 099 ----------------------------------------------------------
##salvar modelo CS
saveRDS(xgb_fit, "xgb_fit_CS.rds")
saveRDS(svm_rad_fit_CI, "svm_rad_fit_CI.rds")

#test data cs
write.table(cs_test,"Conjunto_testeo.txt", sep = "\t",row.names = FALSE,
            quote = FALSE)

#salvar modelo CI
write.table(ci_test,"Conjunto_testeo_CI.txt", sep = "\t",row.names = FALSE,
            quote = FALSE)


# --- BLOQUE R 100 ----------------------------------------------------------
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
    file = "session_info_seleccion de modelos.txt", append = FALSE)


# Tiempo de ejecución aproximado: ~1800-2000 segundos (30-35 minutos)
