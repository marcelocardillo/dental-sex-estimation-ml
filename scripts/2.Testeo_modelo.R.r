# ==============================================================================
# SCRIPT INTEGRADO: TESTEO DE MODELO PREDICTIVO Y ROBUSTEZ (VER 2)
# Descripción: Evaluación del rol de datos ausentes en modelos predictivos 
#              (XGBoost para CS y SVM Radial para CI) mediante escenarios 
#              contrafactuales y bootstrap.
# ==============================================================================

# --- 1. CARGA DE PAQUETES -----------------------------------------------------
# Se cargan todas las librerías necesarias para el flujo de tidymodels, 
# manipulación de datos y visualización.
library(tidymodels)
library(tidyverse)
library(recipes)
library(parsnip)
library(yardstick)
library(purrr)
library(tibble)
library(ggplot2)
library(dplyr)
library(boot)

# --- 2. CARGA DE MODELOS Y DATOS ----------------------------------------------
# Se cargan los objetos modelo entrenados previamente y los conjuntos de testeo.
# NOTA: Asegúrese de que estos archivos existan en el directorio de trabajo.
xgb_fit_CS      <- readRDS("xgb_fit_CS.rds")
svm_rad_fit_CI  <- readRDS("svm_rad_fit_CI.rds")

# Carga de conjuntos de testeo originales
cs_test <- read.table("Conjunto_testeo_CS.txt", header = TRUE)
ci_test <- read.table("Conjunto_testeo_CI.txt", header = TRUE)

# Asegurar que la variable objetivo 'Sexo' sea factor
cs_test$Sexo <- as.factor(cs_test$Sexo)
ci_test$Sexo <- as.factor(ci_test$Sexo)

# --- 3. EVALUACIÓN INICIAL (BASELINE) -----------------------------------------
# Verificación rápida del desempeño original antes de introducir datos faltantes.

# Modelo CS (XGBoost)
pred_cs <- predict(xgb_fit_CS, cs_test) %>%
  bind_cols(cs_test["Sexo"])
ba_cs_init <- bal_accuracy(pred_cs, truth = Sexo, estimate = .pred_class)

# Modelo CI (SVM Radial)
pred_ci <- predict(svm_rad_fit_CI, ci_test) %>%
  bind_cols(ci_test["Sexo"])
ba_ci_init <- bal_accuracy(pred_ci, truth = Sexo, estimate = .pred_class)

# --- 4. PREPARACIÓN DE DATOS (LIMPIEZA DE NA) ---------------------------------
# Se filtran los datos para mantener solo observaciones completas en las variables
# clave para la simulación. Esto genera el baseline limpio.
vars_simulacion <- c("MDCo", "BLCo", "MDCu", "BLCu")

cs_test_complete <- cs_test %>%
  filter(if_all(all_of(vars_simulacion), ~ !is.na(.)))

ci_test_complete <- ci_test %>%
  filter(if_all(all_of(vars_simulacion), ~ !is.na(.)))

# --- 5. FUNCIONES AUXILIARES --------------------------------------------------

# Función para generar combinaciones de variables faltantes (escenarios)
make_na_scenarios <- function(vars) {
  combos_1 <- combn(vars, 1, simplify = FALSE)
  combos_2 <- combn(vars, 2, simplify = FALSE)
  combos_3 <- combn(vars, 3, simplify = FALSE)
  # Se podrían incluir combos_4 si se desea probar todas faltantes
  all_combos <- c(combos_1, combos_2, combos_3)
  names(all_combos) <- sapply(all_combos, paste, collapse = "_")
  return(all_combos)
}

# Función principal para evaluar robustez introduciendo NA sistemáticos
evaluar_robustez_full <- function(modelo_fit, test_data, nombre_modelo) {
  vars <- c("MDCo", "BLCo", "MDCu", "BLCu") # Sin espacios trailing
  
  # --- Baseline (Sin NA) ---
  baseline_pred <- predict(modelo_fit, test_data) %>%
    bind_cols(test_data["Sexo"])
  
  baseline_BA   <- yardstick::bal_accuracy(baseline_pred, truth = Sexo, estimate = .pred_class)$.estimate
  baseline_Sens <- yardstick::sens(baseline_pred, truth = Sexo, estimate = .pred_class)$.estimate
  baseline_Spec <- yardstick::spec(baseline_pred, truth = Sexo, estimate = .pred_class)$.estimate
  
  baseline_row <- tibble::tibble(
    modelo    = nombre_modelo,
    escenario = "None",
    n_na      = 0,
    BA        = baseline_BA,
    Sens      = baseline_Sens,
    Spec      = baseline_Spec
  )
  
  # --- Escenarios con NA ---
  scenarios <- make_na_scenarios(vars)
  
  resultados_na <- purrr::imap_dfr(scenarios, function(miss_vars, scenario_name) {
    test_mod <- test_data
    test_mod[, miss_vars] <- NA # Introducir missing
    
    pred <- predict(modelo_fit, test_mod) %>%
      bind_cols(test_mod["Sexo"])
    
    tibble::tibble(
      modelo    = nombre_modelo,
      escenario = scenario_name,
      n_na      = length(miss_vars),
      BA        = yardstick::bal_accuracy(pred, truth = Sexo, estimate = .pred_class)$.estimate,
      Sens      = yardstick::sens(pred, truth = Sexo, estimate = .pred_class)$.estimate,
      Spec      = yardstick::spec(pred, truth = Sexo, estimate = .pred_class)$.estimate
    )
  })
  
  # --- Unir Resultados ---
  resultados_all <- dplyr::bind_rows(baseline_row, resultados_na) %>%
    mutate(dBA = BA - baseline_BA) # Calcular degradación
  
  return(resultados_all)
}

# --- 6. EJECUCIÓN DE ROBUSTEZ (CS y CI) ---------------------------------------

# Modelo CS
robustez_CS <- evaluar_robustez_full(
  modelo_fit  = xgb_fit_CS,
  test_data   = cs_test_complete,
  nombre_modelo = "CS_XGB"
)
write.table(robustez_CS, "robustez_CS.txt", sep = "\t", quote = FALSE)

# Modelo CI
robustez_CI <- evaluar_robustez_full(
  modelo_fit  = svm_rad_fit_CI,
  test_data   = ci_test_complete,
  nombre_modelo = "CI_SVM"
)
write.table(robustez_CI, "robustez_CI.txt", sep = "\t", quote = FALSE)

# --- 7. VISUALIZACIÓN DE ROBUSTEZ ---------------------------------------------

# Función para preparar tablas de gráfico (excluye baseline)
preparar_tabla_grafico <- function(robustez_df) {
  robustez_df %>%
    filter(n_na > 0) %>%
    mutate(
      BA   = round(BA, 2),
      Sens = round(Sens, 2),
      Spec = round(Spec, 2),
      dBA  = round(dBA, 2)
    ) %>%
    arrange(n_na)
}

robustez_CS_tabla <- preparar_tabla_grafico(robustez_CS)
robustez_CI_tabla <- preparar_tabla_grafico(robustez_CI)

# Gráfico Boxplot CS
ggplot(robustez_CS_tabla, aes(x = factor(n_na), y = dBA)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") + # Corregido geom_hline
  theme_minimal(base_size = 14) +
  labs(x = "Número de variables faltantes", y = "Degradación (Δ Balanced Accuracy)",
       title = "Robustez del modelo CS (XGBoost)")

# Gráfico Boxplot CI
ggplot(robustez_CI_tabla, aes(x = factor(n_na), y = dBA)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(x = "Número de variables faltantes", y = "Degradación (Δ Balanced Accuracy)",
       title = "Robustez del modelo CI (SVM radial)")

# Gráfico Comparativo
robustez_total <- bind_rows(robustez_CS_tabla, robustez_CI_tabla)

ggplot(robustez_total, aes(x = factor(n_na), y = dBA, fill = modelo)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(x = "Número de variables faltantes", y = "Degradación (Δ Balanced Accuracy)",
       title = "Comparación de robustez: CS (XGBoost) vs CI (SVM radial)")

# --- 8. REGLAS DE CONFIANZA ---------------------------------------------------

# Clasificación de confianza basada en Balanced Accuracy
rules_confianza_CS <- robustez_CS_tabla %>%
  mutate(nivel_confianza = case_when(
    BA >= 0.75 ~ "Alta",
    BA >= 0.70 ~ "Moderada",
    TRUE       ~ "Baja"
  )) %>%
  arrange(n_na, BA)
write.table(rules_confianza_CS, "rules_confianza_CS.txt", sep = "\t", quote = FALSE)

rules_confianza_CI <- robustez_CI_tabla %>%
  mutate(nivel_confianza = case_when(
    BA >= 0.75 ~ "Alta",
    BA >= 0.70 ~ "Moderada",
    TRUE       ~ "Baja"
  )) %>%
  arrange(n_na, BA)
write.table(rules_confianza_CI, "rules_confianza_CI.txt", sep = "\t", quote = FALSE)

# --- 9. ANÁLISIS BOOTSTRAP ----------------------------------------------------

# Función para calcular dBA en cada réplica de bootstrap
boot_dBA_general <- function(data, indices, miss_vars, modelo_fit) {
  d <- data[indices, ]
  
  # Baseline en la muestra bootstrap
  pred_base <- predict(modelo_fit, d) %>% bind_cols(d["Sexo"])
  BA_base   <- yardstick::bal_accuracy(pred_base, truth = Sexo, estimate = .pred_class)$.estimate
  
  # Modelo con NA en la muestra bootstrap
  d_mod <- d
  d_mod[, miss_vars] <- NA
  pred_mod <- predict(modelo_fit, d_mod) %>% bind_cols(d_mod["Sexo"])
  BA_mod   <- yardstick::bal_accuracy(pred_mod, truth = Sexo, estimate = .pred_class)$.estimate
  
  return(BA_mod - BA_base)
}

# Configuración de escenarios
vars <- c("MDCo", "BLCo", "MDCu", "BLCu")
scenarios <- make_na_scenarios(vars)
set.seed(123)

# Bootstrap CS
bootstrap_results_CS <- purrr::imap_dfr(scenarios, function(miss_vars, scenario_name) {
  boot_obj <- boot(
    data = cs_test_complete,
    statistic = function(data, indices) {
      boot_dBA_general(data, indices, miss_vars, xgb_fit_CS)
    },
    R = 5000
  )
  ci <- boot.ci(boot_obj, type = "perc")
  
  tibble::tibble(
    escenario   = scenario_name,
    n_na        = length(miss_vars),
    mean_dBA    = mean(boot_obj$t),
    lower_95    = ci$percent[4], # Nota: Depende de la estructura de salida de boot.ci
    upper_95    = ci$percent[5],
    incluye_0   = ci$percent[4] <= 0 & ci$percent[5] >= 0
  )
})

# Bootstrap CI
bootstrap_results_CI <- purrr::imap_dfr(scenarios, function(miss_vars, scenario_name) {
  boot_obj <- boot(
    data = ci_test_complete,
    statistic = function(data, indices) {
      boot_dBA_general(data, indices, miss_vars, svm_rad_fit_CI)
    },
    R = 5000
  )
  ci <- boot.ci(boot_obj, type = "perc")
  
  tibble::tibble(
    escenario   = scenario_name,
    n_na        = length(miss_vars),
    mean_dBA    = mean(boot_obj$t),
    lower_95    = ci$percent[4],
    upper_95    = ci$percent[5],
    incluye_0   = ci$percent[4] <= 0 & ci$percent[5] >= 0
  )
})

# --- 10. VISUALIZACIÓN BOOTSTRAP ----------------------------------------------

# Gráfico CS
ggplot(bootstrap_results_CS, aes(x = reorder(escenario, mean_dBA), y = mean_dBA, color = incluye_0)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2) +
  scale_color_manual(values = c("FALSE" = "firebrick", "TRUE" = "gray40"),
                     labels = c("FALSE" = "Significativa", "TRUE" = "No significativa")) +
  coord_flip() + theme_minimal(base_size = 14) +
  labs(x = "Escenarios", y = expression(Delta * " Balanced Accuracy"),
       color = "Resultado", title = "Bootstrap ΔBA – Robustez CS (XGBoost)")

# Gráfico CI
ggplot(bootstrap_results_CI, aes(x = reorder(escenario, mean_dBA), y = mean_dBA, color = incluye_0)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2) +
  scale_color_manual(values = c("FALSE" = "firebrick", "TRUE" = "gray40"),
                     labels = c("FALSE" = "Significativa", "TRUE" = "No significativa")) +
  coord_flip() + theme_minimal(base_size = 14) +
  labs(x = "Escenarios", y = expression(Delta * " Balanced Accuracy"),
       color = "Resultado", title = "Bootstrap ΔBA – Robustez CI (SVM radial)")

# Gráfico Comparativo Bootstrap
bootstrap_results_CS$modelo <- "CS_XGB"
bootstrap_results_CI$modelo <- "CI_SVM"

# Corrección: Usar bootstrap_results_CI (mayúsculas) en lugar de ci (minúsculas)
bootstrap_total <- bind_rows(bootstrap_results_CS, bootstrap_results_CI)

ggplot(bootstrap_total, aes(x = mean_dBA, y = reorder(escenario, mean_dBA), color = modelo)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95),
                position = position_dodge(width = 0.6), width = 0.2) +
  theme_minimal(base_size = 14) +
  labs(x = "dBA (Δ Balanced Accuracy)", y = "Escenario", color = "Modelo",
       title = "Comparación CS (XGBoost) vs CI (SVM radial)")

# --- 11. TEST ESTADÍSTICO PAREADO ---------------------------------------------

# Preparar datos para test pareado
robustez_pareado <- robustez_total %>%
  select(modelo, escenario, dBA) %>%
  pivot_wider(names_from = modelo, values_from = dBA)

# Test de Wilcoxon pareado para comparar pérdida de BA entre modelos
wilcox_result <- wilcox.test(
  robustez_pareado$CS_XGB,
  robustez_pareado$CI_SVM,
  paired = TRUE
)

# Imprimir resumen final
print("Resumen de diferencias promedio dBA:")
print(paste("mean_dBA CS ≈", round(mean(bootstrap_results_CS$mean_dBA), 3)))
print(paste("mean_dBA CI ≈", round(mean(bootstrap_results_CI$mean_dBA), 3)))
print("Test de Wilcoxon:")
print(wilcox_result)


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
    file = "session_info_estabilidad_modelos.txt", append = FALSE)

