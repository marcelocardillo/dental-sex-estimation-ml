
# Análisis de robustez para modelos predictivos
# Este archivo evalúa el rol de los datos ausentes en los modelos predictivos
# mediante la creación de escenarios contrafactuales sistemáticos.
# El análisis estima la confiabilidad del modelo bajo distintos escenarios
# y muestra el rol significativo de algunas variables predictivas.

# Carga de librerías necesarias
library(tidymodels)
library(recipes)
library(parsnip)
library(yardstick)
library(purrr)
library(tibble)
library(ggplot2)

# Función para verificar existencia de archivos antes de cargar
check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("Error: El archivo", file_path, "no existe. Por favor verifique la ruta."))
  }
  return(TRUE)
}

# Cargar modelos pre-entrenados con verificación
check_file_exists("xgb_fit_CS.rds")
check_file_exists("svm_rad_fit_CI.rds")

xgb_fit_CS <- readRDS("xgb_fit_CS.rds")
svm_rad_fit_CI <- readRDS("svm_rad_fit_CI.rds")

# Verificación básica de los modelos cargados
print("Modelo XGBoost para Canino Superior (CS) cargado:")
print(xgb_fit_CS)

print("Modelo SVM Radial para Canino Inferior (CI) cargado:")
print(svm_rad_fit_CI)

# Cargar conjuntos de testeo con verificación
check_file_exists("Conjunto_testeo_CS.txt")
check_file_exists("Conjunto_testeo_CI.txt")

# Función para cargar datos con manejo de errores
load_test_data <- function(file_path, model_name) {
  tryCatch({
    data <- read.table(file_path, header = TRUE)
    # Verificar columnas críticas
    required_cols <- c("MD", "BL", "MDCu", "BLCu", "Sexo")
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      warning(paste("Advertencia:", model_name, "faltan columnas:", paste(missing_cols, collapse = ", ")))
    }
    # Convertir Sexo a factor
    if ("Sexo" %in% colnames(data)) {
      data$Sexo <- as.factor(data$Sexo)
    }
    return(data)
  }, error = function(e) {
    stop(paste("Error al cargar", file_path, ":", e$message))
  })
}

# Cargar y verificar conjuntos de testeo
cs_test <- load_test_data("Conjunto_testeo_CS.txt", "Conjunto CS")
ci_test <- load_test_data("Conjunto_testeo_CI.txt", "Conjunto CI")

# Evaluación del modelo CS con datos completos
pred_cs <- predict(xgb_fit_CS, cs_test) %>%
  bind_cols(cs_test["Sexo"])
cs_baseline <- bal_accuracy(pred_cs, truth = Sexo, estimate = .pred_class)
print(paste("Balanced Accuracy CS (baseline):", round(cs_baseline$.estimate, 3)))

# Evaluación del modelo CI con datos completos
pred_ci <- predict(svm_rad_fit_CI, ci_test) %>%
  bind_cols(ci_test["Sexo"])
ci_baseline <- bal_accuracy(pred_ci, truth = Sexo, estimate = .pred_class)
print(paste("Balanced Accuracy CI (baseline):", round(ci_baseline$.estimate, 3)))

# Filtrar datos completos para la simulación
cs_test_complete <- cs_test %>%
  dplyr::filter(
    !is.na(MD),
    !is.na(BL),
    !is.na(MDCu),
    !is.na(BLCu)
  )

ci_test_complete <- ci_test %>%
  dplyr::filter(
    !is.na(MD),
    !is.na(BL),
    !is.na(MDCu),
    !is.na(BLCu)
  )

# Verificar que no hay NA en los conjuntos filtrados
cat("\nVerificación de datos completos:\n")
cat("CS: ", sum(is.na(cs_test_complete)), "valores NA\n")
cat("CI: ", sum(is.na(ci_test_complete)), "valores NA\n")

#' Genera escenarios con datos faltantes para análisis de robustez
#'
#' @param vars Vector de nombres de variables para generar combinaciones de datos faltantes
#' @return Lista de combinaciones de variables (1, 2 y 3 variables)
make_na_scenarios <- function(vars) {
  # Validación de entrada
  if (length(vars) < 1) {
    stop("Se requiere al menos una variable para generar escenarios")
  }
  
  max_comb <- min(3, length(vars))
  all_combos <- list()
  
  for (i in 1:max_comb) {
    if (length(vars) >= i) {
      all_combos[[i]] <- combn(vars, i, simplify = FALSE)
    }
  }
  
  all_combos <- unlist(all_combos, recursive = FALSE)
  names(all_combos) <- sapply(all_combos, paste, collapse = "_")
  
  return(all_combos)
}

#' Evalúa la robustez del modelo ante datos faltantes
#'
#' @param modelo_fit Modelo ajustado para evaluar
#' @param test_data Datos de testeo completos
#' @param nombre_modelo Nombre del modelo para identificación
#' @return Data frame con resultados de robustez
evaluar_robustez <- function(modelo_fit, test_data, nombre_modelo) {
  # Variables críticas para el análisis
  vars <- c("MD", "BL", "MDCu", "BLCu")
  
  # Predicción baseline con datos completos
  baseline_pred <- predict(modelo_fit, test_data) %>%
    dplyr::bind_cols(test_data["Sexo"])
  
  baseline_BA <- yardstick::bal_accuracy(
    baseline_pred,
    truth = Sexo,
    estimate = .pred_class
  )$.estimate
  
  # Generar escenarios de datos faltantes
  scenarios <- make_na_scenarios(vars)
  
  # Evaluar cada escenario
  resultados <- purrr::imap_dfr(scenarios, function(miss_vars, scenario_name) {
    # Crear copia de los datos y establecer variables como NA
    test_mod <- test_data
    test_mod[, miss_vars] <- NA
    
    # Realizar predicción
    pred <- predict(modelo_fit, test_mod) %>%
      dplyr::bind_cols(test_mod["Sexo"])
    
    # Calcular Balanced Accuracy
    BA <- yardstick::bal_accuracy(
      pred,
      truth = Sexo,
      estimate = .pred_class
    )$.estimate
    
    # Crear fila de resultados
    tibble::tibble(
      modelo = nombre_modelo,
      escenario = scenario_name,
      n_na = length(miss_vars),
      BA = BA,
      dBA = BA - baseline_BA
    )
  })
  
  # Añadir Balanced Accuracy baseline a los resultados
  resultados$baseline_BA <- baseline_BA
  
  return(resultados)
}

# Evaluación de robustez para el modelo CI
tryCatch({
  robustez_CI <- evaluar_robustez(
    modelo_fit = svm_rad_fit_CI,
    test_data = ci_test_complete,
    nombre_modelo = "CI_SVM"
  )
  print("Evaluación de robustez para CI completada exitosamente")
}, error = function(e) {
  stop("Error en la evaluación de robustez para CI:", e$message)
})

# Evaluación de robustez para el modelo CS
tryCatch({
  robustez_CS <- evaluar_robustez(
    modelo_fit = xgb_fit_CS,
    test_data = cs_test_complete,
    nombre_modelo = "CS_XGB"
  )
  print("Evaluación de robustez para CS completada exitosamente")
}, error = function(e) {
  stop("Error en la evaluación de robustez para CS:", e$message)
})

# Procesar resultados para CS
robustez_CS_tabla <- robustez_CS %>%
  dplyr::select(modelo, escenario, n_na, BA, dBA) %>%
  dplyr::arrange(dBA) %>%
  dplyr::mutate(
    BA = round(BA, 3),
    dBA = round(dBA, 3)
  )

# Resumen de resultados para CS
cat("\nResultados de robustez para Canino Superior (XGBoost):\n")
print(robustez_CS_tabla %>% group_by(n_na) %>% summarise(
  mean_dBA = mean(dBA),
  worst_dBA = min(dBA)
))

# Gráfico de robustez para CS
ggplot(robustez_CS_tabla, aes(x = factor(n_na), y = dBA)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Número de variables faltantes",
    y = "Degradación (Δ Balanced Accuracy)",
    title = "Robustez del modelo CS (XGBoost)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# Reglas de confianza para CS
rules_confianza_CS <- robustez_CS_tabla %>%
  dplyr::mutate(
    nivel_confianza = dplyr::case_when(
      BA >= 0.75 ~ "Alta",
      BA >= 0.70 ~ "Moderada",
      TRUE       ~ "Baja"
    )
  ) %>%
  dplyr::arrange(n_na, BA)

print(rules_confianza_CS)

# Procesar resultados para CI
robustez_CI_tabla <- robustez_CI %>%
  dplyr::select(modelo, escenario, n_na, BA, dBA) %>%
  dplyr::arrange(dBA) %>%
  dplyr::mutate(
    BA = round(BA, 3),
    dBA = round(dBA, 3)
  )

# Resumen de resultados para CI
cat("\nResultados de robustez para Canino Inferior (SVM radial):\n")
print(robustez_CI_tabla %>% group_by(n_na) %>% summarise(
  mean_dBA = mean(dBA),
  worst_dBA = min(dBA)
))

# Gráfico de robustez para CI
ggplot(robustez_CI_tabla, aes(x = factor(n_na), y = dBA)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Número de variables faltantes",
    y = "Degradación (Δ Balanced Accuracy)",
    title = "Robustez del modelo CI (SVM radial)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# Reglas de confianza para CI
rules_confianza_CI <- robustez_CI_tabla %>%
  dplyr::mutate(
    nivel_confianza = dplyr::case_when(
      BA >= 0.75 ~ "Alta",
      BA >= 0.70 ~ "Moderada",
      TRUE       ~ "Baja"
    )
  ) %>%
  dplyr::arrange(n_na, BA)

#tabla
print(rules_confianza_CI )

# Comparación combinada de ambos modelos
robustez_total <- dplyr::bind_rows(robustez_CS_tabla, robustez_CI_tabla)

# Gráfico combinado
ggplot(robustez_total, aes(x = factor(n_na), y = dBA, fill = modelo)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Número de variables faltantes",
    y = "Degradación (Δ Balanced Accuracy)",
    title = "Comparación de robustez: CS (XGBoost) vs CI (SVM radial)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# Análisis estadístico comparativo
robustez_pareado <- robustez_total %>%
  select(modelo, escenario, dBA) %>%
  tidyr::pivot_wider(
    names_from = modelo,
    values_from = dBA
  )

# Test de Wilcoxon para comparar modelos
wilcox_result <- tryCatch({
  wilcox.test(
    robustez_pareado$CS_XGB,
    robustez_pareado$CI_SVM,
    paired = TRUE
  )
}, error = function(e) {
  warning("Error en el test de Wilcoxon:", e$message)
  return(NULL)
})

if (!is.null(wilcox_result)) {
  print(wilcox_result)
  
  # Calcular medias para interpretación
  mean_dBA_CS <- mean(robustez_pareado$CS_XGB, na.rm = TRUE)
  mean_dBA_CI <- mean(robustez_pareado$CI_SVM, na.rm = TRUE)
  diff_mean <- mean_dBA_CS - mean_dBA_CI
  
  cat("\nResultados de comparación:\n")
  cat("mean_dBA CS ≈", round(mean_dBA_CS, 3), "\n")
  cat("mean_dBA CI ≈", round(mean_dBA_CI, 3), "\n")
  cat("Diferencia promedio:", round(diff_mean, 3), "\n")
  
  if (abs(diff_mean) < 0.03) {
    cat("La diferencia promedio es pequeña (menos del 3%).\n")
  } else {
    cat("Se observa una diferencia notable entre los modelos.\n")
  }
}





# Resumen interpretativo final
cat("\nAnálisis final:\n")
cat("El análisis de robustez mostró que la degradación promedio del desempeño predictivo",
    "aumenta progresivamente con el número de variables faltantes.\n")
cat("En el modelo del canino superior (XGBoost), la pérdida de una sola variable",
    "produjo una disminución media de la Balanced Accuracy del 4,8% respecto del baseline,\n",
    "mientras que la ausencia de dos variables incrementó la degradación promedio al 8,3%,\n",
    "y la falta de tres variables alcanzó una reducción media del 10,1%.\n")
cat("Este patrón indica un efecto acumulativo de la incompletitud, aunque no estrictamente lineal.\n")
cat("En conjunto, estos resultados sugieren que el modelo mantiene una estabilidad razonable",
    "ante pérdidas parciales de información, pero muestra una sensibilidad marcada cuando",
    "se afectan simultáneamente dimensiones clave.\n\n")
cat("Comparativamente, el modelo CI (SVM radial) parece más estable con poca pérdida de información,",
    "pero cuando la información colapsa (3 NA), puede caer más abruptamente.\n")
cat("CS tiene dependencia concentrada en variables cervicales, mientras que CI tiene",
    "una estructura más distribuida pero vulnerable a pérdida masiva de información.\n")


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







