# Estimación sexual a partir de métricas dentales
# Sex Estimation from Dental Metrics
## Descripción general | Overview

### Español
Este repositorio contiene el flujo de análisis y los modelos de aprendizaje automático entrenados utilizados para la estimación del sexo en muestras bioarqueológicas a partir de métricas dentales.
El proyecto evalúa el desempeño y la robustez de distintos métodos de clasificación supervisada utilizando mediciones de caninos superiores e inferiores. El flujo de trabajo pone especial énfasis en la reproducibilidad, la prevención de fuga de información (*data leakage*) y la interpretabilidad de las predicciones de los modelos.

### English

This repository contains the analytical pipeline and trained machine learning models used for sex estimation in bioarchaeological samples based on dental metrics.
The project evaluates the performance and robustness of multiple supervised classification methods using measurements from upper and lower canines. The workflow emphasizes reproducibility, avoidance of data leakage, and interpretability of model predictions.

---

## Objetivos principales | Main objectives

### Español

* Estimar el sexo biológico a partir de métricas dentales.
* Comparar el desempeño predictivo de distintos algoritmos de aprendizaje automático.
* Evaluar la robustez de los modelos frente a escenarios de datos faltantes.
* Identificar las variables métricas dentales más informativas.
* Proporcionar un flujo de análisis reproducible para aplicaciones bioarqueológicas.

### English

* Estimate biological sex from dental metrics.
* Compare the predictive performance of multiple machine learning algorithms.
* Evaluate model robustness under missing-data scenarios.
* Identify the most informative dental variables.
* Provide a reproducible analytical workflow for bioarchaeological applications.

---

## Flujo de análisis | Analytical workflow

### Español

El flujo de análisis incluye:

1. Limpieza y validación de los datos.
2. Análisis exploratorio mediante análisis de componentes principales (PCA).
3. División de los datos en conjuntos de entrenamiento y prueba a nivel de individuo.
4. Imputación de datos faltantes mediante KNN, evitando la fuga de información.
5. Estandarización de las variables predictoras.
6. Validación cruzada estratificada.
7. Comparación del desempeño de los modelos.
8. Evaluación de la robustez frente a escenarios simulados de datos faltantes.
9. Interpretación de los modelos e identificación de la importancia de las variables.

### English

The analytical workflow includes:

1. Data cleaning and validation.
2. Exploratory analysis using principal component analysis (PCA).
3. Train/test splitting at the individual level.
4. KNN imputation of missing data while avoiding data leakage.
5. Predictor scaling.
6. Stratified cross-validation.
7. Model performance comparison.
8. Robustness analysis under simulated missing-data scenarios.
9. Model interpretation and assessment of variable importance.

---

## Modelos evaluados | Models evaluated

### Español

Se evaluó un conjunto diverso de modelos de clasificación, seleccionados para representar distintos enfoques metodológicos y supuestos estadísticos:

* Regresión logística (*Logistic Regression*)
* Análisis Discriminante Lineal (LDA)
* Elastic Net
* Máquinas de Vectores de Soporte (SVM), con kernels lineal y radial
* Random Forest
* XGBoost
* Naive Bayes

La combinación de estos algoritmos permite comparar modelos lineales y paramétricos, métodos de regularización, modelos capaces de representar relaciones no lineales e interacciones, y un enfoque generativo probabilístico. Todos los modelos se evaluaron bajo un mismo esquema de validación para permitir una comparación consistente de su desempeño predictivo.

### English

A diverse set of classification models was evaluated to represent different methodological approaches and statistical assumptions:

* Logistic Regression
* Linear Discriminant Analysis (LDA)
* Elastic Net
* Support Vector Machines (SVM), using linear and radial kernels
* Random Forest
* XGBoost
* Naive Bayes

This combination allows comparison of linear and parametric models, regularized methods, models capable of representing nonlinear relationships and interactions, and a probabilistic generative approach. All models were evaluated under the same validation framework to enable consistent comparison of predictive performance.

---

## Evaluación del desempeño | Performance evaluation

### Español

El desempeño predictivo de los modelos se evaluó principalmente mediante la exactitud balanceada (*Balanced Accuracy*), complementada por sensibilidad y especificidad.

La exactitud balanceada permite otorgar igual peso al desempeño de ambas categorías y proporciona una medida simétrica de la capacidad de clasificación. Los modelos fueron evaluados mediante validación cruzada estratificada, preservando la distribución de las categorías durante el procedimiento de evaluación.

Las diferencias entre modelos se evaluaron utilizando las estimaciones obtenidas en las mismas particiones de validación cruzada, permitiendo comparar directamente su desempeño bajo condiciones equivalentes.

### English

Model predictive performance was primarily evaluated using Balanced Accuracy, complemented by sensitivity and specificity.

Balanced Accuracy gives equal weight to the performance of both classes and provides a symmetric measure of classification performance. Models were evaluated using stratified cross-validation, preserving the class distribution throughout the evaluation procedure.

Differences between models were assessed using estimates obtained from the same cross-validation resamples, allowing their performance to be directly compared under equivalent conditions.

---

## Robustez frente a datos faltantes | Robustness to missing data

### Español

La robustez de los modelos se evaluó mediante escenarios simulados de información incompleta. Estos análisis permitieron examinar cómo cambia el desempeño predictivo cuando algunas de las métricas dentales no están disponibles.

La imputación de valores faltantes se realizó mediante KNN dentro del flujo analítico, evitando la utilización de información del conjunto de prueba durante el entrenamiento y reduciendo así el riesgo de fuga de información.

### English

Model robustness was evaluated under simulated missing-data scenarios. These analyses examined how predictive performance changes when some dental measurements are unavailable.

Missing values were imputed using KNN within the analytical workflow, avoiding the use of information from the test set during model training and thereby reducing the risk of data leakage.

---

## Interpretación de los modelos | Model interpretation

### Español

Además de evaluar el desempeño predictivo, se analizaron la importancia de las variables y la contribución de los predictores a las decisiones de clasificación. Para XGBoost se utilizó el método SHAP (*SHapley Additive exPlanations*) para caracterizar la contribución individual de las variables a las predicciones del modelo.

### English

In addition to predictive performance, variable importance and predictor contributions to classification decisions were assessed. For XGBoost, SHAP (*SHapley Additive exPlanations*) was used to characterize the individual contribution of predictors to model predictions.

---

## Reproducibilidad | Reproducibility

### Español

Todos los análisis fueron realizados en **R**. El repositorio contiene los scripts y recursos asociados al flujo analítico y los modelos de aprendizaje automático utilizados en el estudio.

El análisis fue diseñado para minimizar el riesgo de fuga de información y mantener una separación entre los procedimientos de entrenamiento, validación y prueba.

### English

All analyses were conducted in **R**. The repository contains the scripts and associated resources for the analytical workflow and machine learning models used in the study.

The analytical workflow was designed to minimize the risk of data leakage and maintain separation between training, validation, and testing procedures.

---

## Estructura del repositorio | Repository structure

La organización del repositorio sigue la lógica general del flujo analítico:

```text
data/       Datos utilizados en los análisis
scripts/    Scripts del flujo analítico
models/     Modelos de aprendizaje automático entrenados
results/    Resultados y salidas del análisis
```

The repository is organized according to the general analytical workflow:

```text
data/       Data used in the analyses
scripts/    Analytical workflow scripts
models/     Trained machine learning models
results/    Analysis outputs and results
```

> **Nota:** Los nombres y contenidos de las carpetas deben mantenerse sincronizados con la estructura efectiva del repositorio.

---

## Información sobre el estudio | Study information

### Español

Este repositorio acompaña el estudio sobre estimación sexual mediante métricas de caninos permanentes y proporciona los recursos necesarios para documentar y reproducir el flujo de análisis computacional.

### English

This repository accompanies the study on sex estimation using permanent canine metrics and provides the resources required to document and reproduce the computational analytical workflow.

---

## Licencia y citación | License and citation

Para información sobre la licencia, versión del conjunto de datos y forma de citar este repositorio, consulte los metadatos y la versión correspondiente depositada en Zenodo.

For information on licensing, dataset version, and how to cite this repository, please refer to the metadata and corresponding version deposited in Zenodo.
