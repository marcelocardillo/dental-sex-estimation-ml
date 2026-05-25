Overview

This repository contains the analytical pipeline and trained machine learning models used for sex estimation in bioarchaeological samples based on dental metrics.

The project evaluates the performance and robustness of multiple supervised classification methods using measurements from upper and lower canines. The workflow emphasizes reproducibility, avoidance of data leakage, and interpretability of model predictions.

Main objectives
Estimate biological sex from dental metrics
Compare multiple machine learning algorithms
Evaluate robustness under missing-data scenarios
Identify the most informative dental variables
Provide a reproducible analytical workflow for bioarchaeological applications
Methods

The analytical workflow includes:

Data cleaning and validation
PCA exploratory analyses
Train/test split by individual
KNN imputation without data leakage
Feature scaling
Stratified cross-validation
Model comparison using:
Logistic Regression
LDA
Elastic Net
SVM
Random Forest
XGBoost
Naive Bayes
Performance evaluation using:
Balanced Accuracy
Sensitivity
Specificity

Reproducibility: All analyses were conducted in R.


Robustness analysis under simulated missing-data scenarios
Variable importance and model interpretation
