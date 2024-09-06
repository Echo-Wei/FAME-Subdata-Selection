# FAME-Subdata-Selection

Reference: An Efficient Filtering Approach for Model Estimation in Sparse Regression

## Introduction
This package contains several filtering approach for model estimation. 
1. Dataset: Generating sample datasets that can be used as input for functions 
2. Function: 
    - Proposed and existing filtering and estimating approach for subdata selection
    - Implementation file with example on dataset simulation and application of above approaches
    - Score index calculation for input datasets with continuous predictors and continuous or binary response varaible

## Example
```
x <- matrix(rnorm(2000 * 100), nrow = 2000)
y <- rnorm(2000)
FAME_fit <- FAME(0.7, x, y, 200, 5, random_seed, index_beta, tr)
FAME_MSE <- FAME_fit$mse
```

