
# TRES

<!-- badges: start -->
<!-- badges: end -->

The package **TRES** implements the least squares and envelope estimation under the framework of tensor regression models. The general model-free envelope models can also be flexibly handled by the package via three types of envelope estimation algorithms: 
- Full Grassmannian (FG) algorithm.
- 1D algorithm.
- Envelope coordinate descent (ECD) algorithm
- Partial least squares (PLS) type algorithm.

## Installation

You can install the released version of TRES from [CRAN](https://CRAN.R-project.org) with:

``` r
# Install the latest released version from CRAN
install.packages("TRES")

# Or the development version from GitHub:
remotes::install_github("jerryfsu3333/TRES")
```

## Example

This is a basic example which shows you how to use function `TRR.fit` in Tensor Response Regression (TRR) model with least square.

``` r
library(TRES)
## Load data "bat"
data("bat")
x <- bat$x
y <- bat$y

## Fitting with OLS and 1D envelope method.
fit_ols <- TRR.fit(x, y, method="standard")
fit_1D <- TRR.fit(x, y, u = c(14,14), method="1D") # pass envelope rank (14,14)

## Print cofficient
coef(fit_1D)

## Print the summary
summary(fit_1D)

## Make the prediction on the original dataset
predict(fit_1D, x)

## Draw the plots of two-way coefficient tensor (i.e., matrix) and p-value tensor.
plot(fit_ols)
plot(fit_1D)
```

 <!-- badges: start -->
  [![Travis build status](https://travis-ci.org/jerryfsu3333/TRES.svg?branch=master)](https://travis-ci.org/jerryfsu3333/TRES)
  <!-- badges: end -->

