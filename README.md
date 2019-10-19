
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
# Install devtools from CRAN
install.packages("TRES")

# Or the development version from GitHub:
devtools::install_github("jerryfsu3333/TRES_code")
```

## Example

This is a basic example which shows you how to use function `TRR` in Tensor Response Regression (TRR) model with least square.

``` r
library(TRES)
## Load data "bat"
data("bat")
Xn <- bat$Xn
Yn <- bat$Yn
fit <- TRR(Xn, Yn, method="standard")

## Print cofficient
coef(fit)

## Print the summary
summary(fit)

## Make the prediction on the original dataset
predict(fit, Xn)

## Draw the plot of two-way coefficient tensor (or matrix)
plot(fit, ask=FALSE)
```

 <!-- badges: start -->
  [![Travis build status](https://travis-ci.org/jerryfsu3333/TRES_code.svg?branch=master)](https://travis-ci.org/jerryfsu3333/TRES_code)
  <!-- badges: end -->

