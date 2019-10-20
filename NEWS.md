# TRES 1.0.0

## Major changes
- Switch the arguments `Xn` and `Yn` in functions: `TensEnv_dim`, `TensPLS_cv2d3d`, `Tenv`, `Tenv_Pval`, `TPR` and `TRR`.
- Removed package `mvtnorm` in `DESCRIPTION`, it is not used in the package.
- Move all the packages in `Depends` in `DESCRIPTION` to `Imports` except for `ManifoldOptim`.
- Add examples to some functions: `PMSE`, `TensEnv_dim`, `TensPLS_cv2d3d`, `Tenv_Pval`.
- Add simulated data "bat" and "square" with which users can quickly verify the functions in the package.
```
# Usage
> data("bat")
> data("square")
```
- `TPR`, `TRR`: 
  - Change the output names of some functions. E.g., `Bhat` in the object returned from `TPR` and `TRR` is renamed as `coefficients`.
  - Append `fitted.values` and `residuals` into the output.
  - More robust to the input type: `Yn` can be vector or matrix and `Xn` can be matrix, array or tensor for `TPR`, `Xn` can be vector or
  matrix and `Yn` can be matrix, array or tensor for `TRR`.
  - Incorporate `FG_TRR` and `FG_TPR` into `TRR` and `TPR`. Add one more option for `method`, "FG". Also add argument `Gamma_init` for "FG" method. 
- `Tenv_Pval`: Change the data type of outputs from array to tensor. Change the name `se_mat` to `se`.
- `PMSE`: Rewrite `PMSE`. `Xn` can be matrix, array or tensor, `Yn` can be vector or matrix, and `Bhat` can be vector, matrix, array or tensor
as long as the dimensions match the ones of `Xn` and `Yn`.
- `TRR_sim`, `TPR_sim`: Add two simulation functions to generate data used in TRR and TPR which can help user quickly test the functions.


## S3 methods
Construct S3 object for `TRR` and `TPR` with `class` attribute "Tenv"
```
> data("bat")
> fit <- TRR(bat$Xn, bat$Yn, method = "standard")
> class(fit)
> [1] "Tenv"
```
- `print.Tenv`: Prints the coefficients from `TPR` and `TRR`
- `predict.Tenv`: Make predictions of new data.
- `summary.Tenv`: Append dimensions of X,Y, sample size, mse, p_val and s.e. to the output object from `TPR` and `TRR`.
- `print.summary.Tenv`: Print call, dimensions of X, Y, sample size, mse, the coefficient and p_value. (invoked implicitly when there is 
no assignment of `summary.Tenv`)
- `plot.Tenv`: Draw the plot of coefficients from `TRR` and `TPR`, and draw p_val plot from `TRR`. 
- `vcov`: No covariance for coefficients. But for `TRR`, we print the standard error for coefficients.
- `fitted.default`: Calculates the fitted Y for `TPR` and `TRR` separately.
- `residuals.default`: Calculate Y minus fitted Y for `TPR` and `TRR`
- `coef.default`: Print the coefficients for `TPR` and `TRR`.

## Deprecated functions
- `FG_TPR`, `FG_TRR`
>  Help pages for deprecated functions are available at `help("TRES-deprecated")`.

## Bug fixed
- Default value `u=NULL` for `TPR` and `TRR`.
- Default value `maxdim=10` and `nfold=5` for `TensPLS_cv2d3d`
- `PMSE`: now accept tensor `Xn`
