## Test environments
* macOS 10.13.6 High Sierra, R 3.6.3, R 4.0.2
* Ubuntu 16.04.6 LTS (on travis-ci), (devel, release and oldrelease)
* win-builder (devel, release and oldrelease)

## R CMD check results
There were no ERRORs and WARNINGs. But there is a NOTE as following:

* New submission. Package was archived on CRAN. CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2021-05-05 as requires archived package
    'rTensor'.

`TRES 1.1.3` was archived on 2021-05-05 as the import package `rTensor` was archived. Since the `rTensor` passes its build checks on CRAN and was released on 2021-05-15, I aim to submit `TRES` package to bring it back to life.

## Downstream dependencies
No reverse dependence packages.

## Updates
* This is the updated version of old package TRES 1.1.3. 