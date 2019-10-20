## Test environments
* local OS X install, R 3.5.2 and devel
* win-builder (devel, release and oldrelease)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking dependencies in R code ... NOTE

installed size is  6.9Mb
sub-directories of 1Mb or more:
data   6.6Mb

The dataset square.rda in /data consists of a 64 by 64 by 200 tensor. Users can use this dataset check the functions very quickly, instead of build the envelope-structure tensor date from the scratch.

## Downstream dependencies
No reverse dependence packages.

## Updates
* This is the updated version of old package TRES_0.1.0. Since we made many important updates of the package, like rewrite some functions into S3 methods, deprecated some functions, we use the version number 1.0.0.
* We host the project on github: https://github.com/jerryfsu3333/TRES
* The new README.md and NEWS.md is included.