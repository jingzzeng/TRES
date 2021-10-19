## Test environments
* macOS 10.13.6 High Sierra, R 3.6.3, R 4.1.0
* Ubuntu 16.04.6 LTS (on travis-ci), (devel, release and oldrelease)
* win-builder (devel, release and oldrelease)

## R CMD check results
There were no ERRORs and WARNINGs. But there is a NOTE as following:
* Possibly misspelled words in DESCRIPTION:

  Zeng (10:1045)

*Respones*: The spell is correct.

* Found the following (possibly) invalid DOIs:

  DOI: 10.18637/jss.v099.i12

  From: DESCRIPTION
  
  inst/CITATION

  Status: Not Found

  Message: 404

*Response*: The DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN.


## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

## Updates
This is the updated version of old package TRES 1.1.4. 