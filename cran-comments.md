## Test environments
* macOS 10.13.6 High Sierra, R 3.6.3, R 4.0.5, R 4.1.0
* Ubuntu 16.04.6 LTS (on travis-ci), (devel, release and oldrelease)
* win-builder (devel, release and oldrelease)

## R CMD check results
There were no ERRORs and WARNINGs. But there is a NOTE as following:

- New submission. Package was archived on CRAN. CRAN repository db overrides:

- Possibly mis-spelled words in DESCRIPTION:
  - ECD (10:464, 10:939)
  - FG (10:286)
  - Grassmannian (10:272)
  - SIMPLS (10:557)
  - TPR (10:110, 10:763)
  - TRR (10:72, 10:681)
  - Zhang (10:696, 10:772, 10:872, 10:967)

- X-CRAN-Comment: Archived on 2021-05-05 as requires archived package
    'rTensor'.

**Responses**:
- `TRES 1.1.3` was archived on 2021-05-05 as the import package `rTensor` was archived. Since the `rTensor` passes its build checks on CRAN and was released on 2021-05-15, we aim to submit `TRES` package to bring it back online.
- We checked the words and verified that all the spellings are legitimate.

## Downstream dependencies
No reverse dependence packages.

## Updates
This is the updated version of old package TRES 1.1.3. 