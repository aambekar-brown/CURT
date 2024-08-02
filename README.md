# CURT

## Description
This R package implements the censoring unbiased random forest algorithm developed in Steingrimsson et al. (2019) using both a Buckley-James loss function and a doubly robust loss function. The outputs are survival probabilities on a test set and variable importance measures were higher values imply that the variable is more important.

Steingrimsson, J. A., Diao, L., & Strawderman, R. L. (2019). Censoring unbiased regression trees and ensembles. *Journal of the American Statistical Association*, 114(525), 370-383.

## Installation

To install the latest version of CURT from GitHub, use:
```R
# If devtools is not already installed
if(!require(devtools)) install.packages("devtools")
devtools::install_github("aambekar-brown/CURT")