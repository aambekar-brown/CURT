# CURT

## Description
This R package implements the censoring unbiased random forest algorithm developed in Steingrimsson et al. (2019) using both a Buckley-James loss function and a doubly robust loss function. The outputs are survival probabilities on a test set and variable importance measures were higher values imply that the variable is more important.

Reference: Steingrimsson, J. A., Diao, L., & Strawderman, R. L. (2019). Censoring unbiased regression trees and ensembles. *Journal of the American Statistical Association*, 114(525), 370-383.

## Installation

Install CURT directly from GitHub using:

# If devtools is not already installed
if(!require(devtools)) install.packages("devtools")
devtools::install_github("aambekar-brown/CURT")

## Usage

```R

# Load the CURT package
library(CURT)

# Access the function documentation
?curt

# Example function call
curt(
  formula,
  train_data,
  test_data,
  n.tree = 1000,
  tau1 = NULL,
  type = "drl",
  mtry = NULL
)

```

### Arguments

**formula**	an object of class "formula", e.g., `Surv(obs, delta) ~ .`.

**train_data** A dataframe containing the training dataset.

**test_data** A dataframe containing the test dataset for which survival predictions will be made for.

**n.tree** Number of trees used in the forest

**tau1** At what time-point survival predictions should be calculated (e.g, if time is measured in years and tau=1 the algorithm would predict one year survival probabilities).

**type** The type of loss function to be used in the random forest algorithm. Type = "bjl" uses the Buckley-James loss function and is the default and type= "drl" uses the doubly robust loss function.

**mtry** Number of variables randomly sampled as candidates at each split.

### Return Value 
Returns a list containing survival predictions on the test set and variable importance measures

## Examples

### Example 1: Using Built-in Data
```R
train_data <- read.csv(system.file("extdata", "train_data.csv", package = "CURT"))
test_data <- read.csv(system.file("extdata", "test_data.csv", package = "CURT"))
result <- curt(Surv(obs, delta) ~ ., train_data, test_data, n.tree = 500, tau1 = 0.6, type = "drl")
```

### Example 2: Importing Data from Local Drive
```R
# Specify the correct paths to your data files
train_data_path <- "path/to/your/train_data.csv"
test_data_path <- "path/to/your/test_data.csv"
train_data <- read.csv(train_data_path)
test_data <- read.csv(test_data_path)

# Apply the CURT function
result <- curt(Surv(obs, delta) ~ ., train_data, test_data)
```
