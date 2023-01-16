
# changomics

<!-- badges: start -->
<!-- badges: end -->

The goal of changomics is to normalize a dataset of metabolites that has been distorted during the 
data gathering. It detects automatically through dynamic programming the change points in time series and then normalize, residualize and detrend the time series of metabolites.

## Installation

You can install the development version of changomics from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("antoinejeanrenaud/changomics")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(changomics)
# Load the data of metabolites
data(met.df)
# Use changomics to correct the metabolites
corrected<-changomics(as.data.frame(met.df))

```

