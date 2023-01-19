
# changomics

<!-- badges: start -->
<!-- badges: end -->

The goal of changomics is to correct a dataset of metabolites that has been distorted during the 
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
corrected<-changomics(as.data.frame(met.df),graph=FALSE)

```
There is a parameter "graph" which allow to plot the change points found in the data sequence and show how it is detrended. When set to TRUE it will produce a plot for the first metabolite of the data set to be corrected.
``` r
library(changomics)
# Load the data of metabolites
data(met.df)
# Use changomics to correct the metabolites
met1<-as.data.frame(met.df$met1)
corrected<-changomics(met1,graph=TRUE)

```
Here we have already a metabolite that is not distorted no the function changomics() is not really correcting it.
``` r
library(changomics)
# Load the data of metabolites
data(met.df)
# Use changomics to correct the metabolites
met2<-as.data.frame(met.df$met2)
corrected<-changomics(met2,graph=TRUE)

```
This time the signal was distorted and we see the change points estimated and the correction applied.
