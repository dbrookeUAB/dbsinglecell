
# dbsinglecell

<!-- badges: start -->
[![R-CMD-check](https://github.com/dbrookeUAB/dbsinglecell/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dbrookeUAB/dbsinglecell/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

A collection of functions for processing single-cell RNAseq data that I am using constantly. I created this package to make these functions portable for myself. Use at your own risk.  

## Installation

If you so desire to use this package, install  by using

``` r
remotes::install_github("dbrookeUAB/dbsinglecell")
```

To use `HDSCAN` or `umap` functions, you must install their respective python libraries by

``` r
library(dbsinglecell)
install_python_packages()
```
Cheers! 
