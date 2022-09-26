# popdynmodel

<!-- badges: start -->
<!-- badges: end -->
*popdynmodel* provides tools to write and implement population dynamics models in a hierarchical Bayesian framework. Models are intended for occupancy, abundance and/or biomass time series and allow for estimations of distribution change rates and population growth rates at various spatial and temporal scales. This package has two main functions : 
* Item *popbay_model* writes population dynamics models, selects and prepares input model data, and fits models to selected data 
* Item *popbay_data* checks and standardises a dataset in the format requested by *popbay_model*

*popdynmodel* was initially developed to analyse the river fish monitoring data from the ASPE database managed by the French Biodiversity Agency (OFB), but it can be used for a wide variety of taxa.

## Installation

*popdynmodel* requires the following packages to be installed: *rjags*, *coda*, *R2jags*, *stringr*, *magrittr* and *dplyr*. You can install the missing packages from CRAN using the following command :

``` r 
install.packages(c("coda","rjags","R2jags","stringr","magrittr","dplyr"))
```

*popdynmodel* uses the interface from R to JAGS library provides by the *rjags* package to fit models. JAGS is a separate program that must be independently installed ([https://mcmc-jags.sourceforge.io/](https://mcmc-jags.sourceforge.io/)). Then you can install the development version of *popdynmodel* from Github using :

``` r
library(devtools)
install_github("manue6/popdynmodel")
```

You can also install from the downloadable tar archive [*popdynmodel_1.0.tar.gz*](https://github.com/manue6/popdynmodel_archive) using :

``` r
install.packages("popdynmodel_1.0.tar.gz", repos = NULL)
```
