# SolarClusGnr
Solar irradiance time-series clustering and down-scaling

``` SolarClusGnr ``` is an R-package allows a reproducible research for nonparametric clustering and down-scaling of daily solar irradiation time-series. This currently version includes a set of generic functions and S3 methods: 
   ``` SIR_Data ``` a constructor function of objects of SIRData class. Once the user creates a SIRData object, he no longer need other inputs, ALL the rest work will done automatically.
   ``` DPGMMclus ``` a S3 Method for nonparametric Bayesian Dirichlet-Gaussian mixture model clustering of daily clearness index distributions. It can be also used to perform any data clustering of class matrix other than irradiance data of class SIRData.
   ``` parClusGena ``` a constructor function of objects of genData class, needed to perform the generation procedure of hight resolution solar irradiance data.
   ``` GenData ``` a function to generate high resolution solar irradiance data. 
   ``` clPlot ``` a function to generate plot of classes.
   ``` extEHIcalcGeneric ``` a function to generate the extraterrestrial horizontal solar irradiance data at any time resolution.
   
# Installation

You can install ``` SolarClusGnr ``` with the [remotes](https://install-github.me/r-lib/remotes) package:

```
remotes::install_github("frimane/SolarClusGnr")
```

or with [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

```
devtools::install_github("frimane/SolarClusGnr")
```

## License

This package is free and open source software, licensed under GPLv3.

