# SolarClusGnr
Solar irradiance time-series clustering and down-scaling

``` SolarClusGnr ``` is an R-package allows a reproducible research for nonparametric clustering and down-scaling of daily solar irradiation time-series. This currently version includes: 

   ``` SIR_Data ``` constructor function of objects of SIRData class. Once the user creates a 'SIRData' object from his data, he no longer need other inputs, ALL the rest work will done automatically.
   
   ``` DPGMMclus ``` S3 Method for nonparametric Bayesian Dirichlet-Gaussian mixture model clustering of daily clearness index distributions. It can be also used to perform any data clustering of class matrix other than irradiance data. It generate an object of class 'clusData' containing the clustering outputs.
   
   ``` parClusGena ``` constructor function of objects of 'genData' class, needed for the generation of hight resolution solar irradiance data.
   
   ``` GenData ``` function to generate high resolution solar irradiance time-series. It requires object of 'genData' class as input.
   
   ``` clPlot ``` function to generate plots of the resulting classes.
   
#### For any queries, please e-mail: azeddine.frimane@yahoo.com, I will be happy to answer your questions.
   
# Installation

You can install ``` SolarClusGnr ``` with the [remotes](https://install-github.me/r-lib/remotes) package:

```
remotes::install_github("frimane/SolarClusGnr")
```

or with [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

```
devtools::install_github("frimane/SolarClusGnr")
```

# License

This package is free and open source software under MIT-license.
