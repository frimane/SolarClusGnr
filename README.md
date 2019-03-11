[![Build Status](https://travis-ci.org/frimane/SolarClusGnr.svg?branch=master)](https://travis-ci.org/frimane/SolarClusGnr)

# SolarClusGnr: Solar irradiance time-series clustering and down-scaling

## `SolarClusGnr `
It is an R-package allows a reproducible research for non-parametric clustering and downscaling of daily solar irradiation time-series. The current version includes: 

   - `SIR_Data` Constructor function of objects of `SIRData` class. Once the user creates a `SIRData` object from his data, he no longer need other inputs, **ALL** the rest work will done automatically.
   
  - `DPGMMclus` S3 Method for non-parametric Bayesian Dirichlet-Gaussian mixture model clustering of daily clearness index distributions. It can be also used to perform any data clustering of class matrix other than irradiance data. It generate an object of class `clusData` containing the clustering outputs.
   
  - `parClusGen` Constructor function of objects of `genData` class, needed for the generation of hight resolution solar irradiance data.
   
   - `GenData` Function to generate high resolution solar irradiance time-series. It requires object of `genData` class as input.
   
   - `clPlot` Function to generate plots of the resulting classes.
   
## Further inquiry

The author is happy to answer your questions and is open to future collaborations on this topic.
Please contact: azeddine.frimane@yahoo.com or azeddine.frimane@uit.ac.ma.
   
## Installation

Users can install the development version of `SolarClusGnr`:

with either the [remotes](https://install-github.me/r-lib/remotes) package:

```
remotes::install_github("frimane/SolarClusGnr")
```

or with the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package:

```
devtools::install_github("frimane/SolarClusGnr")
```

## Citation

The original paper that describes the methods implemented in `SolarClusGnr` is:

Frimane, Â., Soubdhan, T., Bright, J.M., Aggour, M., 2019. Nonparametric bayesian-based recognition of solar irradiance conditions: Application to the generation of high temporal resolution synthetic solar irradiance data. Solar Energy 182, 462-479. URL:http://www.sciencedirect.com/science/article/pii/S0038092X19301781, doi:https://doi.org/10.1016/j.solener.2019.02.052.

The BibTex entry:

@article{Frimane2019,
title = "Nonparametric Bayesian-based recognition of solar irradiance conditions: Application to the generation of high temporal resolution synthetic solar irradiance data",
journal = "Solar Energy",
volume = "182",
pages = "462-479",
year = "2019",
issn = "0038-092X",
doi = "https://doi.org/10.1016/j.solener.2019.02.052",
url = "http://www.sciencedirect.com/science/article/pii/S0038092X19301781",
author = "{\^A}zeddine Frimane and Ted Soubdhan and Jamie M. Bright and Mohammed Aggour",
keywords = "Solar irradiance, Clustering, Clearness index, Bayesian nonparametric, Synthetic irradiance",
}

## License

This package is free and open source software under MIT-license.
