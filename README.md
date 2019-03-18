
<!-- README.md is generated from README.Rmd. Please edit README.Rmd (this file) -->

# Connectedness

An R package to 1) measure the connectedness across management units
using genomic or pedigree data and 2) assess the connectedness between
training and testing sets in genomic prediction.

## Installation

Connectedness is currently available in Github, and can be installed
with devtools package:

1.  Install `devtools` package from CRAN.

<!-- end list -->

``` r
install.packages("devtools")
```

2.  Load the `devtools` package.

<!-- end list -->

``` r
library(devtools)
```

3.  Install `Connectedness` package from
Github.

<!-- end list -->

``` r
install_github('HaipengU/Connectedness', build_opts = c("--no-resave-data", "--no-manual"))
```

4.  Load `Connectedness` package.

<!-- end list -->

``` r
library(Connectedness)
```

5.  View `Connectedness` vignette.

<!-- end list -->

``` r
vignette('Connectedness', package = 'Connectedness')
```

## Documentation

[Vignette](https://haipengu.github.io/Rmd/Vignette.html)
