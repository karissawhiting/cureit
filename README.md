
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Cureit

<!-- badges: start -->
<!-- badges: end -->

This package allows you to easily fit and report results from cure
mixture models using a tidy framework. The package includes functions
to:

- Fit mixture cure models
- Summarize coefficients using tidiers and gtsummary
- Vizualize the data using nomograms
- Assess model results with Brier scores and K-indicies

## Installation

You can install {cureit} with the following code:

``` r
remotes::install_github("karissawhiting/cureit")
```

Load the package:

``` r
library(cureit)
```

# Fit Mixture Cure Models

Functions to fit the models are wrappers for the `smcure()` function
from the {smcure} package with the additional capability of passing a
formula and directly passing categortical variables without first
creating a model matrix:

``` r
p <- cureit(surv_formula = Surv(ttdeath, death) ~ age,
   cure_formula = ~ age,
   data = trial)
#> Warning: 0 of 100 did not converge.

p$surv_coefs
#> age, Survival model 
#>        -0.001010027

p$cure_coefs
#> (Intercept), Cure model         age, Cure model 
#>             -0.32294763              0.01067634
```

# Contributing

Please note that the cureit project is released with a [Contributor Code
of Conduct](http://www.karissawhiting.com/cureit/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

Thank you to [Sabrina Lin (@stl2137)](https://github.com/stl2137) for
package contributions!
