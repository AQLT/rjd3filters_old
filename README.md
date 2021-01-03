
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rjdfilters

rjdfilters is R package on linear filters. It allows to create symmetric
and asymmetric moving averages with:

-   local polynomial filters, as defined by Proietti and Luati (2008);

-   the FST approach of Grun-Rehomme, Guggemos, and Ladiray (2018),
    based on the optimization of the three criteria Fidelity, Smoothness
    and Timeliness.

Some quality criteria defined by Wildi and McElroy (2019) can also be
computed.

## Installation

rjdfilters relies on the
[rJava](https://CRAN.R-project.org/package=rJava) package and Java SE 8
or later version is required.

``` r
# Install development version from GitHub
# install.packages("devtools")
devtools::install_github("palatej/rjdfilters")
```

If you have troubles with the installation, check the [installation
manual](https://github.com/jdemetra/rjdemetra/wiki/Installation-manual)
of RJDemetra.

## Basic example

To seasonally adjust a time series with a pre-defined specification you
can either use the `x13()` function for the X-13ARIMA method or the
`tramoseats()` function for the TRAMO-SEATS method.

``` r
library(rjdfilters)
#> Java requirements fullfilled, found version 11.0.5
y <- window(retailsa$WomensClothingStores,start = 2000)
```

## Bibliography

Grun-Rehomme, Michel, Fabien Guggemos, and Dominique Ladiray (2018).
“Asymmetric Moving Averages Minimizing Phase Shift”. In: *Handbook on
Seasonal Adjustment*. URL:
<https://ec.europa.eu/eurostat/web/products-manuals-and-guidelines/-/KS-GQ-18-001>.

Proietti, Tommaso and Alessandra Luati (Dec. 2008). “Real time
estimation in local polynomial regression, with application to
trend-cycle analysis”. In: *Ann. Appl. Stat.* 2.4, pp. 1523–1553. DOI:
10.1214/08-AOAS195. URL: <https://doi.org/10.1214/08-AOAS195>.

Wildi, Marc and Tucker McElroy (2019). “The trilemma between accuracy,
timeliness and smoothness in real-time signal extraction”. In:
*International Journal of Forecasting* 35.3, pp. 1072–1084. URL:
<https://EconPapers.repec.org/RePEc:eee:intfor:v:35:y:2019:i:3:p:1072->
1084.
