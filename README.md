
# localPCA

<!-- badges: start -->
<!-- badges: end -->

localPCA is an R-package that uses the R-package [lostruct](https://github.com/petrelharp/local_pca/tree/master) to perform genome-wide local PCA scans (following [Li & Ralph 2019](https://doi.org/10.1534/genetics.118.301747) method) through a reproducible pipeline. The pipeline is divided into three main functions in the R-package: 
```r
# actual lostruct() analysis run across a defined window-size.
loPCA()

# Extract significantly divergent genomic windows (from genome-wide pattern)
# and cluster them into outlier genomic regions.
getRegion()

# Perform K-means clustering and compute Ho (observed hterozygosity) for each outlier genomic regions.
# Output summary table and plots for visual inspection of candidate inverted regions.
detectInv()
```

## Installation

You can install the development version of localPCA either directly in R:

``` r
devtools::install_github("stella-huynh/localPCA")
```
or by cloning the repository:
``` r
git clone git@github.com:stella-huynh/localPCA.git #local bash command-line
devtools::build("<path/to/package>") #in R
```

## Example dataset

This package comes with example results from the three main functions loPCA(), getRegion() and detectInv(), respectively accessible as follows:

``` r
library(localPCA)
data("winList")
data("regList")
data("invList")
```

