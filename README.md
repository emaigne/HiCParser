
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HiCParser

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/emaigne/HiCParser)](https://github.com/emaigne/HiCParser/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/emaigne/HiCParser)](https://github.com/emaigne/HiCParser/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/HiCParser.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/HiCParser)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/HiCParser.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/HiCParser)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/HiCParser.svg)](http://bioconductor.org/packages/stats/bioc/HiCParser/)
[![Bioc
support](https://bioconductor.org/shields/posts/HiCParser.svg)](https://support.bioconductor.org/tag/HiCParser)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/HiCParser.svg)](https://bioconductor.org/packages/release/bioc/html/HiCParser.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/HiCParser.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/HiCParser/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/HiCParser.svg)](https://bioconductor.org/packages/release/bioc/html/HiCParser.html#since)
[![R-CMD-check-bioc](https://github.com/emaigne/HiCParser/actions/workflows/R-CMD-check-bioc.yaml/badge.svg)](https://github.com/emaigne/HiCParser/actions/workflows/R-CMD-check-bioc.yaml)
[![Codecov test
coverage](https://codecov.io/gh/emaigne/HiCParser/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/emaigne/HiCParser?branch=devel)
<!-- badges: end -->

The goal of `HiCParser` is to …

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `HiCParser` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("HiCParser")
```

And the development version from
[GitHub](https://github.com/emaigne/HiCParser) with:

``` r
BiocManager::install("emaigne/HiCParser")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("HiCParser")
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub!

## Citation

Below is the citation output from using `citation('HiCParser')` in R.
Please run this yourself to check for any updates on how to cite
**HiCParser**.

``` r
print(citation('HiCParser'), bibtex = TRUE)
#> Pour citer le package 'HiCParser' dans une publication, utilisez :
#> 
#>   Maigné E, Zytnicki M (2024). _Parse HiC data to import them in R_.
#>   doi:10.18129/B9.bioc.HiCParser
#>   <https://doi.org/10.18129/B9.bioc.HiCParser>,
#>   https://github.com/emaigne/HiCParser/HiCParser - R package version
#>   0.1.0, <http://www.bioconductor.org/packages/HiCParser>.
#> 
#> Une entrée BibTeX pour les utilisateurs LaTeX est
#> 
#>   @Manual{,
#>     title = {Parse HiC data to import them in R},
#>     author = {Elise Maigné and Matthias Zytnicki},
#>     year = {2024},
#>     url = {http://www.bioconductor.org/packages/HiCParser},
#>     note = {https://github.com/emaigne/HiCParser/HiCParser - R package version 0.1.0},
#>     doi = {10.18129/B9.bioc.HiCParser},
#>   }
```

Please note that the `HiCParser` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `HiCParser` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.17/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://emaigne.github.io/HiCParser) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.17/biocthis)*.
