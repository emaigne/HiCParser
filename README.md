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

The goal of `HiCParser` is to parse Hi-C data (`HiCParser` supports serveral formats), and import them in R, as an `InteractionSet` object.

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

## Supported formats

So far, `HiCParser` supports:

  - [cool and mcool](https://github.com/open2c/cooler) formats
  - [hic](https://github.com/aidenlab/hictools) format
  - [HiC-Pro](https://github.com/nservant/HiC-Pro) format
  - A tabular format, where

      - the first column is named "chromosome"
      - the second column is named "position 1" or "position.1"
      - the third column is named "position 2" or "position.2"
      - the fourth column is named "*x*.R*y*", and *x* is the id of the condition ("1", or "2", usually), *y* is the id of the replicate ("1", "2", "3", etc.)

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("HiCParser")
## basic example code
```

## Citation

Below is the citation output from using `citation('HiCParser')` in R.
Please run this yourself to check for any updates on how to cite
**HiCParser**.

To cite the ‘HiCParser’ HiCParser in a publication, use :

  Maigné E, Zytnicki M (2024). _A multiple format Hi-C data parser_.
  doi:10.18129/B9.bioc.HiCParser <https://doi.org/10.18129/B9.bioc.HiCParser>,
  https://github.com/emaigne/HiCParser/HiCParser - R package version 0.1.0,
  <http://www.bioconductor.org/packages/HiCParser>.

As a BibTeX entry :

    @Manual{hicparser,
      title = {A multiple format Hi-C data parser},
      author = {Elise Maigné and Matthias Zytnicki},
      year = {2024},
      url = {http://www.bioconductor.org/packages/HiCParser},
      note = {https://github.com/emaigne/HiCParser/HiCParser - R package version 0.1.0},
      doi = {10.18129/B9.bioc.HiCParser},
    }

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
