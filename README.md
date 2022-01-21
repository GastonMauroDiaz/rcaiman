
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rcaiman

<!-- badges: start -->
<!-- badges: end -->

Its main strength is to classify hemispherical photographs of the plant
canopy with algorithms specially developed for such a task and well
documented in [IEEE Geoscience and Remote Sensing
Letters](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=7103294)
and [Canadian Journal of Forest
Research](http://www.nrcresearchpress.com/doi/full/10.1139/cjfr-2018-0006).

It also supports non-circular hemispherical photography, as the ones
usually obtained with auxiliary fisheye lenses attached to smartphones.
For digital cover photography, please refer to the [coverR
package](https://www.biorxiv.org/content/10.1101/2022.01.13.475850v1),
although you might find some functions from *rcaiman* useful, such as
*enhance_caim()*.

It integrates well with other software, such as [HSP
software](https://www.to.ee/eng/services/research_services/software),
[CIMES](http://jmnw.free.fr/), and
[Hemisfer](https://www.schleppi.ch/patrick/hemisfer/). The former is
ideal for carrying on the preprocessing of photographs acquired in
diffuse light conditions, by that I mean: extracting authentic blue
pixels, reproject to equidistant, and vignetting correction. After the
processing, you may obtain good results with *ootb_mblt()*. **Hemisfer**
and **CIMES** can be used to calculate canopy metrics, such as the leaf
area index, from the binarized images.

## Installation

You can install the development version of rcaiman from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GastonMauroDiaz/rcaiman")
```

You can install the release version from CRAN with:

``` r
install.packages("rcaiman")
```
