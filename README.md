
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rcaiman

<!-- badges: start -->
<!-- badges: end -->

Its main strength is to classify hemispherical photographs of the plant
canopy with algorithms specially developed for such a task and well
documented in [IEEE Geoscience and Remote Sensing
Letters](https://ieeexplore.ieee.org/document/7103294?arnumber=7103294),
[Canadian Journal of Forest
Research](https://cdnsciencepub.com/doi/full/10.1139/cjfr-2018-0006),
and [Methods in Ecology and
Evolution](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14059).
It supports non-circular hemispherical photography such as the ones
usually obtained with auxiliary fisheye lenses attached to smartphones.

It offers a simple method to calibrate fisheye lense, more details
[here](https://www.sciencedirect.com/science/article/abs/pii/S0168192324001357?via%3Dihub).

For digital cover photography, please refer to the [coverR
package](https://link.springer.com/article/10.1007/s00468-022-02338-5),
although you might find some functions from *rcaiman* useful, such as
*enhance_caim()*.

It integrates well with other software, such as 
[CIMES](http://jmnw.free.fr/) and
[Hemisfer](https://www.schleppi.ch/patrick/hemisfer/).

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
