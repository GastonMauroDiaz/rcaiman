% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/masking.R
\name{masking}
\alias{masking}
\title{Image masking}
\usage{
masking(r, m, RGB = c(1, 0, 0))
}
\arguments{
\item{r}{\linkS4class{SpatRaster}. The image. Values should be normalized,
see \code{\link[=normalize_minmax]{normalize_minmax()}}. Only methods for images with one or three
layers have been implemented.}

\item{m}{\linkS4class{SpatRaster}. A mask. For hemispherical photographs,
check \code{\link[=select_sky_vault_region]{select_sky_vault_region()}}.}

\item{RGB}{Numeric vector of length three. RGB color code. Red is the default
color.}
}
\value{
An object of class \linkS4class{SpatRaster} that is
\code{r} but with areas where \code{m} is equal to zero painted in a solid
color. If \code{r} is a single layer image, then the layer is triplicated
to allow the use of color.
}
\description{
Image masking
}
\examples{
\dontrun{
 r <- read_caim()
 z <- zenith_image(ncol(r), lens())
 a <- azimuth_image(z)
 m <- select_sky_vault_region(z, 20, 70) & select_sky_vault_region(a, 90, 180)

 masked_caim <-  masking(normalize_minmax(r), m)
 plotRGB(masked_caim * 255)

 masked_bin <- masking(apply_thr(r$Blue, 125), m)
 plotRGB(masked_bin * 255)
 }

}
\seealso{
\code{\link[=select_sky_vault_region]{select_sky_vault_region()}}
}
