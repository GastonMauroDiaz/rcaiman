% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fisheye_to_pano.R
\name{fisheye_to_pano}
\alias{fisheye_to_pano}
\title{Fisheye to panoramic}
\usage{
fisheye_to_pano(r, z, a, fun = mean, angle_width = 1)
}
\arguments{
\item{r}{\linkS4class{SpatRaster}. A fish-eye image.}

\item{z}{\linkS4class{SpatRaster} built with \code{\link[=zenith_image]{zenith_image()}}.}

\item{a}{\linkS4class{SpatRaster} built with \code{\link[=azimuth_image]{azimuth_image()}}.}

\item{fun}{A \code{function} that takes a vector as input and returns a
one-length numeric or logical vector as output (e.g. mean).}

\item{angle_width}{Numeric vector of length one. It should be \verb{30, 15, 10, 7.5, 6, 5, 3.75, 3, 2.5, 1.875, 1} or \code{0.5} degrees. This constrain is
rooted in the requirement of a value able to divide both the \code{0} to \code{360}
and \code{0} to \code{90} ranges into a whole number of segments.}
}
\description{
Fisheye to panoramic (cylindrical projection)
}
\details{
An early version of this function was used in
\insertCite{Diaz2021;textual}{rcaiman}.
}
\examples{
\dontrun{
caim <- read_caim()
z <- zenith_image(ncol(caim), lens())
a <- azimuth_image(z)
pano <- fisheye_to_pano(caim, z, a)
plotRGB(pano \%>\% normalize_minmax() \%>\% multiply_by(255))
}
}
\references{
\insertAllCited{}
}
