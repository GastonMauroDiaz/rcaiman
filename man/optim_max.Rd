% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optim_max.R
\name{optim_max}
\alias{optim_max}
\title{Optimize a parameter of the function \code{\link[=normalize_minmax]{normalize_minmax()}}}
\usage{
optim_max(caim, bin, method = "BFGS")
}
\arguments{
\item{caim}{\linkS4class{SpatRaster}. The output of \code{\link[=read_caim]{read_caim()}}.}

\item{bin}{\linkS4class{SpatRaster}. This should be a preliminary binarization of
\code{r} useful for masking pixels that are very likely pure sky pixels.}

\item{method}{The method to be used. See \sQuote{Details}.  Can be abbreviated.}
}
\value{
Numeric vector of length one. The values for using as \code{mx} argument
with \code{\link[=normalize_minmax]{normalize_minmax()}}.
}
\description{
Wrapper function for \code{\link[stats:optim]{stats::optim()}}. Optimize the \code{mx} argument of the
function \code{\link[=normalize_minmax]{normalize_minmax()}} by maximizing colorfulness and minimizing
saturation.
}
\examples{
\dontrun{
caim <- read_caim()
z <- zenith_image(ncol(caim), lens())
a <- azimuth_image(z)
m <- !is.na(z)

mn <- quantile(caim$Blue[m], 0.01)
mx <- quantile(caim$Blue[m], 0.99)
r <- normalize_minmax(caim$Blue, mn, mx, TRUE)

bin <- apply_thr(caim$Blue, thr_isodata(caim$Blue[m]))

mx <- optim_max(caim, bin)
ncaim <- normalize_minmax(caim, mx = mx, force_range = TRUE)
plotRGB(ncaim*255)
plotRGB(normalize_minmax(caim)*255)
percentage_of_clipped_highlights(ncaim$Blue, m)
}
}
