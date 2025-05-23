% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_caim.R
\name{write_caim}
\alias{write_caim}
\title{Write canopy image}
\usage{
write_caim(caim, path, bit_depth)
}
\arguments{
\item{caim}{\linkS4class{SpatRaster}.}

\item{path}{Character vector of length one. Path for writing the image.}

\item{bit_depth}{Numeric vector of length one.}
}
\value{
No return value. Called for side effects.
}
\description{
Wrapper function for \code{\link[terra:writeRaster]{terra::writeRaster()}}.
}
\examples{
\dontrun{
caim <- read_caim() \%>\% normalize_minmax(., 0, 255)
write_caim(caim * 2^8-2, file.path(tempdir(), "test_8bit"), 8)
write_caim(caim * 2^16-2, file.path(tempdir(), "test_16bit"), 16)
# Note: the normalized values are scaled by multiplying by 2^bit_depth-2
# to avoid storing in the maximum bin because those values will be
# interpreted as NA by read_caim(), and that is undesired.
}
}
