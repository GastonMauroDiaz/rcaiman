% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_bin.R
\name{write_bin}
\alias{write_bin}
\title{Write binarized images}
\usage{
write_bin(bin, path)
}
\arguments{
\item{bin}{\linkS4class{SpatRaster}.}

\item{path}{Character vector of length one. Path for writing the image.}
}
\value{
No return value. Called for side effects.
}
\description{
Wrapper functions for \code{\link[terra:writeRaster]{terra::writeRaster()}}.
}
\examples{
\dontrun{
z <- zenith_image(1000, lens())
m <- !is.na(z)
my_file <- file.path(tempdir(), "mask")
write_bin(m, my_file)
my_file <- as.filename(my_file) \%>\%
  insert(., ext = "tif", replace = TRUE) \%>\%
  as.character()
m_from_disk <- read_bin(my_file)
plot(m - m_from_disk)
}
}
