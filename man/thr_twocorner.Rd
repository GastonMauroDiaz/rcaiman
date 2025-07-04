% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thr_twocorner.R
\name{thr_twocorner}
\alias{thr_twocorner}
\title{Calculate the thresholds of the two-corner method}
\usage{
thr_twocorner(
  x,
  diagnose = FALSE,
  window_length = NULL,
  slope_reduction = TRUE
)
}
\arguments{
\item{x}{Numeric vector or a single-column \emph{matrix} or \emph{data.frame}
able to be coerced to numeric.}

\item{diagnose}{Logical vector of length one. If \code{TRUE}, then a figure will
be send to the graphical device showing the geometrical contruction made to
find the threholds.}

\item{window_length}{Numeric vector. One or many odd numbers indicating the
size of the window to be used to smooth the histogram.}

\item{slope_reduction}{Logical vector of length one. If \code{TRUE}, then the
slope reduction method by Macfarlane (2011) will be applied.}
}
\value{
a list.
}
\description{
Calculate the thresholds of the two-corner method
}
\examples{
\dontrun{
}
}
