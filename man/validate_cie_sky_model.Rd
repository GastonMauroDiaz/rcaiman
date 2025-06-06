% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_cie_sky_model.R
\name{validate_cie_sky_model}
\alias{validate_cie_sky_model}
\title{Validate CIE sky models}
\usage{
validate_cie_sky_model(model, rr, k = 10)
}
\arguments{
\item{model}{An object of the class \emph{list}. The output of
\code{\link[=fit_cie_sky_model]{fit_cie_sky_model()}}.}

\item{rr}{An object of class \emph{list}. The output of \code{\link[=extract_rel_radiance]{extract_rel_radiance()}}
or an object with same structure and names.}

\item{k}{Numeric vector of length ones. Number of folds.}
}
\value{
A \emph{list} with the following components:
\itemize{
\item An object of class \code{lm} (see \code{\link[stats:lm]{stats::lm()}}).
\item predicted values.
\item observed vales.
\item the coefficient of determination (\eqn{r^2}).
\item The root mean squared error (RMSE).
\item The median absolute error (MAE).
\item A logical vector indicating is_outlier within the sky points set attached to the 'model' argument.
}
}
\description{
Validate CIE sky models fitted with \code{\link[=fit_cie_sky_model]{fit_cie_sky_model()}} and
\code{\link[=ootb_fit_cie_sky_model]{ootb_fit_cie_sky_model()}}.
}
\details{
The validation of the CIE sky model is done with a k-fold approach with k
equal to 10, following \insertCite{Kohavi1995;textual}{rcaiman}, and for the
values of relative radiance. The regression of the predicted vs. observed
values was done following \insertCite{Pineiro2008;textual}{rcaiman}. The
is_outlier were determined following \insertCite{Leys2013;textual}{rcaiman}
and with threshold equal to 3.
}
\examples{
\dontrun{
caim <- read_caim() \%>\% normalize_minmax()
z <- zenith_image(ncol(caim), lens())
a <- azimuth_image(z)

# See fit_cie_sky_model() for details on below file
path <- system.file("external/sky_points.csv",
                    package = "rcaiman")
sky_points <- read.csv(path)
sky_points <- sky_points[c("Y", "X")]
colnames(sky_points) <- c("row", "col")
head(sky_points)
plot(caim$Blue)
points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)

# x11()
# plot(caim$Blue)
# sun_zenith_azimuth <- click(c(z, a), 1) \%>\% as.numeric()
sun_zenith_azimuth <- c(49.5, 27.42481) #taken with above lines then hardcoded

rr <- extract_rel_radiance(caim$Blue, z, a, sky_points)

set.seed(7)
model <- fit_cie_sky_model(rr, sun_zenith_azimuth,
                           general_sky_type = "Clear",
                           twilight = FALSE,
                           method = "CG")
model_validation <- validate_cie_sky_model(model, rr)
model_validation$r_squared
model_validation$rmse
}
}
\references{
\insertAllCited{}
}
