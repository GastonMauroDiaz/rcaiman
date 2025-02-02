% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ootb_sky_reconstruction.R
\name{ootb_sky_reconstruction}
\alias{ootb_sky_reconstruction}
\title{Out-of-the-box sky reconstruction}
\usage{
ootb_sky_reconstruction(
  r,
  z,
  a,
  m,
  bin,
  g = NULL,
  sor_filter_cv = FALSE,
  sor_filter_dn = FALSE,
  refine_sun_coord = FALSE,
  interpolate = TRUE,
  input_sky_points = NULL
)
}
\arguments{
\item{r}{\linkS4class{SpatRaster}. A normalized greyscale image. Typically, the
blue channel extracted from a canopy photograph. Please see \code{\link[=read_caim]{read_caim()}}
and \code{\link[=normalize]{normalize()}}.}

\item{z}{\linkS4class{SpatRaster} built with \code{\link[=zenith_image]{zenith_image()}}.}

\item{a}{\linkS4class{SpatRaster} built with \code{\link[=azimuth_image]{azimuth_image()}}.}

\item{m}{\linkS4class{SpatRaster}. A mask, check \code{\link[=mask_hs]{mask_hs()}}.}

\item{bin}{\linkS4class{SpatRaster}. This should be a preliminary binarization of
\code{r} useful for masking pixels that are very likely pure sky pixels.}

\item{g}{\linkS4class{SpatRaster} built with \code{\link[=sky_grid_segmentation]{sky_grid_segmentation()}} or
\code{\link[=chessboard]{chessboard()}}.}

\item{sor_filter_cv}{Logical vector of length one}

\item{sor_filter_dn}{Logical vector of length one}

\item{refine_sun_coord}{Logical vector of length one}

\item{interpolate}{Logical vector of length one. If \code{TRUE},
\code{\link[=interpolate_sky_points]{interpolate_sky_points()}} will be used.}

\item{input_sky_points}{An object of class \emph{data.frame} with the same
structure than the result of a call to \code{\link[=extract_sky_points]{extract_sky_points()}}. The
\href{https://imagej.net/ij/}{ImageJ} software package can be used to manually
digitize points. See \code{\link[=extract_dn]{extract_dn()}} for details.}
}
\value{
An object from the class \emph{list} that includes the following: (1) the
reconstructed sky (\linkS4class{SpatRaster}), (2) the output produced by
\code{\link[=fit_cie_sky_model]{fit_cie_sky_model()}}, (3) an object from the class \emph{list} that includes an
object from the class \code{lm} (see \code{\link[stats:lm]{stats::lm()}}) and the RMSE, both being the
result of validating (2) with a k-fold approach and following
\insertCite{Pineiro2008;textual}{rcaiman},(4) the \code{dist_to_plant} argument
used when \code{\link[=fit_cie_sky_model]{fit_cie_sky_model()}} was called, (5) the \code{sky_points} argument
used when \code{\link[=extract_rl]{extract_rl()}} was called,  and (6) the out-of-range index (see
details).
}
\description{
Build an above canopy image from a single below canopy image
}
\details{
This function is a hard-coded version of a pipeline that uses these main
functions \code{\link[=fit_cie_sky_model]{fit_cie_sky_model()}} and \code{\link[=interpolate_sky_points]{interpolate_sky_points()}}.

The pipeline is an automatic version of the
\insertCite{Lang2010;textual}{rcaiman} method. A paper for thoroughly
presenting and testing this pipeline is under preparation.

The \strong{out-of-range index} is calculated as foollow:

\eqn{\sum_{i = 1}^{N}(r_i/sky_i)^2},

where \eqn{r} is the \code{r} argument, \eqn{sky} is the
raster obtained from the fitted model with \code{\link[=cie_sky_model_raster]{cie_sky_model_raster()}} and
\code{zenith_dn}, \eqn{i} is the index that represents the position of a given
pixel on the raster grid, and \eqn{N} is the total number of pixels that
satisfy: \eqn{r_i/sky_i<0} or \eqn{r_i/sky_i>1}.
}
\examples{
\dontrun{
caim <- read_caim()
z <- zenith_image(ncol(caim), lens())
a <- azimuth_image(z)
m <- !is.na(z)

r <- caim$Blue

bin <- regional_thresholding(r, rings_segmentation(z, 30),
                             method = "thr_isodata")
g <- sky_grid_segmentation(z, a, 10)
sun_coord <- extract_sun_coord(r, z, a, bin, g)
sun_coord$zenith_azimuth

.a <- azimuth_image(z, orientation = sun_coord$zenith_azimuth[2]+90)
seg <- sectors_segmentation(.a, 180) * rings_segmentation(z, 30)
bin <- regional_thresholding(r, seg, method = "thr_isodata")
plot(bin)


mx <- optim_normalize(caim, bin)
caim <- normalize(caim, mx = mx, force_range = TRUE)
ecaim <- enhance_caim(caim, m, HSV(239, 0.85, 0.5))
bin <- apply_thr(ecaim, thr_isodata(ecaim[m]))

set.seed(7)
sky <- ootb_sky_reconstruction(r, z, a, m, bin,
                               sor_filter_cv = TRUE, sor_filter_dn = TRUE)

sky$sky
plot(sky$sky)
sky$model_validation$rmse
plot(r/sky$sky>1.15)
plot(sky$model_validation$reg$model$x, sky$model_validation$reg$model$y)
abline(0,1)
error <- sky$model_validation$reg$model$x - sky$model_validation$reg$model$y
plot(sky$sky_points$z, error,
                     xlab = "zenith angle", ylab = "relative radiance error")

plot(bin)
points(sky$sky_points$col, nrow(caim) - sky$sky_points$row, col = 2, pch = 10)

}
}
\references{
\insertAllCited{}
}
\seealso{
Other Sky Reconstruction Functions: 
\code{\link{cie_sky_model_raster}()},
\code{\link{fit_cie_sky_model}()},
\code{\link{fit_coneshaped_model}()},
\code{\link{fit_trend_surface}()},
\code{\link{fix_reconstructed_sky}()},
\code{\link{interpolate_sky_points}()}
}
\concept{Sky Reconstruction Functions}
