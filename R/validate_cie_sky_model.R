#' Validate CIE sky models
#'
#' Validate CIE sky models fitted with [fit_cie_sky_model()] and
#' [ootb_fit_cie_sky_model()].
#'
#' The validation of the CIE sky model is done with a k-fold approach with k
#' equal to 10, following \insertCite{Kohavi1995;textual}{rcaiman}, and for the
#' values of relative luminance. The regression of the predicted vs. observed
#' values was done following \insertCite{Pineiro2008;textual}{rcaiman}. The
#' is_outlier were determined following \insertCite{Leys2013;textual}{rcaiman} and
#' with threshold equal to 3.
#'
#' @inheritParams fit_cie_sky_model
#' @param model An object of the class _list_. The output of
#'   [fit_cie_sky_model()].
#' @param k Numeric vector of length ones. Number of folds.
#'
#' @returns A _list_ with the following components:
#' \itemize{
#'   \item An object of class `lm` (see [stats::lm()]).
#'   \item predicted values.
#'   \item observed vales.
#'   \item the coefficient of determination (\eqn{r^2}).
#'   \item The root mean squared error (RMSE).
#'   \item The median absolute error (MAE).
#'   \item A logical vector indicating is_outlier within the sky points set attached to the 'model' argument.
#' }
#' @export
#' @family Tool Functions Functions
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' # See fit_cie_sky_model() for details on below file
#' path <- system.file("external/sky_points.gpkg",
#'                     package = "rcaiman")
#' sky_points <- terra::vect(path)
#' sky_points <- terra::extract(caim, sky_points, cells = TRUE)
#' sky_points <- terra::rowColFromCell(caim, sky_points$cell) %>% as.data.frame()
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(caim$Blue)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' xy <- c(210, 451) #taken with click() after x11(), then hardcoded here
#' sun_coord <- zenith_azimuth_from_row_col(z, a, c(nrow(z) - xy[2],xy[1]))
#' points(sun_coord$row_col[2], nrow(caim) - sun_coord$row_col[1],
#'        col = 3, pch = 1)
#'
#' rr <- extract_rel_radiance(caim$Blue, z, a, sky_points)
#'
#' set.seed(7)
#' model <- fit_cie_sky_model(rr, sun_coord,
#'                            general_sky_type = "Clear",
#'                            twilight = FALSE,
#'                            method = "CG")
#' model_validation <- validate_cie_sky_model(model, rr)
#' model_validation$r_squared
#' model_validation$rmse
#' }
validate_cie_sky_model <- function(model, rr, k = 10) {

  stopifnot(length(k) == 1)
  stopifnot(.is_whole(k))
  stopifnot(k >= 3)

  # START K-fold method ####
  ## k=10 based on https://dl.acm.org/doi/10.5555/1643031.1643047
  folds <- seq_along(rr$sky_points$row)
  folds <- split(folds, 1:k) %>% suppressWarnings()
  x <- c()
  y <- c()
  for (i in 1:k) {
    rr.2 <- rr
    rr.2$sky_points <- rr$sky_points[-folds[[i]],]
    model.2 <- fit_cie_sky_model(rr.2,  model$sun_coord,
                                 custom_sky_coef = model$coef + .noise(0.1),
                                 twilight = 90,
                                 method = model$method)

    x <- c(x, .cie_sky_model(AzP = rr$sky_points[folds[[i]], "a"] %>%
                               .degree2radian(),
                             Zp = rr$sky_points[folds[[i]], "z"] %>%
                               .degree2radian(),
                             AzS = model$sun_coord$zenith_azimuth[2] %>%
                               .degree2radian(),
                             Zs =  model$sun_coord$zenith_azimuth[1] %>%
                               .degree2radian(),
                             model$coef[1], model$coef[2], model$coef[3],
                             model$coef[4], model$coef[5]))
    y <- c(y, rr$sky_points[folds[[i]], "rr"])
  }

  #following Leys2013 10.1016/j.jesp.2013.03.013
  error <- y - x
  u <- abs((error - stats::median(error)) / stats::mad(error)) < 3
  x <- x[u]
  y <- y[u]

  reg <- lm(x~y) #following Pineiro2008 10.1016/j.ecolmodel.2008.05.006
  # END K-fold method ####

  list(lm = reg,
       pred = reg$model$x,
       obs = reg$model$y,
       r_squared = summary(reg) %>% .$r.squared,
       rmse = .calc_rmse(reg$model$y - reg$model$x),
       mae = stats::median(abs(reg$model$y - reg$model$x)),
       is_outlier = !u)
}
