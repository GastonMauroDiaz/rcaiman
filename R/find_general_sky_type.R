#' Find the general sky type attributable to a given set of coefficients
#'
#' @param model An object of the class _list_. The result of calling
#'   [fit_cie_sky_model()].
#'
#' @returns A character vector of length one.
#' @export
#' @family Tool Functions
#'
#' @examples
#' \dontrun{
#' caim <- read_caim() %>% normalize()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' # Manual method following Lang et al. (2010) using QGIS
#' path <- system.file("external/sky_points.gpkg",
#'                     package = "rcaiman")
#' sky_points <- terra::vect(path)
#' sky_points <- terra::extract(caim, sky_points, cells = TRUE)
#' sky_points <- terra::rowColFromCell(caim, sky_points$cell) %>% as.data.frame()
#' colnames(sky_points) <- c("row", "col")
#' xy <- c(210, 451) #taken with click() after x11(), then hardcoded here
#' sun_coord <- zenith_azimuth_from_row_col(z, a, c(nrow(z) - xy[2],xy[1]))
#'
#' rr <- extract_rel_radiance(caim$Blue, z, a, sky_points)
#'
#' set.seed(7)
#' model <- fit_cie_sky_model(rr, sun_coord,
#'                            general_sky_type = "Clear",
#'                            twilight = FALSE,
#'                            method = "CG")
#' find_general_sky_type(model)
#' }
find_general_sky_type <- function(model) {
  path <- system.file("external", package = "rcaiman")
  skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))

  .find_std_sky_no <- function(model) {
    .median <- apply(skies[,1:5], 2, stats::median)
    .mad <- apply(skies[,1:5], 2, stats::mad)
    coef <- model$coef
    delta <- skies[, 1:5]
    for (i in 1:nrow(skies)) {
      delta[i,] <- delta[i,] - coef
      delta[i,] <- (delta[i,] - .median) / .mad #z-score
    }
    apply(delta^2, 1, sum) %>% which.min()
  }

  std_sky_no <- .find_std_sky_no(model)
  skies[std_sky_no, "general_sky_type"]
}
