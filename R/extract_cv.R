#' Extract local CV from sky points
#'
#' @description
#' Compute a robust coefficient of variation (CV) at each sky point using a
#' \eqn{3 \times 3} neighbourhood in a raster image.
#'
#' @param sky_points `data.frame` with columns `row` and `col` (raster
#'   coordinates).
#'
#' @inheritParams extract_dn
#'
#' @details
#' For each sky point, the function extracts the radiance values in a
#' \eqn{3 \times 3} window centered on the corresponding pixel of `r`, and
#' computes a robust CV defined as:
#'
#' \deqn{\mathrm{CV}_\mathrm{robust} = \mathrm{MAD} / \mathrm{median}}
#'
#' where MAD is the median absolute deviation.
#' This measure is suited for detecting local fluctuations caused by bright
#' canopy elements or mixed pixels, while being less sensitive to outliers than
#' a mean–sd based CV.
#'
#' The extraction procedure mirrors the spatial logic used in
#' [extract_dn()], but instead of averaging the neighbourhood, it applies a
#' robust dispersion measure (MAD/median) to quantify local radiance
#' variability around each point.
#'
#' @return Numeric vector of length equal to `nrow(sky_points)`, containing one
#'   robust CV value per sky point.
#'
#' @seealso
#' [extract_dn()] for extracting digital numbers and [tune_sky_sampling()] for
#' an application of this metric.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' # See fit_cie_model() for details on below file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(caim$Blue)
#' points(sky_points$col, nrow(caim) - sky_points$row, col = 2, pch = 10)
#'
#' extract_cv(caim$Blue, sky_points)
#'
#' }
extract_cv <- function(r, sky_points) {
  .assert_single_layer(r)
  .check_sky_points(sky_points)

  cells <- terra::cellFromRowCol(r, sky_points$row, sky_points$col)
  xy <-  terra::xyFromCell(r, cells)

  .fn <- function(x) stats::mad(x) / stats::median(x)
  Map(
    function(x, y) {
      ma <- expand.grid(c(-1,0,1) + x, c(-1,0,1) + y)
      .fn(terra::extract(r, ma, method = "simple", na.rm = TRUE)[,-1])
    },
    xy[,1],
    xy[,2]) %>% unlist()
}
