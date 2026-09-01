#' Defuzzify a fuzzy classification
#'
#' Converts fuzzy membership values into a binary classification using a
#' regional approach that preserves aggregation consistency between the fuzzy
#' and binary representations.
#'
#' The conversion is applied within segments defined by `seg`,
#' ensuring that, in each segment, the aggregated Boolean result matches the
#' aggregated fuzzy value. This approach is well suited for converting subpixel
#' estimates, such as gap fraction, into binary outputs.
#'
#' @param mem numeric [terra::SpatRaster-class] of one layer. Degree of
#'   membership in a fuzzy classification.
#'
#' @inheritParams compute_mean
#'
#' @return Logical [terra::SpatRaster-class] of the same dimensions as `mem`,
#'   where each pixel value represents the binary version of `mem` after
#'   applying the regional defuzzification procedure.
#'
#' @note This method is also available in the HSP software package. See
#'   [hsp_compat()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' plot(r)
#'
#' # Manual method following Lang et al. (2010)
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(r)
#' points(sky_points$col, nrow(r) - sky_points$r, col = 2, pch = 10)
#'
#' # x11()
#' # plot(r)
#' # sun_angles <- click(c(z, a), 1) %>% as.numeric()
#' sun_angles <- c(zen_deg = 49.5, azi_deg = 27.4) #taken with above lines then hardcoded
#'
#' sun_row_col <- row_col_from_zen_azi(z, a,
#'                                     sun_angles["zen_deg"] * pi/180,
#'                                     sun_angles["azi_deg"] * pi/180)
#' points(sun_row_col[2], nrow(r) - sun_row_col[1],
#'        col = "yellow", pch = 8, cex = 3)
#'
#' sky_points <- extract_dn(r, sky_points)
#' sky_points <- add_directions(sky_points, z, a)
#' rr <- calc_relative_radiance(sky_points, 3)
#'
#' set.seed(7)
#' model <- fit_cie_model(rr, sun_angles,
#'                        general_sky_type = "Clear")
#'
#' plot(model$rr$sky_points$rr, model$rr$sky_points$pred)
#' abline(0,1)
#' lm(model$rr$sky_points$pred~model$rr$sky_points$rr) %>% summary()
#'
#' sky <- cie_image(z, a, model$sun_angles, model$coef) * model$rr$zenith_dn
#' names(sky) <- "DN"
#'
#' ratio <- r/sky
#' seg <- equal_area_segmentation(z, a, 50)
#' bin2 <- defuzzify(ratio, seg)
#' plot(bin2) # unsatisfactory results due to light conditions
#' }
defuzzify <- function (mem, seg) {
  .assert_single_layer(mem)
  if (any(.get_max(mem) > 1, .get_min(mem) < 0)) {
    warning("Check if 'mem' values are degree of membership in a fuzzy classification.")
  }
  .assert_single_layer(seg)
  .assert_same_geom(mem, seg)

  mem[is.na(mem)] <- 0

  .fun <- function(x) {
    no_of_pixels <- round(mean(x) * length(x))
    if (no_of_pixels > 0 &
        no_of_pixels != length(x)) {
      indices <- order(x, decreasing = TRUE)[1:no_of_pixels]
      x[indices] <- 1
      x[x != 1] <- 0
    } else {
      x <- as.numeric(x > 0.5)
    }
    x
  }

  cells <- mem
  terra::values(cells) <- 1:ncell(mem)
  cells <- tapply(terra::values(cells),
                   terra::values(seg), function(x) x)
  cells <- unlist(cells)
  bin <- tapply(terra::values(mem), terra::values(seg), .fun)
  mem[cells] <- unlist(bin)
  as.logical(mem)
}
