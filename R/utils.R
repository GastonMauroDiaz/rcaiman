.radian2degree <- function(x) x * 180 / pi
.degree2radian <- function(x) x * pi / 180

.get_max <- function(r) max(r[], na.rm = TRUE)
.get_min <- function(r) min(r[], na.rm = TRUE)

.is_even <- function(x) round(x/2) == x/2
.is_whole <- function(x) round(x) == x
.is_integerish <- function(x) x == round(x)

.decode_label <- function(label) {
  sector_ID <- trunc(label / 1000)
  rings_ID <- label - sector_ID * 1000
  data.frame(sector_ID, rings_ID)
}

.calc_rmse <- function(x) sqrt(mean(x^2))

.check_if_r_was_normalized <- function(r, name = "r") {
  if (max(r[], na.rm = TRUE) > 1)
    warning(paste("Please check if \"", name, "\" was correctly normalized."))
}

.check_if_r_z_and_a_are_ok <- function(r, z, a) {
  stopifnot(class(r) == "SpatRaster")
  .check_if_r_was_normalized(r)
  stopifnot(class(z) == "SpatRaster")
  stopifnot(class(a) == "SpatRaster")
  stopifnot(.get_max(z) <= 90)
  terra::compareGeom(r, z)
  terra::compareGeom(z, a)
}

.is_class_from_colorspace <- function(x) {
  error_msn <- paste("\"target_color\" must be a subclass of the",
                     "virtual class named \"color\" from",
                     "\"colorspace\" package")
  if (!isS4(x)) stop(error_msn)
  if (attr(class(x), "package") != "colorspace") stop(error_msn)
}

.is_logic_and_NA_free <- function(x) {
  stopifnot(!any(is.na(x[])))
  stopifnot(is.logical(x[1][,]))
}
