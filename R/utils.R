.radian2degree <- function(x) x * 180 / pi
.degree2radian <- function(x) x * pi / 180

.get_max <- function(r) max(r[], na.rm = TRUE)
.get_min <- function(r) min(r[], na.rm = TRUE)

.is_even <- function(x) round(x/2) == x/2

.is_whole <- function(x) round(x) == x

.check_if_r_was_normalized <- function(r, name = "r") {
  if (max(r[], na.rm = TRUE) > 1)
    warning(paste("Please check if \"", name, "\" was correctly normalized."))
}

.check_if_r_z_and_a_are_ok <- function(r, z, a) {
  stopifnot(class(r) == "RasterLayer")
  .check_if_r_was_normalized(r)
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(a) == "RasterLayer")
  stopifnot(.get_max(z) <= 90)
  compareRaster(r, z)
  compareRaster(z, a)
}

.is_class_from_colorspace <- function(x) {
  error_msn <- paste("\"target_color\" must be a subclass of the",
                     "virtual class named \"color\" from",
                     "\"colorspace\" package")
  if (!isS4(x)) stop(error_msn)
  if (attr(class(x), "package") != "colorspace") stop(error_msn)
}

