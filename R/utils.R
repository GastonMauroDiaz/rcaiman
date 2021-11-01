radian2degree <- function(x) x * 180 / pi
degree2radian <- function(x) x * pi / 180

.get_max <- function(r) max(r[], na.rm = TRUE)
.get_min <- function(r) min(r[], na.rm = TRUE)

.is_even <- function(x) round(x/2) == x/2

.is_whole <- function(x) round(x) == x

.is_integerish <- function(x) x == round(x)

.check_if_r_was_normalized <- function(r) {
  if (max(r[], na.rm = TRUE) > 1)
    warning("Please check if \"r\" was correctly normalized")
}

.check_if_r_z_and_a_are_ok <- function(r, z, a) {
  stopifnot(class(r) == "RasterLayer")
  .check_if_r_was_normalized(r)
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(a) == "RasterLayer")
  stopifnot(.get_max(z) < 90)
  compareRaster(r, z)
  compareRaster(z, a)
}

.this_requires_EBImage <- function() {
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop(paste("Package \"EBImage\" needed for this function to work.",
               "Please install it. Find instructions here:",
               "https://bioconductor.org/packages/release/bioc/html/EBImage.html"
               ),
         call. = FALSE)
  }
}

.calc_angular_distance <- function(z1, a1, z2, a2) {
  #https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  acos(pmax(pmin(cos(z1) * cos(z2) + sin(z1) * sin(z2) * cos(abs(a2 - a1)), 1), -1))
}


