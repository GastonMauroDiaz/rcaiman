.is_even <- function(x) round(x/2) == x/2
.is_whole <- function(x) round(x) == x

.was_normalized <- function(r, name = "r") {
  if (.get_max(r) > 1)
    warning(paste("Please check if \"", name, "\" was correctly normalized."))
}

.is_single_layer_raster <- function(r, name = "r") {
  msn <- paste0("The argument \"", name,
               "\" should be a single-layer object",
               "from the SpatRaster class.")
  if (!is(r, "SpatRaster")) {
      stop(msn)
    } else {
      if (terra::nlyr(r) != 1)  stop(msn)
    }
}

.check_if_r_z_and_a_are_ok <- function(r, z, a) {
  .is_single_layer_raster(r, "r")
  # .was_normalized(r)
  .is_single_layer_raster(z, "z")
  .is_single_layer_raster(a, "a")
  stopifnot(ncol(z) == nrow(z))
  stopifnot(.get_max(z) <= 90)
  stopifnot(.get_max(a) <= 360)
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

.this_requires_EBImage <- function() {
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop(paste("Package \"EBImage\" needed for this function to work.",
               "Please install it. Find instructions here:",
               "https://bioconductor.org/packages/release/bioc/html/EBImage.html"
    ),
    call. = FALSE)
  }
}

.is_logic_and_NA_free <- function(x, name = "x") {
  stopifnot(!any(is.na(x[])))
  stopifnot(is.logical(x[1][,]))
}

.this_requires_conicfit <- function() {
  if (!requireNamespace("conicfit", quietly = TRUE)) {
    stop(paste("Package \"conicfit\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }
}
