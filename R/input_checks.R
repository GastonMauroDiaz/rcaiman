
# .check_*() verify structure. Classes of the base pkgs. ------------------

.check_sun_angles <- function(sun_angles) {
  # Structural check:
  if (!is.numeric(sun_angles) || length(sun_angles) != 2 ||
      !identical(names(sun_angles), c("zen_deg", "azi_deg"))) {
    stop("`sun_angles` must be a named numeric vector of length ",
         "two with names 'zen_deg' and 'azi_deg', in that order.")
  }
  # range check
  if (sun_angles["zen_deg"] < 0) {
    stop("Sun disk zenith angle cannot be negative.", call. = FALSE)
  }
  if (sun_angles["zen_deg"] > 96) {
    warning("Sun disk zenith angle should not be greater than 96 deg.")
  }
  if (any(sun_angles["azi_deg"] < 0 | sun_angles["azi_deg"] >= 360)) {
    stop("Sun disk azimuth angle must lie within [0, 360) deg.", call. = FALSE)
  }

  invisible(TRUE)
}




# .assert_*(), for semantic objects -------------------------------------

.assert_rgb3 <- function(x, squared = FALSE, name = deparse(substitute(x))) {
  .assert_spatraster(x, squared = squared,  name = name)
  if (terra::nlyr(x) != 3)
    stop(sprintf("`%s` must have three layers named `Red`, `Green`, `Blue`.", name), call. = FALSE)
  nm <- names(x)
  if (!identical(nm, c("Red","Green","Blue")))
    stop(sprintf("`%s` layers must be named `Red`, `Green`, `Blue`.", name), call. = FALSE)
  invisible(TRUE)
}

.assert_azimuth_image <- function(a) {
  .assert_single_layer(a, squared = TRUE)
  orientation <- attr(a, "orientation")
  if (is.null(orientation) || !is.numeric(orientation)) {
    stop("`orientation` attribute is missing or not numeric.", call. = FALSE)
  }
  .assert_raster_range(a, 0, 360, open_upper = TRUE, label = "Azimuth values")
}

.assert_area_weight_image <- function(w) {
  .assert_single_layer(w, squared = TRUE)
  projection_coef <- attr(w, "projection_coef")
  if (is.null(projection_coef) || !is.numeric(projection_coef)) {
    stop("`projection_coef` attribute is missing or not numeric.", call. = FALSE)
  }
  invisible(TRUE)
}

.assert_sky_grid <- function(g, name = deparse(substitute(g))) {
  .assert_single_layer(g, squared = TRUE, name)

  angle_width <- attr(g, "angle_width")
  ring_flag   <- attr(g, "first_ring_different")

  if (!is.numeric(angle_width) || length(angle_width) != 1 ||
      !is.finite(angle_width)) return(FALSE)
  if (!is.logical(ring_flag) || length(ring_flag) != 1 ||
      anyNA(ring_flag)) return(FALSE)

  nm <- names(g)
  if (is.null(nm)) return(FALSE)

  if (!all(nm == paste0("Sky grid, ", angle_width, " degrees"))) {
    stop(sprintf("`%s` must be the output of `sky_grid_segmentation()`.", name),
         call. = FALSE)
  }
  invisible(TRUE)
}

.assert_equal_area_seg <- function(equal_area_seg) {
  .assert_single_layer(equal_area_seg, squared = TRUE)
  if (names(equal_area_seg) != "Sky segments of equal area") {
    stop("`equal_area_seg` must be the output of `equal_area_segmentation()")
  }
}

.assert_cie_model <- function(model) {
  if (!is.list(model)) {
    stop("`model` must be a list returned by `fit_cie_model()`. ")
  }
  required_names <- c("rr", "opt_result", "coef", "sun_angles",
                      "method", "start" , "metric")
  if (!all(required_names %in% names(model))) {
    stop(sprintf("`model` must contain %s.",
                 paste(sprintf('"%s"', required_names), collapse = ", ")))
  }
}

.assert_neighbors <- function(x) {
  if (!is.list(x)) {
    stop("`neighbors` must be a list.")
  }

  ok <- vapply(x, function(e) {
    if (!is.list(e)) {
      return(FALSE)
    }

    if (!all(c("index", "distance") %in% names(e))) {
      return(FALSE)
    }

    if (!is.integer(e$index)) {
      return(FALSE)
    }

    if (!is.numeric(e$distance)) {
      return(FALSE)
    }

    if (length(e$index) != length(e$distance)) {
      return(FALSE)
    }

    TRUE

  }, logical(1))

  if (!all(ok)) {
    stop(
      "`neighbors` must be a list of lists with ",
      "`index` and `distance` components of equal length."
    )
  }

  invisible(TRUE)
}
