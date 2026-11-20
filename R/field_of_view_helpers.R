.internal_coverage_angles <- function(g, name) {
  .assert_spatraster(g, name = name)
  if (!identical(names(g), c("off_axis", "internal_azimuth"))) {
    stop("`", name, "` must have layers off_axis and internal_azimuth, ",
         "in that order.", call. = FALSE)
  }
  if (!terra::hasValues(g))
    stop("`", name, "` must contain direction values.", call. = FALSE)
  angles <- terra::values(g)
  .check_vector(angles[, 1], "numeric", allow_na = TRUE, name = "off_axis")
  .check_vector(angles[, 2], "numeric", allow_na = TRUE, name = "internal_azimuth")
  if (any(angles[, 1] < 0 | angles[, 1] > 180, na.rm = TRUE) ||
      any(angles[, 2] < 0 | angles[, 2] >= 360, na.rm = TRUE)) {
    stop("Internal angles must be in [0, 180] and [0, 360) degrees.", call. = FALSE)
  }
  angles
}

.check_internal_directions <- function(directions, name, allow_na = FALSE) {
  .check_columns(directions, required_columns = c("off_axis", "internal_azimuth"))
  .check_vector(directions$off_axis, "numeric", allow_na = allow_na,
                name = paste0(name, "$off_axis"))
  .check_vector(directions$internal_azimuth, "numeric", allow_na = allow_na,
                name = paste0(name, "$internal_azimuth"))
  if (any(directions$off_axis < 0 | directions$off_axis > pi, na.rm = TRUE) ||
      any(directions$internal_azimuth < 0 | directions$internal_azimuth >= 2 * pi,
          na.rm = TRUE)) {
    stop("Internal angles must be in [0, pi] and [0, 2 * pi) radians.", call. = FALSE)
  }
  invisible(TRUE)
}

.empty_field_of_view <- function() {
  data.frame(off_axis = numeric(), internal_azimuth = numeric())
}

.spherical_field_of_view <- function(field_of_view, name) {
  .check_internal_directions(field_of_view, name)
  if (!nrow(field_of_view)) return(s2::s2_geog_from_text("POLYGON EMPTY"))
  if (nrow(unique(field_of_view[c("off_axis", "internal_azimuth")])) < 3L) {
    stop("`", name, "` must delimit a nonzero area with at least three vertices.",
         call. = FALSE)
  }
  # Longitude increases with clock azimuth. This reflection maps the clockwise
  # image contour to S2's convention: the interior is on the left of each edge.
  polygon <- tryCatch(s2::s2_make_polygon(
    field_of_view$internal_azimuth * 180 / pi,
    90 - field_of_view$off_axis * 180 / pi, oriented = TRUE),
    error = function(e) stop("Invalid contour in `", name, "`: ",
                             conditionMessage(e), call. = FALSE))
  if (s2::s2_area(polygon, radius = 1) == 0) {
    stop("`", name, "` must delimit a nonzero spherical area.", call. = FALSE)
  }
  polygon
}

.field_of_view_from_sphere <- function(polygon) {
  if (s2::s2_is_empty(polygon)) return(.empty_field_of_view())
  # terra is used only to read the vertices, without planar geometric operations.
  vertices <- terra::geom(terra::vect(s2::s2_as_text(polygon, precision = 17)),
                           df = TRUE)
  if (length(unique(vertices$part)) != 1L || any(vertices$hole != 0)) {
    stop("The intersection requires more than one contour; ",
         "only a single contour without holes is supported.", call. = FALSE)
  }
  off_axis <- pmax(0, pmin(pi, (90 - vertices$y) * pi / 180))
  internal_azimuth <- (vertices$x * pi / 180) %% (2 * pi)
  internal_azimuth[internal_azimuth >= 2 * pi] <- 0
  data.frame(off_axis = off_axis, internal_azimuth = internal_azimuth)
}

.directions_in_field_of_view <- function(directions, polygon) {
  selected <- rep(FALSE, nrow(directions))
  valid <- !is.na(directions$off_axis) & !is.na(directions$internal_azimuth)
  if (!any(valid) || s2::s2_is_empty(polygon)) return(selected)
  points <- s2::s2_geog_point(directions$internal_azimuth[valid] * 180 / pi,
                             90 - directions$off_axis[valid] * 180 / pi)
  # A contour-to-angle round trip must not remove its own boundary vertices.
  selected[valid] <- s2::s2_dwithin(polygon, points,
    distance = 64 * .Machine$double.eps, radius = 1)
  selected
}
