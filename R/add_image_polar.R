#' Add polar predictors to image coordinates
#'
#' @description
#' Enrich continuous image coordinates with relative radius and image-plane
#' angle for model fitting and prediction.
#'
#' @details
#' The default center is `c(width / 2, height / 2)` in the `coord`
#' representation described in [y_down_from_y()]. Supply `center_of_vision`
#' to use a declared optical center. Offsets from this center
#' retain the image axes: horizontal offsets increase rightward, and vertical
#' offsets increase downward, without converting to the ordinary `xy`
#' representation.
#'
#' Radial distances are divided by the image semidiagonal:
#' \deqn{r_{rel} = \frac{\sqrt{(x-c_x)^2 + (y_{down}-c_y)^2}}
#'   {\sqrt{(w/2)^2 + (h/2)^2}}}
#' where \eqn{(c_x, c_y)} is the chosen center. It has relative radius zero.
#' The image corners have relative radius one only with the default center.
#' This scale remains the image semidiagonal when the center is moved.
#'
#' The image angle is measured from the rightward horizontal axis and increases
#' clockwise because the vertical image axis points downward. It is an
#' image-plane predictor, not internal azimuth. At the exact chosen
#' center, the angle is zero following the `atan2(0, 0)` convention used by
#' R; this value does not identify a unique direction.
#'
#' @note
#' Dimensions must describe the same image frame as the input coordinates.
#' Points outside that frame are retained without clipping.
#'
#' @param sampling_coord `data.frame` with finite numeric columns `"x"` and
#'   `"y_down"`. Additional columns are preserved. The output column names
#'   `"relative_radius"` and `"image_angle"` must not already exist.
#' @param width,height integer-like numeric vectors of length one. Positive
#'   image width and height in pixels, respectively.
#' @param center_of_vision finite numeric vector of length two, `c(x, y_down)`.
#'   Center from which offsets are measured. If `NULL`, use the image center.
#'   For fitting an internal geometry, supply the same definitive integer
#'   center that will be used by [set_internal_geometry()].
#'
#' @return A `data.frame` with all input columns and the following added columns:
#' \describe{
#'   \item{`relative_radius`}{Distance to the chosen center divided by the
#'     image semidiagonal.}
#'   \item{`image_angle`}{Image-plane angle in radians, in `[-pi, pi]`.}
#' }
#' Row order and duplicates are preserved. A zero-row input with the required
#' columns are returned with additional columns.
#'
#' @seealso [internal_directions_from_coord()], [coord_from_indices()]
#'
#' @md
#' @export
#' @examples
#' sampling_coord <- data.frame(
#'   x = c(400, 700, 400),
#'   y_down = c(100, 300, 500),
#'   id = c("top", "right", "bottom")
#' )
#' add_image_polar(sampling_coord, width = 800, height = 600)
#' # Use the definitive center when preparing an internal geometry model.
#' add_image_polar(sampling_coord, width = 800, height = 600,
#'                 center_of_vision = c(395, 302))
#'
#' sampling_coord <- coord_from_indices(col = c(1, 800), row = c(1, 600))
#' add_image_polar(sampling_coord, width = 800, height = 600)
add_image_polar <- function(sampling_coord, width, height, center_of_vision = NULL) {
  .check_columns(
    sampling_coord,
    required_columns = c("x", "y_down"),
    reserved_columns = c("relative_radius", "image_angle")
  )
  .check_vector(width, "integerish", 1, sign = "positive")
  .check_vector(height, "integerish", 1, sign = "positive")
  if (is.null(center_of_vision)) center_of_vision <- c(width / 2, height / 2)
  .check_vector(center_of_vision, "numeric", 2)
  for (column in c("x", "y_down")) {
    .check_vector(
      sampling_coord[[column]],
      "numeric",
      name = paste0("sampling_coord$", column)
    )
  }
  if (any(!is.finite(c(sampling_coord$x, sampling_coord$y_down)))) {
    stop("`x` and `y_down` must contain finite values.", call. = FALSE)
  }

  polar <- .image_polar_predictors(
    sampling_coord$x - center_of_vision[1],
    sampling_coord$y_down - center_of_vision[2], width, height
  )
  sampling_coord$relative_radius <- polar$relative_radius
  sampling_coord$image_angle <- polar$image_angle
  sampling_coord
}

.image_polar_predictors <- function(dx, dy, width, height) {
  data.frame(relative_radius = sqrt(dx^2 + dy^2) /
               sqrt((width / 2)^2 + (height / 2)^2),
             image_angle = atan2(dy, dx))
}
