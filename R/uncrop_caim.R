#' @rdname crop_caim
#' @export
uncrop_caim <- function(r, upper_left_indices, original_dim) {
  .assert_spatraster(r)
  .check_vector(upper_left_indices, "integerish", 2, sign = "positive")
  .check_vector(original_dim, "integerish", 2, sign = "positive")

  last_col <- upper_left_indices[1] + terra::ncol(r) - 1
  last_row <- upper_left_indices[2] + terra::nrow(r) - 1
  if (last_col > original_dim[1] || last_row > original_dim[2]) {
    stop("The requested placement extends outside `original_dim`.",
         call. = FALSE)
  }

  r <- set_rcaiman_geometry(r)
  r <- terra::shift(
    r,
    dx = upper_left_indices[1] - 1,
    dy = original_dim[2] - last_row
  )
  r <- terra::extend(r, terra::ext(0, original_dim[1], 0, original_dim[2]))
  set_rcaiman_geometry(r)
}
