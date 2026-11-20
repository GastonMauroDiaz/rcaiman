#' Aggregate point observations by segment
#'
#' @description
#' Locate continuous image coordinates in a segmentation raster and summarize
#' their measured values separately for each segment.
#'
#' @param sampling_coord data frame with numeric `x` and `y_down` columns in
#'   the image-coordinate convention of [coord_from_indices()], and a numeric
#'   or logical column of observations. Coordinates must refer to the grid of
#'   `seg`. Additional columns are ignored.
#' @param seg single-layer [terra::SpatRaster-class] containing nonnegative
#'   integer segment labels. A label can occupy multiple cells, including
#'   disconnected regions. Missing labels are excluded.
#'   A logical mask (see [direction_image()] for examples) is also accepted:
#'   `TRUE` is label 1 and `FALSE` is label 0. All selected pixels belong to
#'   one segment, even when disconnected. Use [label_components()] if connected
#'   regions should have separate labels.
#' @param col_id character or positive integer-like numeric vector of length
#'   one. Name or position of the observation column. Defaults to `"dn"`.
#' @param fun function applied to the observations within each segment.
#'   Must return one numeric or logical value, possibly missing. Defaults
#'   to `mean`. Named scalar results, such as those from `quantile`, are allowed.
#' @param return character vector of length one: `"raster"` (default) or
#'   `"vector"`, following [compute_mean()].
#' @param ignore_label_0 logical vector of length one. If `TRUE` (default),
#'   observations on cells labeled zero are ignored.
#' @param ... Additional arguments passed to `fun`, such as `na.rm = TRUE`
#'   or `probs = 0.99`.
#'
#' @details
#' Each coordinate is assigned to its containing cell using
#' [indices_from_coord()], and that cell's label determines the segment.
#' No nearest-neighbor search or interpolation is performed. Image coordinates
#' start at the upper-left corner, independently of the raster's stored extent.
#' Cells include their upper and left boundaries, excluding their lower and
#' right boundaries. Consequently, `x = ncol(seg)` and `y_down = nrow(seg)`
#' lie outside the image. Missing or nonfinite coordinates and positions
#' outside the image do not contribute.
#'
#' Every retained row is an observation. Repeated positions and multiple
#' observations in the same cell all contribute; observations are not first
#' averaged by cell. Original row order is retained within each segment.
#' Missing or nonfinite measured values are passed to `fun` unchanged: choose
#' their treatment through the function and its arguments. A segment with no
#' observations receives `NA`, without calling `fun` on an empty vector.
#'
#' [compute_mean()] summarizes values already stored on a raster. This function
#' instead summarizes values in a table whose positions identify their segments.
#' To aggregate separately within each raster cell, give every cell a distinct
#' label; see [panoramic()] for an angular aggregation workflow.
#'
#' @return With `return = "raster"`, a single-layer raster on the grid of
#'   `seg`, named after the observation column. Each segment receives its
#'   aggregate value; excluded cells and segments without observations are `NA`.
#'   With `return = "vector"`, a numeric vector named by segment label, sorted
#'   in increasing label order. It includes every nonmissing label in `seg`,
#'   except zero when ignored, with `NA` for segments without observations.
#'   Logical summaries are represented as 0 or 1 in either output.
#' @seealso [compute_mean()], [feature_image()], [panoramic()]
#' @export
#' @examples
#' seg <- set_rcaiman_geometry(terra::rast(nrows = 2, ncols = 3))
#' terra::values(seg) <- c(1, 1, 2, 1, 2, 3)
#' samples <- data.frame(x = c(0.2, 1.4, 2.6, 1.8),
#'                       y_down = c(0.3, 0.7, 0.4, 1.2),
#'                       dn = c(2, 6, 10, 20))
#' aggregate_by_segment(samples, seg)
#' aggregate_by_segment(samples, seg, fun = quantile, probs = 0.99,
#'                      return = "vector")
aggregate_by_segment <- function(sampling_coord, seg, col_id = "dn",
                                 fun = mean, return = "raster",
                                 ignore_label_0 = TRUE, ...) {
  .check_columns(sampling_coord, required_columns = c("x", "y_down"))
  if (anyDuplicated(names(sampling_coord))) {
    stop("`sampling_coord` must have unique column names.", call. = FALSE)
  }
  if (is.numeric(col_id)) .check_vector(col_id, "integerish", 1, sign = "positive")
  col_id <- .resolve_required_column(col_id, sampling_coord)
  .assert_single_layer(seg)
  .assert_choice(return, c("raster", "vector"))
  .check_vector(ignore_label_0, "logical", 1)
  fun <- match.fun(fun)
  for (column in c("x", "y_down")) {
    .check_vector(sampling_coord[[column]], "numeric", allow_na = TRUE,
                  allow_inf = TRUE, name = column)
  }
  observations <- sampling_coord[[col_id]]
  if ((!is.numeric(observations) && !is.logical(observations)) ||
      is.complex(observations) || !is.null(dim(observations))) {
    stop("The observation column must be a numeric or logical vector.", call. = FALSE)
  }
  labels <- .segment_values(seg)
  ids <- sort(unique(labels[!is.na(labels) & (!ignore_label_0 | labels != 0)]))
  result <- stats::setNames(rep(NA_real_, length(ids)), format(ids, scientific = FALSE, trim = TRUE))

  x <- sampling_coord$x
  y <- sampling_coord$y_down
  inside <- is.finite(x) & is.finite(y) & x >= 0 & x < ncol(seg) &
    y >= 0 & y < nrow(seg)
  indices <- indices_from_coord(x[inside], y[inside])
  cells <- terra::cellFromRowCol(seg, indices$row, indices$col)
  group <- match(labels[cells], ids)
  groups <- split(observations[inside], group)
  for (name in names(groups)) {
    value <- fun(groups[[name]], ...)
    if (length(value) != 1L || (!is.numeric(value) && !is.logical(value)) ||
        is.complex(value)) {
      stop("`fun` must return one numeric or logical value per segment ",
           "(label ", ids[as.integer(name)], ").", call. = FALSE)
    }
    result[as.integer(name)] <- as.numeric(value)
  }
  if (return == "vector") return(result)
  out <- terra::rast(seg)
  terra::values(out) <- unname(result[match(labels, ids)])
  names(out) <- col_id
  out
}
