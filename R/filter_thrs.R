#' Title
#'
#' @param thrs Character vector of length one.
#' @inheritParams sor_filter
#'
#' @returns a list.
#' @export
#'
#' @examples
#' \dontrun{
#' }
filter_thrs <- function(thrs, r, z, a,
                        g = NULL,
                        laxity = 2.5,
                        per_direction = TRUE,
                        n_min = 20) {

  if (per_direction) {

    # filter undersampled directions
    i <- terra::cellFromRowCol(r, thrs$row, thrs$col)
    n <- tapply(thrs[,3], i, length)
    thrs <- cbind(thrs, n = n[match(i, names(n))])
    thrs <- thrs[thrs$n >= n_min,]

    if (nrow(thrs) == 0) {
      stop(
        paste0(
         "No threshold value remains after filtering with ´n_min´ equal to ",
         n_min, "."
        )
      )
    }

    # Filter outliers
    i <- terra::cellFromRowCol(r, thrs$row, thrs$col)
    M <- tapply(thrs[,3], i, stats::median)
    MAD <- tapply(thrs[,3], i, stats::mad)
    thrs <- cbind(thrs, M = M[match(i, names(M))], MAD = MAD[match(i, names(M))])
    deviation <- (thrs[,3] - thrs$M) / thrs$MAD
    i <- abs(deviation) <= laxity
    i[is.na(i)] <- TRUE
    thrs <- thrs[i,]


  } else {
    M <- stats::median(thrs[,3])
    MAD <- stats::mad(thrs[,3])
    deviation <- (thrs[,3] - M) / MAD
    i <- abs(deviation) <= laxity
    i[is.na(i)] <- TRUE
    thrs <- thrs[i,]

  }

  # Built rr
  if (ncol(thrs) == 3) thrs <- cbind(thrs, n = 1)
  i <- terra::cellFromRowCol(r, thrs$row, thrs$col)
  row <- tapply(thrs[, "row"], i, max)
  col <- tapply(thrs[, "col"], i, max)
  n <- tapply(thrs[, "n"], i, max)
  dn <- tapply(thrs[,3], i,  mean) # average <<<<<<<< IMPORTANT!!!! <<<<<<<
  sky_points <- data.frame(row, col, dn, n)

  if (!is.null(g)) {
    cell_ids <- extract_dn(g,
                           sky_points[, c("row", "col")],
                           use_window = FALSE)[,3]
    .sample_per_cell <- function(cell_id) {
      sky_points <- sky_points[cell_ids == cell_id,]
      row.names(sky_points[sample(1:nrow(sky_points),
                                  size = 1,
                                  prob = sky_points$n), ]
                )
    }
    row_names <- Map(.sample_per_cell, unique(cell_ids)) %>% unlist()
    sky_points <- sky_points[row_names,]
  }

  r_dummy <- interpolate_planar(sky_points, r, k = 1, p = 1, rmax = 3,
                                col_id = "dn")

  rr <- extract_rel_radiance(r_dummy, z, a,
                             sky_points[, c("row", "col")],
                             no_of_points = 3,
                             use_window = FALSE)

  list(thrs = thrs,
       rr = rr)
}
