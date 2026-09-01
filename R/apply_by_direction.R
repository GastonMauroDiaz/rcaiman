#' Apply methods by direction using a constant field of view
#'
#' Apply methods to each set of pixels defined by a direction and a constant
#' field of view (FOV).
#'
#' @details
#' Neighboring pixels are selected around a set of query directions distributed
#' over the hemispherical domain. For each directional neighborhood, `fun` is
#' evaluated using either:
#'
#' \describe{
#'   \item{`method = "vector"`}{
#'   argument `fun` receives a numeric vector extracted from neighboring pixels.
#'   }
#'
#'   \item{`method = "single_layer_image"`}{
#'   `fun` receives a temporary single-layer raster and its corresponding
#'   logical mask identifying valid pixels.
#'   }
#' }
#'
#' The function supplied to `fun` must return either a single numeric value or
#' a named numeric vector. When multiple values are returned, each element is
#' processed independently and stored as a separate raster layer in the output.
#'
#' Multiple directional estimates may contribute to the same raster cell.
#' Contributions are grouped by cell, filtered using a robust MAD-based outlier
#' criterion, and averaged to obtain the final pixel value. An additional layer
#' named `"n"` stores the number of directional estimates associated with each
#' cell prior to outlier removal.
#'
#' To calculates sphercial distances, directions are first mapped to unit
#' Cartesian vectors and distances are computed as the arc-cosine of the
#' dot-product. For details, see [calc_spherical_distance()].
#'
#' @param spacing numeric vector of length one. Angular spacing (in degrees)
#'   between directions to process.
#' @param fov numeric vector of length one. Field of view in degrees.
#' @param method character vector of length one. Evaluation framework used to
#'   apply `fun`. Must be `"vector"` or `"single_layer_image"`.
#' @param fun function used to compute one or more directional estimators.
#'   For `method = "vector"`, the expected signature is `function(x)`; for
#'   `method = "single_layer_image"`, `function(r, m)`.
#'
#' @inheritParams binarize_with_thr
#' @inheritParams add_directions
#' @inheritParams compute_canopy_openness
#' @inheritParams rem_outliers
#'
#' @return Numeric [terra::SpatRaster-class].
#'
#' Output layers correspond to the directional estimators returned by `fun`.
#' An additional layer named `"n"` contains the number of valid directional
#' estimates contributing to each pixel.
#'
#' @note
#' This function is part of a manuscript currently under preparation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' thr <- apply_by_direction(r, z, a, m, spacing = 10, fov = 40)
#' plot(thr$v1)
#' plot(thr$n)
#' bin <- binarize_with_thr(r, thr$v1)
#' plot(bin)
#'
#' thr <- apply_by_direction(r, z, a, m, spacing = 10, fov = 40,
#'                           method = "vector",
#'                           fun = function(x) {
#'                             thr_two_corner(x)$tm
#'                           })
#' plot(thr$v1)
#' bin <- binarize_with_thr(r, thr$v1)
#' plot(bin)
#'
#' my_fun <- function(x) {
#'   m <- !is.na(x)
#'   autothresholdr::auto_thresh(x[m], "Otsu")[1]
#' }
#' nr <- normalize_min_max(r) %>% multiply_by(255) %>% round()
#' my_fun(nr)
#' thr <- apply_by_direction(nr, z, a, m, spacing = 10, fov = 40,
#'                           method = "vector",
#'                           fun = my_fun)
#' plot(thr$v1)
#' bin <- binarize_with_thr(r, thr$v1)
#' plot(bin)
#'
#' # Severals at once
#' my_fun <- function(x) {
#'   m <- !is.na(x)
#'   c(Huang = autothresholdr::auto_thresh(x[m], "Huang")[1],
#'     Moments = autothresholdr::auto_thresh(x[m], "Moments")[1],
#'     Otsu = autothresholdr::auto_thresh(x[m], "Otsu")[1])
#' }
#' nr <- normalize_min_max(r) %>% multiply_by(255) %>% round()
#' thr <- apply_by_direction(nr, z, a, m, spacing = 10, fov = 40,
#'                           method = "vector",
#'                           fun = my_fun)
#' plot(thr)
#' bin <- lapply(seq(1, nlyr(thr) - 1), function(i) {
#'   bin <- binarize_with_thr(r, thr[[i]])
#'   names(bin) <- names(thr)[i]
#'   bin
#' })
#' plot(rast(bin))
#'
#'
#' # A method that requires images
#' my_fun <- function(r, m, ...) {
#'   r[!m] <- 0
#'   c(
#'     sw0.1 = thr_nbhu(
#'       r,
#'       max(round(min(r[m])), 10),
#'       round(max(r[m])),
#'       support_weight = 0.1,
#'       smooth_curve = TRUE
#'     ),
#'     sw0.2 = thr_nbhu(
#'       r,
#'       max(round(min(r[m])), 10),
#'       round(max(r[m])),
#'       support_weight = 0.2,
#'       smooth_curve = TRUE
#'     ),
#'     sw0.3 = thr_nbhu(
#'       r,
#'       max(round(min(r[m])), 10),
#'       round(max(r[m])),
#'       support_weight = 0.3,
#'       smooth_curve = TRUE
#'     ),
#'     sw0.4 = thr_nbhu(
#'       r,
#'       max(round(min(r[m])), 10),
#'       round(max(r[m])),
#'       support_weight = 0.4,
#'       smooth_curve = TRUE
#'     )
#'   )
#' }
#' thr <- apply_by_direction(r, z, a, m, spacing = 10, fov = 40,
#'                           method = "single_layer_image",
#'                           fun = my_fun)
#' plot(thr)
#' bin <- lapply(seq(1, nlyr(thr) - 1), function(i) {
#'   bin <- binarize_with_thr(r, thr[[i]])
#'   names(bin) <- names(thr)[i]
#'   bin
#' })
#' plot(rast(bin))
#'
#' }
apply_by_direction <- function(r, z, a, m,
                               spacing = 20,
                               laxity = 2.5,
                               fov = 40,
                               method = "vector",
                               fun = thr_isodata
) {

  .vector <- function(df, ...) {

    v <- df[, r_name]

    thr <- try(fun(v), silent = TRUE)

    if (inherits(thr, "try-error"))
      return(NULL)
    if (!is.numeric(thr))
      return(NULL)

    if (is.null(names(thr))) {
      names(thr) <- paste0("v", seq_along(thr))
    }

    do.call(rbind, lapply(seq_along(thr), function(i) {
      data.frame(
        cell = df$idx,
        variable = rep(names(thr)[i], nrow(df)),
        value = rep(unname(thr[i]), nrow(df)),
        row.names = NULL
      )
    }))
  }

  .single_layer <- function(df, ...) {
    r_local <- terra::rast(nrow = diff(range(df$row)) + 1,
                           ncol = diff(range(df$col)) + 1,
                           nlyrs = 1)

    i <- terra::cellFromRowCol(r_local, df$row - min(df$row) + 1, df$col - min(df$col) + 1)

    if (length(i) == 0)
      return(NULL)

    r_local[i] <- df[, r_name]
    r_local <- set_rcaiman_geometry(r_local)
    names(r_local) <- r_name

    m_local <- !is.na(r_local)

    thr <- try(fun(r_local, m_local), silent = TRUE)

    if (inherits(thr, "try-error"))
      return(NULL)
    if (!is.numeric(thr))
      return(NULL)

    if (is.null(names(thr))) {
      names(thr) <- paste0("v", seq_along(thr))
    }

    do.call(rbind, lapply(seq_along(thr), function(i) {
      data.frame(
        cell = df$idx,
        variable = rep(names(thr)[i], nrow(df)),
        value = rep(unname(thr[i]), nrow(df)),
        row.names = NULL
      )
    }))
  }

  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .assert_single_layer(r, allow_na = FALSE)
  .assert_same_geom(r, z)
  .assert_logical_mask(m)
  .assert_same_geom(m, z)

  .check_vector(spacing, "numeric", 1, sign = "positive")
  .check_vector(laxity, "numeric", 1, sign = "positive")
  .check_vector(fov, "numeric", 1, sign = "positive")

  .assert_choice(method, c(
    "vector",
    "single_layer_image"
  ))

  .fun <- switch(
    method,
    vector = .vector,
    single_layer_image = .single_layer
  )

  r_name <- names(r)

  idx <- rast(r[[1]])
  idx[] <- 1:ncell(r)

  names(idx) <- "idx"
  names(m) <- "m"

  sampling_points <- raster_grid_in_circle(z, a, 1)
  sampling_points <- extract_dn(m, sampling_points, use_window = FALSE)
  sampling_points <- sampling_points[sampling_points$m == 1, ]
  sampling_points <- extract_dn(r, sampling_points, use_window = FALSE)
  sampling_points <- extract_dn(idx, sampling_points, use_window = FALSE)

  search_radius <- .degree2radian(fov / 2)

  xyz <- .directions_to_cartesian(sampling_points)

  qp <- .directions_to_cartesian(fibonacci_points(z, a, spacing))

  output <- lapply(seq_len(nrow(qp)), function(i) {
    xyz0 <- c(x = qp$x[i], y = qp$y[i], z = qp$z[i])
    dot <- xyz$x * xyz0["x"] + xyz$y * xyz0["y"] + xyz$z * xyz0["z"]
    dot <- pmin(pmax(dot, -1), 1)

    d <- acos(dot)
    idx <- seq_along(d)[d <= search_radius]

    .fun(sampling_points[idx, ], d[idx], i)
  })

  output <- output[!vapply(output, is.null, logical(1))]

  acc <- do.call(rbind, output)
  acc <- acc[!is.na(acc$value), ]

  splitted <- split(acc, acc$variable)

  results <- lapply(splitted, function(x) {
    robust_grouped_mean_cpp(cell = x$cell,
                            value = x$value,
                            laxity = laxity)

  })

  sky_list <- lapply(results, function(res) {
    rr <- rast(m)
    rr[res$cell] <- res$value
    rr
  })

  n <- rast(m)
  first_result <- results[[1]]
  n[first_result$cell] <- first_result$n
  n[is.na(n)] <- 0
  sky <- c(rast(sky_list), n)
  names(sky) <- c(names(results), "n")
  sky
}
