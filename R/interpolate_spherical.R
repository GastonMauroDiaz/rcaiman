#' Expand sky points
#'
#' Expand sky points using a _k_ nearest neighbors approach and IDW
#' interpolation. Equation 6 from Lang et al. 2006
#'
#' @inheritParams apply_thr
#' @inheritParams sky_grid_segmentation
#' @inheritParams extract_rel_radiance
#' @inheritParams sor_filter
#' @inheritParams interpolate_planar
#' @param chi_max Numeric vector of length one. Maximum spherical distance (in
#'   degrees) to search for nearest neighbors.
#' @param size Numeric vector of length one. Number of rows and colums of the
#'   planar grid that will be the skeleton for intepolating the points.
#' @param filling_source Single-layer [SpatRaster-class] or `NULL`.  A source to
#'   complement the sky pixels detected through the gaps of `r`. See
#'   [ootb_fit_cie_sky_model()].
#' @param rule Character vector of length one. Either `"any"` or `"all"`.
#' @param w Numeric vector of length one. Weight parameter for the sky model.
#' @param return_raster Logical vector of length one. If `TRUE` (default), then
#'   raster will be returned. Otherwise, a _data.frame_ representing points
#'   distributed in a regular grid will be returned.
#'
#' @returns An object of the class _data.frame_. It is the input argument
#'   `sky_points` with the following additional data:
#' \itemize{
#'   \item Rows resulting of interpolating the `sky_points` argument.
#'   \item Columns _a_, _z_, _dn_, and _initial_.
#'      \itemize{
#'        \item _a_, the azimuthal angle.
#'        \item _z_:, the zenithal angle.
#'        \item _dn_, the digital number.
#'        \item _initial_, if 'TRUE' the point is from the `sky_points` argument.
#'      }
#' }
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' r <- caim$Blue
#' bin <- apply_thr(r, thr_isodata(r[]))
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#'
#' path <- system.file("external/ootb_sky.txt", package = "rcaiman")
#' ootb_sky <- read_ootb_sky_model(gsub(".txt", "", path), z, a)
#'
#' sky_points <- ootb_sky$sky_points
#'
#' size <- 200
#' sky_points2 <- interpolate_spherical(r, z, a, sky_points,
#'                                      sky_model =  ootb_sky$sky,
#'                                      k = 10,
#'                                      p = 2,
#'                                      w = 1,
#'                                      rule = "any",
#'                                      chi_max = 20,
#'                                      size = size)
#' plot(sky)
#' }
interpolate_spherical <- function(r, z, a, sky_points, filling_source, k, p,
                                  chi_max = 20,
                                  w = 1,
                                  size = 50,
                                  rule = "any",
                                  use_window = TRUE,
                                  return_raster = TRUE
                              ) {

  .check_if_r_z_and_a_are_ok(r, z, a)
  if (is.null(filling_source)) {
    filling_source <- r
    filling_source[] <- NA
    w <- 0
  } else {
    .is_single_layer_raster(filling_source, "filling_source")
    terra::compareGeom(r, filling_source)
  }

  stopifnot(length(k) == 1)
  stopifnot(is.numeric(k))
  stopifnot(length(p) == 1)
  stopifnot(is.numeric(p))
  stopifnot(length(chi_max) == 1)
  stopifnot(is.numeric(chi_max))

  sky_points2 <- expand.grid(row = seq(1, nrow(r), length.out = size) %>% round(),
                             col = seq(1, nrow(r), length.out = size) %>% round())
  sky_points2 <- extract_dn(z, sky_points2, use_window = use_window)
  sky_points2 <- sky_points2[!is.na(sky_points2[,3]), c("row", "col")]

  # eq6 from Lang et al. 2010
  .eq6 <- function(w_k, b_k, W, M) {
    sum(c(sum(w_k * b_k), W * M), na.rm = TRUE) /
      (sum(w_k) + W)
  }

  sky_points2 <- extract_rel_radiance(r, z, a, sky_points2,
                                      no_of_points = NULL,
                                      use_window = use_window)$sky_points
  sky_points <- extract_rel_radiance(r, z, a, sky_points,
                                     no_of_points = NULL,
                                     use_window = use_window)$sky_points
  sky_points[, c("z", "a")] <- .degree2radian(sky_points[, c("z", "a")])
  sky_points2[, c("z", "a")] <- .degree2radian(sky_points2[, c("z", "a")])
  chi_max <- .degree2radian(chi_max)
  .calculate_dn <- function(i) {
    spherical_distance <- calc_spherical_distance(sky_points$z,
                                                  sky_points$a,
                                                  sky_points2[i , "z"],
                                                  sky_points2[i , "a"])
    u <- terra::cellFromRowCol(r,
                               sky_points2[i , "row"],
                               sky_points2[i,  "col"])

    sky_points <- cbind(sky_points, chi = spherical_distance)
    sky_points <- sky_points[order(spherical_distance),]
    sky_points <- sky_points[1:k, ]
    m <- sky_points$chi <= chi_max
    if (switch(rule, any = any(m), all = all(m))) {
      sky_points <- sky_points[m,]
      w_k <- 1 / sky_points$chi^p
      w_k <- w_k/max(w_k) #interpreting paragraph below eq6
      return(.eq6(w_k, sky_points[, "dn"], w, filling_source[u][,]))
    } else {
      return(filling_source[u][,])
    }
  }

  dn <- Map(.calculate_dn, seq_len(nrow(sky_points2))) %>% unlist()
  sky_points2$dn <- dn
  sky_points2 <- sky_points2[!is.na(dn),]
  initial <- c(rep(TRUE, nrow(sky_points)), rep(FALSE, nrow(sky_points2)))
  sky_points <- rbind(sky_points, sky_points2)
  sky_points <- sky_points[, -ncol(sky_points)] #remove the rl column
  sky_points[, c("z", "a")] <- .radian2degree(sky_points[, c("z", "a")])
  sky_points <- cbind(sky_points, initial)

  if (return_raster) {
    # sky <-  interpolate_planar(sky_points[!sky_points$initial, ], r,
    #                            k = 1, p = 1,
    #                            rmax = ncol(r)/size + 1,
    #                            col_id = "dn")
    # sky[is.na(z)] <- NA
    sky_points <- sky_points[!sky_points$initial, ]
    sky <- terra::rast(ncol = size, nrow = size)
    sky <- terra::rast(r)
    terra::ncol(sky) <- terra::nrow(sky) <- size
    xy <- terra::xyFromCell(r,
                            terra::cellFromRowCol(r, sky_points$row,
                                                  sky_points$col))
    cell <- terra::cellFromXY(sky, xy)
    sky[cell] <- sky_points$dn
    sky <- terra::focal(sky, 3, mean, na.policy = "only")
    sky <- terra::disagg(sky, terra::res(sky)[1], method = "bilinear")
    sky <- terra::resample(sky, r)
    sky[is.na(z)] <- NA

    names(sky) <- "Sky vault"
    return(sky)
  } else {
    return(sky_points)
  }
}
