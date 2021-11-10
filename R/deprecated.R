#' Mask image
#'
#' Given angle restrictions, produce an image for masking fisheye photos.
#'
#'
#' @inheritParams sky_grid_segmentation
#' @param zlim Numeric vector of length two. Angles in degrees. Set the zenith
#'   angle range with inclusive limits.
#' @param alim Numeric vector of length two. Angles in degrees. Set the azimuth
#'   angle range with inclusive limits.
#'
#' @return \linkS4class{RasterLayer}
#' @export
#'
#' @family masking functions
#' @seealso \code{\link{write_bin}}
#'
#' @examples
#' z <- zenith_image(1000, lens())
#' a <- azimuth_image(z)
#'
#' m <- mask_image(z, a, c(20, 70), c(90, 180))
#' plot(m)
#'
#' m1 <- mask_image(z, a, alim = c(90, 180))
#' plot(m1)
#'
#' m2 <- mask_image(z, zlim = c(20, 70))
#' plot(m2)
#'
#' plot(m1 & m2)
#'
#' m <- mask_image(z)
#' plot(m)
mask_image <- function(z,
                       a = NULL,
                       zlim = NULL,
                       alim = NULL) {
  no_data_area <- is.na(z)

  if (all(is.null(zlim), is.null(alim))) {
    m <- !is.na(z)
  } else {
    if (all(!is.null(zlim), !is.null(alim))) {
      stopifnot(length(zlim) == 2)
      stopifnot(length(alim) == 2)

      stopifnot(all(zlim[1] >= 0, zlim[2] <= 90))
      stopifnot(all(alim[1] >= 0, alim[2] <= 360))

      z[is.na(z)] <- 0
      a[is.na(a)] <- 0
      z[z >= zlim[1] & z <= zlim[2]] <- NA
      a[a >= alim[1] & a <= alim[2]] <- NA
      m <- is.na(z) + is.na(a)
      m[m == 2] <- NA
      m <- is.na(m)
    } else {
      if (!is.null(zlim)) {
        stopifnot(length(zlim) == 2)
        stopifnot(all(zlim[1] >= 0, zlim[2] <= 90))

        z[is.na(z)] <- 0

        z[z >= zlim[1] & z <= zlim[2]] <- NA
        m <- is.na(z)
      } else {
        stopifnot(length(alim) == 2)
        stopifnot(all(alim[1] >= 0, alim[2] <= 360))

        a[is.na(a)] <- 0
        a[a >= alim[1] & a <= alim[2]] <- NA
        m <- is.na(a)
      }
    }
  }
  # fix inclusion of the area outside the circle if zmin is 0
  m[no_data_area] <- 0
  m
}

#' Choose a standard CIE sky model
#'
#' Choose one from the 15th standard CIE sky models based on sky DN values.
#'
#' The extraction of sky marks and sun coordinates is automatically done by the
#' function assuming the same that it is assumed for \code{\link{find_sky_dns}},
#' so please refer to that function help page, Details headline, first
#' paragraph.
#'
#' Since a sky grid of 30 degrees is used, the maximum number of sky marks is
#' 36, and the sun location is coarse. So, although it is not enough to estimate
#' model coefficients with \code{\link{fit_cie_sky_model}}, it can be used to
#' choose the closest standard sky.
#'
#' To know more about the standard CIE sky models, please refer to
#' \insertCite{Li2016;textual}{rcaiman}.
#'
#'
#' @inheritParams fit_cie_sky_model
#' @inheritParams extract_sky_marks
#'
#' @references \insertRef{Li2016}{rcaiman}
#'
#' @family  cie sky model functions
#'
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' my_file <- path.expand("~/DSCN5548.JPG")
#' download.file("https://osf.io/kp7rx/download", my_file,
#'               method = "auto", mode = "wb")
#' r <- read_caim(my_file,
#'                c(1280, 960) - 745,
#'                745 * 2,
#'                745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' sky <- choose_std_cie_sky_model(blue, z, a)$relative_luminance
#' plot(sky)
#' }
choose_std_cie_sky_model <- function(r, z, a, bin,
                                     use_window = TRUE,
                                     twilight = TRUE,
                                     rmse = FALSE) {
  .check_if_r_z_and_a_are_ok(r, z, a)
  compareRaster(r, bin)

  g <- sky_grid_segmentation(z, a, 30)
  sky_marks <- extract_sky_marks(r, bin, g)
  sun_coord <- extract_sun_mark(r, bin, z, a, g)$zenith_azimuth

  cells <- cellFromRowCol(a, sky_marks$row, sky_marks$col)
  sky_marks$z <- z[cells]

  if (use_window) {
    xy <-  xyFromCell(r, cells)
    sky_marks$dn <-  extract(r, xy, buffer = 1.5, fun = "mean")
  } else {
    sky_marks$dn <-  r[cells]
  }

  z_thr <- 2
  zenith_dn <- c()
  while (length(zenith_dn) < 20) {
    zenith_dn <- sky_marks$dn[sky_marks$z < z_thr]
    z_thr <- z_thr + 2
  }

  zenith_dn <- mean(zenith_dn)
  attr(zenith_dn, "max_zenith_angle") <- z_thr
  sky_marks$dn <- sky_marks$dn / zenith_dn

  path <- system.file("external", package = "rcaiman")
  skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))

  .fun <- function(i) {
    std_sky <- cie_sky_model_raster(z, a, sun_coord, as.numeric(skies[i,1:5]))

    if (use_window) {
      pred <-  extract(std_sky, xy, buffer = 1.5, fun = "mean")
    } else {
      pred <-  std_sky[cells]
    }

    return(list(std_sky_no = i,
                sun_coord = sun_coord,
                pred = pred,
                obs = sky_marks$dn,
                relative_luminance = std_sky,
                zenith_dn = zenith_dn,
                general_sky_type = skies[i,"general_sky_type"])
    )

  }

  models <- Map(.fun, 1:nrow(skies))

  if (twilight) {
    indices <- match(11:15, as.numeric(rownames(skies)))
    indices <- indices[!is.na(indices)]
    if (length(indices) != 0) {
      skies <- skies[indices,]
      civic_twilight <-  c(seq(90, 96, 1))
      for (i in seq_along(civic_twilight)) {
        sun_coord <-  c(civic_twilight[i], sun_coord[2])
        models <- c(models, suppressWarnings(Map(.fun, 1:nrow(skies))))
      }
    }
  }

  .calc_ratio_squared <- function(x) {
    ratio <- r / (x$relative_luminance * x$zenith_dn)
    sum(ratio[]^2, na.rm = TRUE)
  }

  if (!rmse) {
    error <- Map(.calc_ratio_squared, models) %>% unlist()
  } else {
    error <- Map(function(x) .calc_rmse(x$pred - x$obs), models) %>% unlist()
  }
  models[[which.min(error)]]
}


#' Out-of-the-box model-based local thresholding with CIE sky model included
#'
#' The same as \code{\link{ootb_mblt}}, it is a hard-coded version of a MBLT
#' pipeline. The code can be easily inspected by calling \code{ootb_cie_mblt}
#' --no parenthesis--, so that advanced users could use the code as a template.
#'
#' The pipeline combines \code{\link{find_sky_dns}},
#' \code{\link{fit_cone_shaped_model}}, \code{\link{choose_std_cie_sky_model}},
#' \code{\link{fit_trend_surface}}, and \code{\link{thr_image}}. The conceptual
#' design is the same that for \code{\link{ootb_mblt}}, but the working
#' binarized image produced by \code{\link{fit_cone_shaped_model}} is refined by
#' \code{\link{choose_std_cie_sky_model}}, and the standard CIE sky model is the
#' filling source for \code{\link{fit_trend_surface}}, instead of the cone
#' shaped model.
#'
#' Also, the sky actually used to obtain the binarized image is from the
#' combination of the standard CIE sky model and the trend surface, instead of
#' the cone shaped model and the trend surface.
#'
#' @inheritParams fit_cone_shaped_model
#'
#' @noRd
#' @family mblt functions
#'
#' @return Object of class list with the binarized image (named "bin") and the
#'   reconstructed skies named as follow: "sky_cs" is the cone shaped model,
#'   "std_cie_sky" is the output from \code{\link{choose_std_cie_sky_model}},
#'   "sky_s" is the trend surface, and "sky" is the actually used to obtain
#'   "bin" (please see Details).
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#' @examples
#' \dontrun{
#' my_file <- path.expand("~/DSCN5548.JPG")
#' download.file("https://osf.io/kp7rx/download", my_file,
#'               method = "auto", mode = "wb")
#' r <- read_caim(my_file,
#'                c(1280, 960) - 745,
#'                745 * 2,
#'                745 * 2)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' blue <- gbc(r$Blue)
#' bin <- ootb_cie_mblt(blue, z, a)$bin
#' plot(bin)
#' }
ootb_cie_mblt <- function(r, z, a) {
  .check_if_r_z_and_a_are_ok(r, z, a)

  bin <- find_sky_dns(r, z, a, round((360/5) * (90/5) * 0.3))

  sky_cs <- fit_cone_shaped_model(r, z, a, bin)$image
  bin <- apply_thr(r, thr_image(sky_cs, 0, 0.5))

  std_cie_sky <- choose_std_cie_sky_model(r, z, a, bin, rmse = FALSE)
  sky_cie <- std_cie_sky$relative_luminance * std_cie_sky$zenith_dn
  sky_combo <- (sky_cs + sky_cie) / 2
  suppressWarnings(bin <- apply_thr(r, thr_image(sky_combo, 0, 0.5)))

  flat <- r
  flat[] <- 0
  residu <- r - sky_combo
  residu_s <- fit_trend_surface(residu, bin, filling_source =  flat)
  sky <- sky_combo + residu_s$image
  suppressWarnings(bin <- apply_thr(r, thr_image(sky, 0, 0.5)))

  list(bin = bin,
       sky_cs = sky_cs,
       std_cie_sky = std_cie_sky,
       residu_s = residu_s,
       sky = sky)
}
