#' Extract radiometry
#'
#' Extract radiometry
#'
#' Consult this
#' \href{https://docs.google.com/document/d/178yZDAcfx--Xn1Ye8Js-kUXuPCuYOHQL5fxAH7KBEoY/edit?usp=sharing}{document}
#' for additional details.
#'
#' This code exemplify how to use the function to obtain the data and use base R
#' functions to obtain the vignetting function.
#' \preformatted{
#' up_data <- extract_radiometry("D610/up")
#' down_data <- extract_radiometry("D610/down")
#' right_data <- extract_radiometry("D610/right")
#' left_data <- extract_radiometry("D610/left")
#'
#' ds <- rbind(up_data, down_data, right_data, left_data)
#'
#' plot(ds, xlim = c(0, pi/2), ylim= c(0.5,1.05),
#'       col = c(rep(1,9),rep(2,9),rep(3,9),rep(4,9)))
#' legend("bottomleft", legend = c("up", "down", "right", "left"),
#'        col = 1:4, pch = 1)
#'
#' fit <- lm((1 - ds$radiometry) ~ poly(ds$theta, 3, raw = TRUE) - 1)
#' summary(fit)
#' coef <- fit$coefficients
#' fun <- function(x) 1 - (coef[1] * x + coef[2] * x^2 + coef[3] * x^3)
#' curve(fun, add = TRUE, col= 2)
#' coef
#' }
#'
#' @param path_to_pgms Character vector of length one. Path to a folder with the
#'   PGM files for radiometry sampling. These files must comply with the
#'   equidistant projection.
#' @param z_thr Numeric vector of length one. Maximum angle in degrees used to
#'   calculate the digital number (DN) not affected by the vignetting effect. If
#'   \code{NULL} is provided, DNs are returned instead of relative luminance.
#' @param size_px Numeric vector of length one. Diameter in pixels of the
#'   circular sampling area at the image center. This area is modified
#'   considering the equidistant projection distortion. Therefore, it will be
#'   visualized as an ellipse any other place but the image center.
#' @inheritParams base::dir
#'
#' @family Lens Functions
#'
#' @return An object from the class \code{data.frame} with columns \emph{theta}
#'   (zenith angle in radians) and \emph{radiometry} (relative luminance or
#'   digital number, depending on argument \code{z_thr}.
#' @export
#'
extract_radiometry <- function(path_to_pgms,
                               z_thr = 15,
                               size_px = NULL,
                               pattern = "pgm") {
  .this_requires_conicfit()

  files <- dir(path_to_pgms, full.names = TRUE, pattern = pattern)

  for (i in seq_along(files)) {
    if (i == 1) {
      r <- read_caim(files[i])
    } else {
      r <- c(r, read_caim(files[i]))
    }
  }
  r <- max(r)
  z <- zenith_image(ncol(r), lens())
  a <- azimuth_image(z)

  if (is.null(size_px)) {
    size_px <- round(terra::ncol(r) * 0.02)
  }
  .size_fun <- function(theta, size_px) {
    theta <- theta*pi/180
    size_px * (theta/sin(theta))
  }

  message("Look for instructions on the popup window")

  grDevices::x11()
  plot(r, legend = FALSE, col = grDevices::grey((1:255)/255))
  .show_popup(
   "Please zoom to the area of interest by making a rectangle with two points."
  )
  r <- terra::crop(r, terra::draw())
  z <- terra::crop(z, r)
  a <- terra::crop(a, r)
  if ((terra::ncol(r) / terra::nrow(r)) < 1) {
    r <- terra::t(r)
    z <- terra::t(z)
    a <- terra::t(a)
    .make_elli <- function(x, y, theta, phi, size_px) {
      conicfit::calculateEllipse(x, y,
                                 size_px, .size_fun(theta, size_px),
                                 360 - phi, 50)
    }
  } else {
    .make_elli <- function(x, y, theta, phi, size_px) {
      conicfit::calculateEllipse(x, y,
                                 .size_fun(theta, size_px), size_px,
                                 phi, 50)
    }
  }
  plot(r, legend = FALSE, col = grDevices::grey((1:255)/255))

  .show_popup("Please, mark one point in each central cross.")

  coords <- terra::click(n = 1)
  points(coords, pch = 4, col = "red", cex = 1)

  for (i in 2:length(files)) {
    foo <- terra::click(n = 1)
    points(foo, pch = 4, col = "red", cex = 1)
    coords <- rbind(coords, foo)
  }


  img_points <- data.frame(col = terra::colFromX(r, coords[,1]),
                           row = terra::rowFromY(r, coords[,2]))

  theta <- extract_dn(z, img_points)[,3]
  phi <- extract_dn(a, img_points)[,3]

  z <- matrix(ncol = 5) %>% as.data.frame()
  colnames(z) <- c("object", "part", "x", "y", "hole")
  for (i in seq_along(img_points$row)) {
    elli <- .make_elli(coords[i,1], coords[i,2], theta[i], phi[i], size_px)
    elli <- cbind(i, 1, elli, 0)
    elli2 <- .make_elli(coords[i,1], coords[i,2], theta[i], phi[i], size_px/2)
    elli2 <- cbind(i, 1, elli2, 1)
    colnames(elli) <- c("object", "part", "x", "y", "hole")
    if (i == 1) {
      # z <- elli
      z <- rbind(elli, elli2)
    } else {
      # z <- rbind(z, elli)
      z <- rbind(z, elli, elli2)
    }
  }

  p <- vect(z, "polygons")
  plot(p, add = TRUE)
  dns <- terra::extract(r, p, fun = median)
  if (is.null(z_thr)) {
    radiometry <- dns[, 2]
  } else {
    radiometry <- dns[, 2]/mean(dns[theta < z_thr, 2])
  }
  theta <- theta*pi/180
  data.frame(theta, radiometry)
}
