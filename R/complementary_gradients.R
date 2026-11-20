#' Calculate complementary gradients
#'
#' @description
#' Compute three color-opponent gradients to enhance the visual separation
#' between sky and canopy in photographs, particularly under
#' non-diffuse light or complex cloud patterns.
#'
#' @details
#' Inspired by opponent processing in human color vision, this function
#' combines spectral bands into opposing chromatic signals designed to
#' distinguish sky from canopy. The gradients are tailored to canopy
#' photography rather than intended to reproduce human visual responses.
#'
#' Use the four bands returned by [split_cfa()], with the names described
#' in `caim`. No mapping to RGB is needed. Layers are selected by name,
#' independently of their order or position in the CFA pattern.
#'
#' For `cfa_type = "bayer"`, it averages the two green bands:
#' \eqn{G = (G_1 + G_2)/2}. The gradients are:
#'
#' - `"green_magenta"` = \eqn{(R - G + B) / (R + G + B)}
#' - `"yellow_blue"`   = \eqn{(-R - G + B) / (R + G + B)}
#' - `"red_cyan"`      = \eqn{(-R + G + B) / (R + G + B)}
#'
#' For `cfa_type = "cygm"`, it uses Cyan (C), Yellow (Y), Green (G), and
#' Magenta (M) directly:
#'
#' - `"green_magenta"` = \eqn{(M - G) / (M + G)}
#' - `"yellow_blue"`   = \eqn{(M + C - Y) / (M + C + Y)}
#' - `"red_cyan"`      = \eqn{(C - M - Y) / (C + M + Y)}
#'
#' Higher values favor magenta over green, blue over yellow, and cyan over
#' red, respectively, for both CFA types. These contrasts are intended to
#' place the sky toward higher values, consistently with brightness.
#' The calculation is pixelwise: multiplying all bands by the same positive
#' factor leaves the gradients unchanged. CYGM combinations define chromatic
#' contrasts; they do not reconstruct missing RGB bands. The same polarity
#' does not imply numerically interchangeable gradients between CFA types.
#'
#' If all four input bands have the same positive value, Bayer gradients
#' equal 1/3, -1/3, and 1/3 in the order listed above; CYGM gradients equal
#' 0, 1/3, and -1/3. Thus, even when white is represented by equal band
#' values, it does not generally coincide with zero. Zero is not defined
#' as a common achromatic reference for these gradients.
#'
#' Equal band responses to white light should not be assumed for RAW data
#' from a photographic system without a characterization of its relative band
#' responses. This is a limitation to consider when using photographic
#' equipment as a measuring instrument.
#'
#' @note
#' This function is part of a paper under preparation.
#'
#' @param caim numeric [terra::SpatRaster-class] with exactly four uniquely
#'   named layers: `"Red"`, `"Green1"`, `"Green2"`, and `"Blue"` for Bayer;
#'   `"Cyan"`, `"Yellow"`, `"Green"`, and `"Magenta"` for CYGM. Layer order
#'   is unrestricted. Declare these names in the `cfa_pattern` supplied to
#'   [split_cfa()]. The two Bayer green positions may be named `"Green1"`
#'   and `"Green2"` in either order because their average is used.
#'   Values are sensor-band measurements before RGB mapping. Their relative
#'   scaling affects the gradients; no white balance is estimated or applied.
#' @param cfa_type character vector of length one: `"bayer"` (default) or
#'   `"cygm"`. Selects the expected band names and gradient formulas.
#'
#' @return Numeric [terra::SpatRaster-class] with three layers and the same
#'   geometry as `caim`. The layers (`"green_magenta"`, `"yellow_blue"`,
#'   `"red_cyan"`) are normalized chromatic gradients. A missing value in a
#'   band used by a gradient, or a zero denominator, produces `NA` in that
#'   gradient. For Bayer, all four bands contribute to every gradient.
#'   For CYGM, `green_magenta` uses only Green and Magenta; the other two
#'   gradients use only Cyan, Yellow, and Magenta.
#'
#' @export
#'
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g)
#' sun_angles <- data.frame(zenith = 45 * pi/180, azimuth = 10 * pi/180)
#' m <- !is.na(d$zenith)
#'
#' # Construct synthetic Bayer bands from three different CIE distributions.
#' # Both green positions receive the same values in this example.
#' # Independent rescaling creates the colors of this example; it is not
#' # a required preprocessing step for photographs.
#' r <- cie_image(d, sun_angles, as.numeric(cie_table[12, 1:5]), m = m) %>%
#'   normalize_min_max()
#' green <- cie_image(d, sun_angles, as.numeric(cie_table[4, 1:5]), m = m) %>%
#'   normalize_min_max()
#' b <- cie_image(d, sun_angles, as.numeric(cie_table[1, 1:5]), m = m) %>%
#'   normalize_min_max()
#' caim <- c(r, green, green, b)
#' names(caim) <- c("Red", "Green1", "Green2", "Blue")
#' terra::plotRGB(caim, r = 1, g = 2, b = 4, scale = 1, stretch = NULL)
#'
#' com <- complementary_gradients(caim, cfa_type = "bayer")
#' plot(com)
#'
#' # A synthetic CYGM example: each row supplies four sensor-band responses.
#' cygm <- set_rcaiman_geometry(rast(nrows = 2, ncols = 2, nlyrs = 4))
#' names(cygm) <- c("Cyan", "Yellow", "Green", "Magenta")
#' values(cygm) <- rbind(c(8, 2, 2, 8), c(2, 8, 8, 2),
#'                       c(8, 2, 2, 2), c(5, 5, 5, 5))
#' complementary_gradients(cygm, cfa_type = "cygm")[]
complementary_gradients <- function(caim, cfa_type = c("bayer", "cygm")) {
  cfa_type <- match.arg(cfa_type)
  .assert_spatraster(caim)
  required <- switch(cfa_type,
    bayer = c("Red", "Green1", "Green2", "Blue"),
    cygm = c("Cyan", "Yellow", "Green", "Magenta"))
  if (terra::nlyr(caim) != 4L || anyDuplicated(names(caim)) ||
      !setequal(names(caim), required)) {
    stop("For `cfa_type = \"", cfa_type, "\"`, `caim` must have exactly ",
         "four uniquely named layers: ", paste(required, collapse = ", "),
         ".", call. = FALSE)
  }
  if (cfa_type == "bayer") {
    R <- caim$Red
    G <- (caim$Green1 + caim$Green2) / 2
    B <- caim$Blue
    brightness <- R + G + B
    brightness[brightness == 0] <- NA
    green_magenta <- (R - G + B) / brightness
    yellow_blue <- (-R - G + B) / brightness
    red_cyan <- (-R + G + B) / brightness
  } else {
    C <- caim$Cyan
    Y <- caim$Yellow
    G <- caim$Green
    M <- caim$Magenta
    green_magenta_denominator <- M + G
    green_magenta_denominator[green_magenta_denominator == 0] <- NA
    brightness <- M + C + Y
    brightness[brightness == 0] <- NA
    green_magenta <- (M - G) / green_magenta_denominator
    yellow_blue <- (M + C - Y) / brightness
    red_cyan <- (C - M - Y) / brightness
  }

  names(green_magenta) <- "green_magenta"
  names(yellow_blue) <- "yellow_blue"
  names(red_cyan) <- "red_cyan"
  com <- c(green_magenta, yellow_blue, red_cyan)
  com
}

