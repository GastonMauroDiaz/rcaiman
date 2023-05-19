#' Quantify the colorfulness of an image
#'
#' Quantify the colorfulness of an sRGB image using a bidimensional space formed
#' by the green/red and the blue/yellow axes of the CIE \emph{L*a*b*} space,
#' symbolized with \emph{a*} and \emph{b*}, respectively. The proposed index is
#' defined as the surface of the \emph{a*b*} plane covered by colors from the
#' image relative to the surface that the whole sRGB cube covers in the same
#' plane, expressed in percentage.
#'
#' Pixels from the image covered by pixels from \code{m} with value \code{1}
#' will be taking into account in the computations.
#'
#' If \code{plot = TRUE} is used, a plot is sent to the active graphics device.
#' It shows the color from the image plotted into a bidimensional space made by
#' the axis \emph{a*} and \emph{b*} of the CIE \emph{L*a*b* space}.
#'
#' An early version of this function was used in
#' \insertCite{Martin2020;textual}{rcaiman}.
#'
#' @inheritParams enhance_caim
#' @param m \linkS4class{SpatRaster}. A mask. For hemispherical photographs,
#'   check \code{\link{mask_hs}}. Default (\code{NULL}) is the equivalent to
#'   enter \code{!is.na(caim$Red)}.
#' @param plot Logical vector of length one. If is \code{TRUE}, a plot will be
#'   send to the graphic device, showing the data on the CIE \emph{a*b*} space.
#'
#' @return A numeric vector of length one and, if the argument \code{plot} is
#'   \code{TRUE}, an object returned by \code{\link[base]{plot}} is send to the
#'   graphic device.
#' @family Tool Functions
#' @export
#'
#' @references \insertAllCited{}
#'
#' @examples
#'
#' #Sunlight image
#' caim <- read_caim() %>% normalize()
#' #plotRGB(caim*255)
#' colorfulness(caim)
#'
#' #Diffuse-light image
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1280, 960) - 745, 745 * 2, 745 * 2) %>% normalize()
#' #plotRGB(caim*255)
#' colorfulness(caim)
#'
colorfulness <- function (caim, m = NULL, plot = FALSE)  {
  stopifnot(class(caim) == "SpatRaster")
  if (.get_max(caim) > 1)
    stop("Please check if \"caim\" was correctly normalized.")
  if (is.null(m)) m <- !is.na(caim$Red)
  .is_single_layer_raster(m, "m")
  .is_logic_and_NA_free(m)
  stopifnot(compareGeom(caim, m) == TRUE)
  stopifnot(length(plot) == 1)
  stopifnot(is.logical(plot))

  stopifnot(all(names(caim) == c("Red", "Green", "Blue")))

  .fun <- function(x) {
     colorspace::hex(colorspace::sRGB(x[, 1], x[, 2], x[, 3]))
  }

  hexs <- .fun(caim[m])
  rm(caim)
  hexs <- unique(hexs)
  rgb <- colorspace::hex2RGB(hexs)
  lab <- as(rgb, "LAB")
  lab <- colorspace::coords(lab)
  rm(rgb)
  lab <- terra::vect(lab[,2:3])

  r <- terra::rast()
  terra::ext(r) <- terra::ext(-86,98,-108,94)
  terra::res(r) <- 1

  r <- terra::rasterize(lab, r)

  # full <- expand.grid(1:256, 1:256, 1:256)/256
  # full <- sRGB(full[,1], full[,2], full[,3])
  # full <- as(full, "LAB")
  # full <- colorspace::coords(full)
  # full <- terra::vect(full[,2:3])
  # full <- terra::rasterize(full, r)
  # terra::freq(full, value = 1)[3] %>% as.numeric()
  full <- 22639

  if (plot) {
    full <- expand.grid(1:256, 1:256, 1:256)/256
    full <- sRGB(full[,1], full[,2], full[,3])
    full <- as(full, "LAB")
    full <- colorspace::coords(full)
    i <- grDevices::chull(full[,2:3])
    full <- full[i, 2:3]

    plot(full, type = "n", xlab="a*", ylab="b*", asp = 1)
    lines(rbind(full, full[1,]))
    points(lab, col = hexs)
  }
  terra::freq(r, value = 1)[3] %>% as.numeric() %>%
    magrittr::divide_by(full) %>%
    magrittr::multiply_by(100)
}

