#' Quantify colorfulness
#'
#' Quantify the colorfulness of an image
#'
#' Quantify the colorfulness of an sRGB image using a bidimensional space formed
#' by the green/red and the blue/yellow axes of the \emph{CIE LAB} space,
#' symbolized with \emph{A} and \emph{B}, respectively. The colorfulness index
#' (CI) is defined as
#'
#' \eqn{CI = \dfrac{A_o}{A_p} \cdot 100},
#'
#' where \eqn{A_o} and \eqn{A_p} are the observed and potential area of the
#' \emph{AB} plane. \eqn{A_o} refers to the colors from the image while
#' \eqn{A_p} to the colors from the whole sRGB cube.
#'
#' @note An early version of this function was used in
#'   \insertCite{Martin2020;textual}{rcaiman}.
#'
#' @inheritParams enhance_caim
#'
#' @return A numeric vector of length one.
#'
#' @family Tool Functions
#' @export
#'
#' @references \insertAllCited{}
#'
#' @examples
#' caim <- read_caim() %>% normalize()
#' colorfulness(caim)
colorfulness <- function (caim, m = NULL)  {
  stopifnot(class(caim) == "SpatRaster")
  if (.get_max(caim) > 1)
    stop("Please check if \"caim\" was correctly normalized.")
  if (is.null(m)) m <- !is.na(caim$Red)
  .is_single_layer_raster(m, "m")
  .is_logic_and_NA_free(m)
  stopifnot(compareGeom(caim, m) == TRUE)

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

  if (FALSE) {
    full <- expand.grid(1:256, 1:256, 1:256)/256
    full <- sRGB(full[,1], full[,2], full[,3])
    full <- as(full, "LAB")
    full <- colorspace::coords(full)
    i <- grDevices::chull(full[,2:3])
    full <- full[i, 2:3]

    plot(full, type = "n", xlab="a", ylab="b", asp = 1)
    lines(rbind(full, full[1,]))
    points(lab, col = hexs)
  }
  terra::freq(r, value = 1)[3] %>% as.numeric() %>%
    magrittr::divide_by(full) %>%
    magrittr::multiply_by(100)
}

