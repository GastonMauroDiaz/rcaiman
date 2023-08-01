#' Regional thresholding
#'
#' Regional thresholding of greyscale images.
#'
#' Methods currently implemented are:
#'
#' \itemize{
#'
#' \item \strong{Diaz2018}: method presented in
#' \insertCite{Diaz2018;textual}{rcaiman} applied regionally. If this method is
#' selected, the arguments \code{intercept}, \code{slope}, and \code{prob}
#' should be provided. It works segment-wise extracting the digital numbers
#' (dns) per segment and passing them to \code{quantile(dns, prob)}, which
#' aggregated result (\code{x}) is in turn passed to \code{thr_image(x,
#' intercept, slope)}. Finally, this threshold image is applied to obtain a
#' binarized image.
#'
#' \item \strong{Methods from autothresholdr package}: this function can call
#' methods from \code{\link[autothresholdr]{auto_thresh}}. Use \code{"IsoData"}
#' to use the algorithm by \insertCite{isodata;textual}{rcaiman}, which was
#' recommended by \insertCite{Jonckheere2005;textual}{rcaiman}.
#'
#' \item \strong{Method isodata from this package}: Use \code{"thr_isodata"} to
#' use \code{\link{thr_isodata}}.
#'
#' }
#'
#' @inheritParams ootb_mblt
#' @param segmentation \linkS4class{SpatRaster}. The result of segmenting
#'   \code{r}. Probably, \code{\link{rings_segmentation}} will be the most used
#'   for fisheye images.
#' @param method Character vector of length one. See details for current
#'   options.
#' @inheritParams thr_image
#' @inheritParams fit_trend_surface
#' @param prob Numeric vector of length one. Probability for
#'   \code{\link[stats]{quantile}} calculation.
#'
#' @return An object of class \linkS4class{SpatRaster} with values \code{0} and
#'   \code{1}.
#'
#' @export
#' @family Binarization Functions
#'
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \donttest{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' r <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' blue <- gbc(r$Blue)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' rings <- rings_segmentation(z, 10)
#' bin <- regional_thresholding(blue, rings, "Diaz2018", -8, 0.5, 1)
#' plot(bin)
#' bin <- regional_thresholding(blue, rings, "thr_isodata")
#' plot(bin)
#' }
regional_thresholding <- function(r,
                                  segmentation,
                                  method,
                                  intercept = NULL,
                                  slope = NULL,
                                  prob = NULL) {
  .is_single_layer_raster(r, "r")
  .was_normalized(r)
  .is_single_layer_raster(segmentation, "segmentation")
  stopifnot(class(method) == "character")
  stopifnot(length(method) == 1)
  if (!is.null(intercept)) stopifnot(length(intercept) == 1)
  if (!is.null(slope)) stopifnot(length(slope) == 1)
  if (!is.null(prob)) stopifnot(length(prob) == 1)

  if (method == "Diaz2018") {
    if (any(is.null(intercept), is.null(slope), is.null(prob))) {
      stop("Arguments \"intercept\", \"slope\", and \"prob\" should be provided.")
    }
  }

  fun <- switch(method,
    Diaz2018 = function(dns) {
      dn <- quantile(dns, prob)
      thr_image(dn, intercept, slope)
    },
    thr_isodata = thr_isodata
  )

  if (is.null(fun)) {
    if (!requireNamespace("autothresholdr", quietly = TRUE)) {
      stop(paste(
        "Package \"autothresholdr\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
      )
    }

    fun <- function(dns) {
      autothresholdr::auto_thresh(round(dns * 255),
                                  method,
                                  ignore_black = TRUE,
                                  ignore_white = TRUE,
                                  ignore_na = TRUE)[1] / 255
    }
  }
  bin <- apply_thr(r, .get_min(r))
  .binarize_per_ring <- function(segment_id) {
    indices <- segmentation == segment_id
    thr <- fun(r[indices])
    bin[indices] <<- r[indices] > thr
  }
  segs <- unique(terra::values(segmentation)) %>% as.numeric()
  segs <- segs[!is.na(segs)]
  Map(.binarize_per_ring, segs)
  bin
}
