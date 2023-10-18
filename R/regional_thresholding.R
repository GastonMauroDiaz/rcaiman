#' Regional thresholding
#'
#' Regional thresholding of greyscale images.
#'
#' Methods currently implemented are:
#'
#' * __Diaz2018__: method presented in
#' \insertCite{Diaz2018;textual}{rcaiman} applied regionally. If this method is
#' selected, the arguments `intercept`, `slope`, and `prob` should be provided.
#' It works segment-wise extracting the digital numbers per segment and passing
#' them to [stats::quantile()] along with `prob`, which aggregated result is in
#' turn passed to [thr_mblt()] along with `intercept` and `slope`. Finally, this
#' threshold image is applied to obtain a binarized image.
#' * __Methods from autothresholdr package__: this function can call
#' methods from [autothresholdr::auto_thresh()]. For instance, use `"IsoData"`
#' to use the algorithm by \insertCite{isodata;textual}{rcaiman}, which was
#' recommended by \insertCite{Jonckheere2005;textual}{rcaiman}.
#' * __Method isodata from this package__: Use `"thr_isodata"` to
#' use [thr_isodata()].
#'
#' @inheritParams ootb_mblt
#' @param segmentation [SpatRaster-class]. The result of segmenting `r`.
#'   Arguably, the result of a call to [rings_segmentation()] will be the
#'   preferred choice for fisheye images.
#' @param method Character vector of length one. See details for current
#'   options.
#' @inheritParams thr_mblt
#' @inheritParams fit_trend_surface
#' @param prob Numeric vector of length one. Probability for [stats::quantile()]
#'   calculation.
#'
#' @return An object of class [SpatRaster-class] with values `0` and `1`.
#'
#' @export
#' @family Binarization Functions
#'
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, c(1250, 1020) - 745, 745 * 2, 745 * 2)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' r <- gbc(caim$Blue)
#' r <- correct_vignetting(r, z, c(0.0638, -0.101)) %>% normalize()
#' rings <- rings_segmentation(z, 15)
#' bin <- regional_thresholding(r, rings, "Diaz2018", -7.8, 0.95 * 0.5, 0.99)
#' plot(bin)
#' bin <- regional_thresholding(r, rings, "thr_isodata")
#' plot(bin)
#' #' }
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
      thr_mblt(dn, intercept, slope)
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
