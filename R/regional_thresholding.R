#' Regional thresholding
#'
#' Regional thresholding of greyscale images
#'
#' Methods currently implemented are:
#'
#' \itemize{ \item \strong{Diaz2018}: method presented in
#' \insertCite{Diaz2018;textual}{rcaiman} applied regionally. If this method is
#' selected, the arguments \code{intercept}, \code{slope}, and \code{prob}
#' should be provided. It works segmentwise extracting the digital numbers (dns)
#' per segment and passing them to \code{quantile(dns, prob)}, which aggregated
#' result (x) is in turn passed to \code{thr_image(x, intercept, slope)}.
#' Finally, this threshold image is applied to obtain a binarized image.
#'
#' \item \strong{Methods from autothresholdr package}: this function can call
#' methods from \code{\link[autothresholdr]{auto_thresh}}. Use \code{"IsoData"}
#' to use the algorithm by \insertCite{isodata;textual}{rcaiman}, which is the one
#' recommended by \insertCite{Jonckheere2005;textual}{rcaiman}. }
#'
#' @param r \linkS4class{RasterLayer}. Normalized greyscale image. See
#'   \code{\link{normalize}} and \code{\link{gbc}}
#' @param segmentation \linkS4class{RasterLayer}. The result of segmenting
#'   \code{r}. Probably, \code{\link{rings_segmentation}} will be the most used
#'   for fisheye images.
#' @param method Character vector of length one. See details for current
#'   options.
#' @inheritParams thr_image
#' @inheritParams fit_coneshaped_model
#'
#' @return \linkS4class{RasterLayer}.
#'
#' @export
#' @family Tools functions
#'
#' @seealso \code{\link{thr_image}}
#'
#' @references \insertAllCited{}
#'
#' @examples
#' r <- read_caim()
#' blue <- gbc(r$Blue)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' rings <- rings_segmentation(z, 10)
#' bin <- regional_thresholding(blue, rings, "Diaz2018", -8, 0.5, 0.9)
regional_thresholding <- function(r,
                                  segmentation,
                                  method,
                                  intercept = NULL,
                                  slope = NULL,
                                  prob = NULL) {
  stopifnot(class(r) == "RasterLayer")
  stopifnot(class(segmentation) == "RasterLayer")
  stopifnot(class(method) == "character")
  stopifnot(length(method) == 1)

  .check_if_r_was_normalized(r)

  if (method == "Diaz2018") {
    if (any(is.null(intercept), is.null(slope), is.null(prob))) {
      stop("Arguments \"intercept\", \"slope\", and \"prob\" should be provided.")
    }
  }

  fun <- switch(method,
    Diaz2018 = function(dns) {
      dn <- quantile(dns, prob)
      thr_image(dn, intercept, slope)
    }
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
  Map(.binarize_per_ring, raster::unique(segmentation))
  bin
}
