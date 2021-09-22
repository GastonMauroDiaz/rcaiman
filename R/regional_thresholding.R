#' Regional thresholding
#'
#' Regional thresholding of greyscale images
#'
#' Methods currently implemented are:
#'
#' @section mblt: is the method described in \insertCite{Diaz2018}{rcaiman} but
#'   applied regionally instead of locally. Digital numbers from each segment
#'   are extracted and passed to \code{quantile(x, 0.9)}, which result is in
#'   turn passed to \code{\link{thr_fun}}. The "Generic" coefficients are used.
#'   Use mblt_default to use \code{w} equal to \code{0.5}, or mblt_initial_value
#'   to use \code{w} equal to \code{1}.
#'
#' @section methods from autothresholdr package: this function can call methods
#'   from \code{\link[autothresholdr]{auto_thresh}}. Use "IsoData" to use the
#'   algorithm by \insertCite{isodata}{rcaiman}, which is the one recommended by
#'   \insertCite{Jonckheere2005}{rcaiman}.
#'
#' @param r \linkS4class{RasterLayer}. Normalized greyscale image. See
#'   \code{\link{normalize}} and \code{\link{gbc}}
#' @param segmentation \linkS4class{RasterLayer}. The result of segmenting
#'   \code{r}. Probably, code{\link{rings_segmentation}} will be the most
#'   used for fisheye images.
#' @param method Character vector of length one. See details for current
#'   options.
#'
#' @return \linkS4class{RasterLayer}.
#' @export
#'
#' @references
#' \insertRef{Diaz2018}{rcaiman}
#'
#' \insertCite{Jonckheere2005}{rcaiman}
#'
#' \insertCite{isodata}{rcaiman}
#'
#' @examples
#' r <- read_caim()
#' blue <- gbc(r$Blue)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' rings <- rings_segmentation(z, 10)
#' bin <- regional_thresholding(blue, rings, "IsoData")
regional_thresholding <- function(r, segmentation, method) {
  stopifnot(class(r) == "RasterLayer")
  stopifnot(class(segmentation) == "RasterLayer")
  stopifnot(class(method) == "character")
  stopifnot(length(method) == 1)

  if (method == "IsoData") {
    if (!requireNamespace("autothresholdr", quietly = TRUE)) {
      stop(paste(
        "Package \"autothresholdr\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
      )
    }
  }

  fun <- switch(method,
    mblt_default = function(dns) {
      dn <- quantile(dns, 0.9)
      thr_fun(dn * 255, 0.5) / 255
    },
    mblt_initial_value = function(dns) {
      dn <- quantile(dns, 0.9)
      thr_fun(dn * 255, 1) / 255
    }
  )

  if (is.null(fun)) {
    fun <- function(dns) {
      autothresholdr::auto_thresh(round(dns * 255), method)[1] / 255
    }
  }

  bin <- apply_thr(r, .get_min(r))
  .binarize_per_ring <- function(segment_id) {
    indices <- segmentation == segment_id
    thr <- fun(r[indices])
    bin[indices] <<- r[indices] > thr
  }
  Map(.binarize_per_ring, unique(segmentation))
  bin
}
