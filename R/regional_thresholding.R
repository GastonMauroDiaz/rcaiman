#' Regional thresholding
#'
#' Regional thresholding of greyscale images
#'
#' Methods currently implemented are:
#'
#' \itemize{ \item \strong{Diaz2018}: method presented in
#' \insertCite{Diaz2018}{rcaiman} applied regionally. If this method is
#' selected, the arguments \code{w}, \code{type},and \code{prob} should be
#' provided. It works extracting the digital numbers from each segment and
#' passing them to \code{quantile(x, prob)}, which result is in turn passed to
#' \code{thr_image(dns, w, type)}.
#'
#' \item \strong{Methods from autothresholdr package}: this function can call
#' methods from \code{\link[autothresholdr]{auto_thresh}}. Use \code{"IsoData"}
#' to use the algorithm by \insertCite{isodata}{rcaiman}, which is the one
#' recommended by \insertCite{Jonckheere2005}{rcaiman}. }
#'
#' @param r \linkS4class{RasterLayer}. Normalized greyscale image. See
#'   \code{\link{normalize}} and \code{\link{gbc}}
#' @param segmentation \linkS4class{RasterLayer}. The result of segmenting
#'   \code{r}. Probably, code{\link{rings_segmentation}} will be the most used
#'   for fisheye images.
#' @param method Character vector of length one. See details for current
#'   options.
#' @inheritParams thr_image
#' @inheritParams model_sky_dn
#'
#'
#'
#' @return \linkS4class{RasterLayer}.
#' @export
#'
#' @seealso \code{\link{thr_image}}
#'
#' @references \insertRef{Diaz2018}{rcaiman}
#'
#'   \insertCite{Jonckheere2005}{rcaiman}
#'
#'   \insertCite{isodata}{rcaiman}
#'
#' @examples
#' r <- read_caim()
#' blue <- gbc(r$Blue)
#' z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
#' rings <- rings_segmentation(z, 10)
#' bin <- regional_thresholding(blue, rings, "Diaz2018", 0.5, "Generic", 0.9)
regional_thresholding <- function(r, segmentation, method,
                                  w = NULL, type = NULL, prob = NULL) {
  stopifnot(class(r) == "RasterLayer")
  stopifnot(class(segmentation) == "RasterLayer")
  stopifnot(class(method) == "character")
  stopifnot(length(method) == 1)

  if (method == "Diaz2018") {
    if (any(is.null(w), is.null(type), is.null(prob)))
      stop("Arguments \"w\", \"type\",and \"prob\" should be provided.")

  }

  fun <- switch(method,
    Diaz2018 = function(dns) {
      dn <- quantile(dns, prob)
      thr_image(dn, w, type)
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
