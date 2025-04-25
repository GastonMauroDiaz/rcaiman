#' Regional thresholding
#'
#' Regional thresholding of greyscale images.
#'
#' Methods currently implemented are:
#'
#' * __Methods from autothresholdr package__: this function can call
#' methods from [autothresholdr::auto_thresh()]. For instance, use `"IsoData"`
#' to use the algorithm by \insertCite{isodata;textual}{rcaiman}, which was
#' recommended by \insertCite{Jonckheere2005;textual}{rcaiman}.
#' * __Method isodata from this package__: Use `"thr_isodata"` to
#' use [thr_isodata()].
#'
#' @inheritParams obia
#' @inheritParams sky_grid_segmentation
#' @param segmentation [SpatRaster-class]. The result of segmenting `r`.
#'   Arguably, the result of calling [rings_segmentation()] will be the
#'   preferred choice for fisheye images.
#' @param method Character vector of length one. See details for current
#'   options.
#' @inheritParams thr_mblt
#' @inheritParams fit_trend_surface
#'
#' @return An object of class [SpatRaster-class] with values `0` and `1`.
#'
#' @export
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' zenith_colrow <- c(1276, 980)
#' diameter <- 756*2
#' caim <- read_caim(path, zenith_colrow - diameter/2, diameter, diameter)
#' z <- zenith_image(ncol(caim), lens("Nikon_FCE9"))
#' r <- gbc(caim$Blue)
#' r <- correct_vignetting(r, z, c(0.0638, -0.101)) %>% normalize_minmax()
#' rings <- rings_segmentation(z, 15)
#' bin <- regional_thresholding(r, rings, "thr_isodata")
#' plot(bin)
#' }
regional_thresholding <- function(r, segmentation, method) {
  .is_single_layer_raster(r, "r")
  .is_single_layer_raster(segmentation, "segmentation")
  stopifnot(class(method) == "character")
  stopifnot(length(method) == 1)

  fun <- switch(method,
    thr_isodata = thr_isodata
  )

  if (is.null(fun)) {
    .was_normalized(r, "r")
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
  segs <- segs[segs != 0]
  Map(.binarize_per_ring, segs)
  bin
}
