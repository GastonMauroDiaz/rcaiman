#' Regional thresholding of greyscale images
#'
#' Perform thresholding of greyscale images by applying a method regionally,
#' using a segmentation map. See [apply_by_direction()] for a cutting-edge
#' alternative.
#'
#' @details
#' This function supports several thresholding methods applied within the
#' regions defined by `seg`:
#'
#' \describe{
#'   \item{Methods from the `autothresholdr` package:}{Any method supported by
#'   [autothresholdr::auto_thresh()] can be used by specifying its name.
#'   For example, `"IsoData"` applies the classic iterative intermeans algorithm
#'   \insertCite{isodata;textual}{rcaiman}, which is among the most recommended
#'   for canopy photography \insertCite{Jonckheere2005}{rcaiman}.}
#'    \item{In-package implementation of IsoData:}{Use `"thr_isodata"` to apply
#'   [thr_isodata()], with optional observation weights.}
#'   \item{Two-corner method:}{Use `"thr_two_corner"` to apply
#'   [build_histogram()] with 256 bins, [find_macfarlane_peaks()], and
#'   [thr_two_corner()] without smoothing and with `slope_reduction = TRUE`.
#'   This geometric thresholding strategy was first
#'   introduced to canopy photography by
#'   \insertCite{Macfarlane2011;textual}{rcaiman}. After all segments have been
#'   processed, failed estimates are replaced by the arithmetic mean of the
#'   successful thresholds (`tm`). Each successful segment contributes equally
#'   to this mean, regardless of its size or total weight. If every segment
#'   fails, the function stops because no fallback threshold can be calculated.}
#'   \item{Nobis-Hunziker:}{Use `"thr_nbhu"` to apply
#'   [thr_nbhu()] with `support_weight = 0`, which implements the edge-based criterion proposed by
#'   \insertCite{Nobis2005;textual}{rcaiman}.}
#' }
#'
#' @param r numeric [terra::SpatRaster-class] of one layer. Typically the blue
#'   band of a canopy photograph. Must not contain `NA`.
#' @param seg numeric [terra::SpatRaster-class] of one layer. A labeled
#'   segmentation map defining the regions over which to apply the thresholding
#'   method. Labels zero and `NA` are excluded. Ring segmentation (see [ring_segmentation()]) is often preferred
#'   for fisheye images \insertCite{Leblanc2005}{rcaiman}.
#'   A logical mask (see [direction_image()] for examples) is also accepted:
#'   `TRUE` is label 1 and `FALSE` is label 0. All selected pixels belong to
#'   one segment, even when disconnected. Use [label_components()] if connected
#'   regions should have separate labels.
#' @param method character vector of length one. Name of the thresholding method
#'   to apply. See *Details*.
#' @param w optional single-layer raster of nonnegative, finite weights, or
#'   `NA`, on the same grid as `r`. Supported with `"thr_isodata"` and
#'   `"thr_two_corner"`; supplying weights for another method causes an
#'   error. With `NULL`, every pixel has weight one. For solid-angle weights,
#'   use [area_weight_image()].
#'
#' @section Solid-angle weighting:
#' Within each segment, matching DN and weight values are supplied to
#' [thr_isodata()] or [build_histogram()], according to `method`.
#' Zero or missing weights exclude observations from threshold estimation.
#' The resulting threshold is applied to every pixel in that segment, including
#' pixels excluded from estimation. Weights neither modify DN values nor define
#' the region to classify; `seg` defines that region.
#'
#' Equal-area segments do not imply equal solid angles per pixel. Likewise,
#' within a zenith ring in an inclined acquisition, off-axis angles and pixel
#' solid angles can vary. Weights account for these differences in the DN
#' distribution used to estimate each threshold. Only relative weights matter.
#'
#' At least two observations with positive, nonmissing weights are required
#' to estimate a segment's threshold. With `"thr_isodata"`, a failed estimate
#' stops the calculation and the error identifies the segment. If all retained
#' values are equal, that common value becomes the threshold and a warning is
#' issued. With `"thr_two_corner"`, insufficient or constant data instead use
#' the mean-threshold fallback described in *Details*. Pixels equal to their
#' threshold are classified as `FALSE`.
#'
#' @note
#' For methods from `autothresholdr`, the image is normalized to \eqn{[0, 1]}
#' with [normalize_min_max()] and DN values are rounded to 8-bit integers for
#' threshold estimation.
#'
#' @seealso [apply_by_direction()] to binarize by direction
#'
#' @return Logical [terra::SpatRaster-class] without `NA`, on the same grid as
#'   `r`. Within processed regions, `TRUE` represents sky and `FALSE` non-sky.
#'   Cells labeled zero or `NA` in `seg` retain the initial value `FALSE`.
#'   This is a fill value, not a classification of those cells; use a mask to
#'   delimit the analysis domain.
#'
#' @export
#'
#' @references \insertAllCited{}
#'
#' @examples
#' # Synthetic radiance gradient on an inclined radial geometry.
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(80, 60), geometry)
#' d <- direction_image(g, elevation = 32.5)
#' r <- g$off_axis + 20
#' r[is.na(r)] <- 0
#' rings <- ring_segmentation(d, 30)
#' w <- area_weight_image(g, geometry)
#'
#' bin <- binarize_by_region(r, rings, "thr_isodata")
#' bin_weighted <- binarize_by_region(r, rings, "thr_isodata", w = w)
#' plot(c(bin, bin_weighted), main = c("Equal pixel weights", "Solid-angle weights"))
#' sum(bin[] != bin_weighted[])
#'
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/DSCN4248.tif", package = "rcaiman")
#' metadata <- read.delim(gsub("tif", "tsv", path))
#' sensor_data <- read_caim(path)
#'
#' # Get specs
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' raw <- get_signature(my_entry, "raw")
#' cygm <- get_band_spec(raw, "cygm")
#' raw_geometry <- get_geometry_spec(cygm, "simple_method")
#' raw_radiometry <- get_radiometry_spec(raw_geometry, "photometric_sphere")
#'
#' # Represent image geometry
#' g <- set_internal_geometry(raw$dim, raw_geometry)
#' d <- direction_image(g)
#' m <- !is.na(g$off_axis)
#' w <- area_weight_image(g, raw_geometry)
#'
#' # Read and preprocess the image
#' caim <- split_cfa(sensor_data, cfa_pattern = cygm$cfa_pattern)
#' caim
#' r <- caim$Cyan
#' aperture <- as.character(metadata$aperture)
#' ff <- trend_surface_image(r, raw_radiometry$parameters[[aperture]]$Cyan, m = m)
#' r <- apply_flat_field(r, ff)
#'
#' # Binarize by region
#' seg <- beckers_segmentation(d, 20)
#' plot(seg)
#' bin <- binarize_by_region(r, seg, "thr_isodata", w)
#' plot(bin, main = "thr_isodata with w")
#'
#' bin <- binarize_by_region(r, seg, "thr_isodata")
#' plot(bin, main = "thr_isodata")
#'
#' bin <- binarize_by_region(r, seg, "thr_two_corner")
#' plot(bin, main = "thr_two_corner")
#'
#' bin <- binarize_by_region(r, seg, "thr_nbhu")
#' plot(bin, main = "thr_nbhu")
#'
#' # Method from 'autothresholdr' package:
#' bin <- binarize_by_region(r, seg, "Huang2")
#' plot(bin, main = "Huang2")
#' }
binarize_by_region <- function(r, seg, method, w = NULL) {
  .assert_single_layer(r, allow_na = FALSE)
  labels <- .segment_values(seg)
  .assert_same_geom(r, seg)
  .check_vector(method, "character", 1)
  if (!is.null(w)) {
    if (!method %in% c("thr_isodata", "thr_two_corner")) {
      stop("`w` is supported only with 'thr_isodata' and 'thr_two_corner'.", call. = FALSE)
    }
    .assert_area_weight_image(w)
    .assert_same_geom(r, w)
  }

  .fun <- switch(method,
    thr_isodata = function(x, r, indices) {
      thr_isodata(x, w = if (is.null(w)) NULL else w[indices])
    },
    thr_two_corner = function(x, r, indices) {
      weights <- if (is.null(w)) NULL else w[indices]
      h <- build_histogram(x, weights, n_bins = 256)
      peaks <- find_macfarlane_peaks(h)
      thr_two_corner(h, peaks, slope_reduction = TRUE)$tm
    },
    thr_nbhu = function(x, r, indices) {
      r <- terra::trim(r[indices,, drop = FALSE])
      m <- !is.na(r)
      r[!m] <- 0
      thr_nbhu(r,
               max(round(min(r[m])), 1),
               round(max(r[m])),
               support_weight = 0)
    }
  )

  if (is.null(.fun)) {
    r <- normalize_min_max(r)
    if (!requireNamespace("autothresholdr", quietly = TRUE)) {
      stop(paste(
        "Package \"autothresholdr\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
      )
    }

    .fun <- function(dns, ...) {
      autothresholdr::auto_thresh(round(dns * 255),
                                  method,
                                  ignore_black = TRUE,
                                  ignore_white = TRUE,
                                  ignore_na = TRUE)[1] / 255
    }
  }
  bin <- terra::rast(r)
  bin[] <- FALSE
  .threshold_per_region <- function(segment_id) {
    indices <- !is.na(seg) & seg == segment_id
    tryCatch({
      value <- .fun(r[indices], r, indices)
      .check_vector(value, "numeric", 1, name = "threshold")
      value
    }, error = function(e) {
      if (method == "thr_two_corner") return(NA_real_)
      stop("Threshold estimation failed for segment ", segment_id, ": ",
           conditionMessage(e), call. = FALSE)
    })
  }
  segs <- unique(labels)
  segs <- segs[!is.na(segs)]
  segs <- segs[segs != 0]
  thresholds <- vapply(segs, .threshold_per_region, numeric(1))
  failed <- is.na(thresholds)
  if (any(failed)) {
    if (all(failed)) {
      stop("Two-corner threshold estimation failed for every segment; no successful thresholds are available to calculate the mean.",
           call. = FALSE)
    }
    fallback <- mean(thresholds[!failed])
    thresholds[failed] <- fallback
    message("Using mean threshold ", signif(fallback, 6), " for ", sum(failed),
            " segment(s) where two-corner estimation failed.")
  }
  for (i in seq_along(segs)) {
    indices <- !is.na(seg) & seg == segs[i]
    bin[indices] <- r[indices] > thresholds[i]
  }
  bin
}
