#' Regional thresholding of greyscale images
#'
#' Perform thresholding of greyscale images by applying a method
#' regionally, using a segmentation map.
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
#'   [thr_isodata()], a native implementation of the same algorithm}
#'   \item{Two-corner method:}{Use `"thr_two_corner"` to apply
#'   [thr_two_corner()] with `sigma = NULL`, `slope_reduction = TRUE`,
#'   and `method = macfarlane`, which implements a geometric thresholding
#'   strategy based on identifying inflection points in the histogram, first
#'   introduced to canopy photography by
#'   \insertCite{Macfarlane2011;textual}{rcaiman}. Since this method might fail,
#'   the fallback is `thr_isodata`}
#'   \item{Nobis-Hunziker:}{Use `"thr_nbhu"` to apply
#'   [thr_nbhu()] with `support_weight = 0`, which implements the edge-based criterion proposed by
#'   \insertCite{Nobis2005;textual}{rcaiman}.}
#' }
#'
#' @param r numeric [terra::SpatRaster-class] of one layer. Typically the blue
#'   channel of a canopy photograph.
#' @param seg numeric [terra::SpatRaster-class] of one layer. A labeled
#'   segmentation map defining the regions over which to apply the thresholding
#'   method. Ring segmentation (see [ring_segmentation()]) is often preferred
#'   for fisheye images \insertCite{Leblanc2005}{rcaiman}.
#' @param method character vector of length one. Name of the thresholding method
#'   to apply. See *Details*.
#'
#' @note
#' When methods from the `autothresholdr` package are used, `r` values
#' should be constrained to the range \eqn{[0, 1]}. See [normalize_min_max()].
#'
#' @seealso [apply_by_direction()] to binarize by direction
#'
#' @return Logical [terra::SpatRaster-class] (`TRUE` for sky, `FALSE` for
#'   non-sky) of the same dimensions as `r`.
#'
#' @export
#'
#' @references \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_jpeg.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry, "exif_jpeg", "delta_t")
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' caim <- read_caim(path, zenith_col_row = spec$zenith_col_row ,
#'                   horizon_radius = spec$horizon_radius)
#' z <- zenith_image(ncol(caim), spec$parameters)
#' rspec <- get_radiometry_spec(spec, "wagner2001_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`0`)
#' r <- apply_flat_field(caim$Blue, ff)
#'
#' rings <- ring_segmentation(z, 15)
#' bin <- binarize_by_region(r, rings, "thr_two_corner")
#' plot(bin, main = "Twocorner")
#'
#' a <- azimuth_image(z)
#' seg <- equal_area_segmentation(z, a, 20)
#' bin <- binarize_by_region(r, seg, "thr_isodata")
#' plot(bin, main = "Isodata")
#'
#' bin <- binarize_by_region(r, seg, "thr_nbhu")
#' plot(bin, main = "Nobis and Hunziker")
#'
#' # Method from 'autothresholdr' package:
#' bin <- binarize_by_region(r, seg, "Huang2")
#' plot(bin, main = "Huang2")
#' }
binarize_by_region <- function(r, seg, method) {
  .assert_single_layer(r, allow_na = FALSE)
  .assert_single_layer(seg)
  .assert_same_geom(r, seg)
  .check_vector(method, "character", 1)

  .fun <- switch(method,
    thr_isodata = function(x, ...) thr_isodata(x),
    thr_two_corner = function(x, ...) {
      thr <- tryCatch(thr_two_corner(x,
                                     sigma = NULL,
                                     slope_reduction = TRUE,
                                     method = "macfarlane")$tm, error = function(e) NULL)
      if (is.null(thr)) {
        message("Fallback to 'thr_isodata()'")
        thr <- thr_isodata(x)
      }
      thr
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
  bin <- binarize_with_thr(r, .get_min(r))
  .binarize_per_ring <- function(segment_id) {
    indices <- seg == segment_id
    thr <- .fun(r[indices], r, indices)
    bin[indices] <<- r[indices] > thr
  }
  segs <- unique(terra::values(seg)) %>% as.numeric()
  segs <- segs[!is.na(segs)]
  segs <- segs[segs != 0]
  Map(.binarize_per_ring, segs)
  bin
}
