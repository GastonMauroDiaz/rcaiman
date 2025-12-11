#' Write and read out-of-the-box CIE sky model and raster
#'
#' @description
#' Persist and restore the out-of-the-box CIE sky model, its diagnostics, and
#' related rasters/points. The writer produces a human-readable `.txt` manifest
#' plus CSV and GeoPackage sidecar files. The reader reconstructs a list object
#' compatible with the out-of-the-box pipeline.
#'
#' @section Functions:
#' \describe{
#'   \item{`write_sky_cie`}{No return value. Writes four files to disk with the
#'   prefix `name` (see below).}
#'   \item{`read_sky_cie`}{Returns a `list` similar to the output of
#'   `ootb_sky_cie()` and suitable as input to `ootb_sky_above()`.}
#' }
#'
#' @section Files written by `write_sky_cie`:
#' \itemize{
#'   \item Plain text manifest: `name.txt`
#'   \item CSV with sky radiance samples: `name_rr.csv`
#'   \item GeoPackage with sky sample points: `name_sky_points.gpkg`
#'   \item GeoPackage with the sun disk location: `name_sun_disk.gpkg`
#' }
#'
#' @section Text file keys:
#'
#' \describe{
#'   \item{`format_version:`}{Semantic version of the manifest.}
#'   \item{`sun_theta:`}{Solar zenith (deg).}
#'   \item{`sun_phi:`}{Solar azimuth (deg).}
#'   \item{`method_sun:`}{Method used to optimize sun coordinates.}
#'   \item{`zenith_dn:`}{Reference DN at zenith.}
#'   \item{`start_a:`…`start_e:`}{Initial CIE coefficients.}
#'   \item{`fit_a:`…`fit_e:`}{Fitted CIE coefficients.}
#'   \item{`method:`}{Method used to fit CIE coefficients.}
#' }
#'
#' @details
#' Encoding is UTF-8. Decimal point is `.`. Unknown keys are
#' ignored with a warning. Missing required keys trigger an error. The manifest
#' begins with `format_version:` which is checked for basic compatibility.
#'
#' When `read_sky_cie()` detects manual edits (moved sun disk or changed sky
#' points) and `refit_allowed = TRUE`, it re-fits the CIE model using the
#' current `r`, `z`, and `a`, then revalidates.
#'
#' @param name character vector of length one. File base name without extension.
#'   A path can be included, e.g., `"C:/Users/Doe/Documents/DSCN4500"`.
#' @param sky_cie list. Object holding the fitted CIE model, diagnostics, and
#'   derived rasters, as produced by the out-of-the-box workflow.
#' @param r numeric [terra::SpatRaster-class] of one layer. The canopy image
#'   used in the out-of-the-box workflow (used by `read_sky_cie()` when
#'   refitting).
#' @param refit_allowed logical vector of length one. If `TRUE`, allow automatic
#'   re-fit when manual edits are detected.
#'
#' @inheritParams skygrid_segmentation
#'
#' @return See *Functions*.
#'
#' @name write_sky_cie
#' @rdname write_sky_cie
#' @aliases read_sky_cie
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' # Read a previously written model
#' path <- system.file("external/example.txt", package = "rcaiman")
#' sky_cie <- read_sky_cie(gsub(".txt", "", path), r = caim$Blue, z = z, a = a,
#'                         refit_allowed = TRUE)
#' }
write_sky_cie <- function(sky_cie, name) {
  if (!is.list(sky_cie)) {
    stop("`sky_cie` must be a list returned by `ootb_sky_cie()`. ")
  }
  required_names <- c(
                      "rr",
                      "model"
                      )
  if (!all(required_names %in% names(sky_cie))) {
    stop(sprintf("`sky_cie` must contain %s.",
                 paste(sprintf('"%s"', required_names), collapse = ", ")))
  }
  .check_vector(name, "character", 1)

  # helpers ---------------------------------------------------------------
  .print_line <- function(...) cat(paste0(...), "\n")
  .hr <- function(ch = "-", n = 80) paste0(rep(ch, n), collapse = "")

  # manifest --------------------------------------------------------------
  con <- file(paste0(name, ".txt"), open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  .print_line(.hr("-"))
  .print_line("format_version:", " 1.3")
  .print_line("generated_by:", " rcaiman::write_sky_cie()  # do not edit by hand")
  .print_line(.hr("-"))

  # sun and methods
  .print_line("sun_theta:", unname(sky_cie$model$sun_angles["z"]))
  .print_line("sun_phi:",   unname(sky_cie$model$sun_angles["a"]))
  if (!is.null(sky_cie$model$method_sun))
    .print_line("method_sun:", sky_cie$model$method_sun)
  for (m in sky_cie$tested_methods) .print_line("[Tested:", m, "]")


  .print_line(.hr("."))
  # starts
  .print_line("zenith_dn:", sky_cie$model$rr$zenith_dn)
  starts <- unname(sky_cie$model$start)
  if (length(starts) != 5) starts <- rep(NULL, 5)
  .print_line("start_a:", starts[1]); .print_line("start_b:", starts[2])
  .print_line("start_c:", starts[3]); .print_line("start_d:", starts[4])
  .print_line("start_e:", starts[5])
  .print_line(.hr("."))
  # fit
  co <- unname(sky_cie$model$coef[1:5])
  .print_line("fit_a:", co[1]); .print_line("fit_b:", co[2])
  .print_line("fit_c:", co[3]); .print_line("fit_d:", co[4])
  .print_line("fit_e:", co[5])
  .print_line("method:", unname(sky_cie$model$method))

  .print_line(.hr("-"))
  sink()

  # sidecar files ---------------------------------------------------------
  # sky points to GPKG
  cells <- terra::cellFromRowCol(sky_cie$rr,
                                 sky_cie$sky_points$row,
                                 sky_cie$sky_points$col)
  xy <- terra::xyFromCell(sky_cie$rr, cells)
  p <- terra::vect(xy, "points")
  terra::crs(p) <- terra::crs(sky_cie$rr)
  terra::writeVector(p, paste0(name, "_sky_points.gpkg"), filetype = "GPKG")

  # sun disk to GPKG
  cells <- terra::cellFromRowCol(sky_cie$rr,
                                 sky_cie$sun_row_col$row,
                                 sky_cie$sun_row_col$col)
  xy <- terra::xyFromCell(sky_cie$rr, cells)
  p <- terra::vect(xy, "points")
  terra::crs(p) <- terra::crs(sky_cie$rr)
  terra::writeVector(p, paste0(name, "_sun_disk.gpkg"), filetype = "GPKG")

  # rr CSV
  if (!is.null(sky_cie$model$rr$sky_points)) {
    utils::write.csv2(sky_cie$model$rr$sky_points,
                      paste0(name, "_rr.csv"),
                      row.names = FALSE)
  }

  invisible(NULL)
}

#' @rdname write_sky_cie
#' @export
read_sky_cie <- function(name, r, z, a, refit_allowed = FALSE) {
  .check_vector(name, "character", 1)
  .assert_file_exists(paste0(name, ".txt"))
  .check_r_z_a_m(r, z, a, r_type = "single")
  .check_vector(refit_allowed, "logical", 1)

  # helpers ---------------------------------------------------------------
  .read_manifest <- function(path) {
    if (!file.exists(path)) stop("Manifest not found: ", path)
    readLines(path, encoding = "UTF-8", warn = FALSE)
  }
  .get_token_after <- function(lines, key, cast = c("char","num","log"), required = TRUE) {
    cast <- match.arg(cast)
    pat <- paste0("^", key, "\\s*:\\s*(.*)$")
    hit <- grep(pat, lines)
    if (!length(hit)) {
      if (required) stop("Missing key: ", key) else return(NA)
    }
    val <- sub(pat, "\\1", lines[hit[1]])
    if (cast == "num") return(suppressWarnings(as.numeric(val)))
    if (cast == "log") return(tolower(trimws(val)) %in% c("true","t","1"))
    trimws(val)
  }

  # init structure --------------------------------------------------------
  sky_cie <- list()
  sky_cie$io <- list()
  sky_cie$model <- list()
  sky_cie$model$rr <- list()

  # parse manifest --------------------------------------------------------
  lines <- .read_manifest(paste0(name, ".txt"))

  sky_cie$io$format_version <- .get_token_after(lines, "format_version", "char", required = TRUE)
  if (!grepl("^1\\.", sky_cie$io$format_version)) {
    warning("Unknown manifest format_version: ", sky_cie$io$format_version)
  }

  sky_cie$model$sun_angles <- c(
    z = .get_token_after(lines, "sun_theta", "num", TRUE),
    a = .get_token_after(lines, "sun_phi",   "num", TRUE)
  )
  msun <- .get_token_after(lines, "method_sun", "char", required = FALSE)
  if (!is.na(msun)) sky_cie$model$method_sun <- msun

  sky_cie$model$rr$zenith_dn <- .get_token_after(lines, "zenith_dn", "num", TRUE)

  sky_cie$optimal_start <- stats::setNames(
    c(
      .get_token_after(lines, "start_a", "num", TRUE),
      .get_token_after(lines, "start_b", "num", TRUE),
      .get_token_after(lines, "start_c", "num", TRUE),
      .get_token_after(lines, "start_d", "num", TRUE),
      .get_token_after(lines, "start_e", "num", TRUE)
    ),
    c("a","b","c","d","e")
  )

  sky_cie$model$coef <- c(
    .get_token_after(lines, "fit_a", "num", TRUE),
    .get_token_after(lines, "fit_b", "num", TRUE),
    .get_token_after(lines, "fit_c", "num", TRUE),
    .get_token_after(lines, "fit_d", "num", TRUE),
    .get_token_after(lines, "fit_e", "num", TRUE)
  )
  sky_cie$model$method <- .get_token_after(lines, "method", "char", TRUE)

  # sidecar files ---------------------------------------------------------

  # rr
  sky_cie$model$rr$sky_points <- utils::read.csv2(paste0(name, "_rr.csv"))

  # sky_points from GPKG → row/col
  sp <- terra::vect(paste0(name, "_sky_points.gpkg"))
  ext <- terra::extract(z, sp, cells = TRUE)
  rc  <- terra::rowColFromCell(z, ext$cell)
  sky_cie$sky_points <- data.frame(row = rc[,1], col = rc[,2])

  # sun disk from GPKG → row/col and angles
  sd <- terra::vect(paste0(name, "_sun_disk.gpkg"))
  ext <- terra::extract(z, sd, cells = TRUE)
  rc  <- terra::rowColFromCell(z, ext$cell)
  sun_row_col <- data.frame(row = rc[,1], col = rc[,2])

  # consistency checks and optional refit --------------------------------
  # angular distance between stored sun and GPKG sun
  if (nrow(sun_row_col) != 1) {
    warning("Sun is outside the raster grid; cannot verify manual edit of sun position.")
    delta <- 0
  } else {
    sza_saa <- zenith_azimuth_from_row_col(z, a, sun_row_col$row, sun_row_col$col)
    sza_saa <- as.numeric(sza_saa)
    delta <- calc_spherical_distance(
      .degree2radian(sky_cie$model$sun_angles["z"]),
      .degree2radian(sky_cie$model$sun_angles["a"]),
      .degree2radian(sza_saa[1]),
      .degree2radian(sza_saa[2])
    )
    names(sza_saa) <- c("z","a")
  }

  .norm_points <- function(x) {
    x <- x[, c("row","col"), drop = FALSE]
    x$row <- as.integer(x$row)
    x$col <- as.integer(x$col)
    x <- x[order(x$row, x$col), , drop = FALSE]
    rownames(x) <- NULL
    x
  }

  points_changed <- !identical(
    .norm_points(sky_cie$model$rr$sky_points),
    .norm_points(sky_cie$sky_points)
  )

  if ( (delta > pi/180) || points_changed ) {
    msg <- "Detected manual changes. "
    if (isTRUE(refit_allowed)) {
      rr <- extract_rr(r, z, a, sky_cie$sky_points, use_window = TRUE)
      sky_cie$model <- fit_cie_model(rr, sza_saa,
                                     custom_sky_coef = sky_cie$model$coef)
      sky_cie$model_validation <- validate_cie_model(sky_cie$model, k = 10)
      message(msg, "Automatically re-fitting the CIE model.")
    } else {
      warning(msg, "Re-fitting the CIE model is recommended.")
    }
  }

  # regenerate rr for convenience
  sky_cie$rr <- cie_image(z, a, sky_cie$model$sun_angles, sky_cie$model$coef)

  sky_cie$io$source <- name
  sky_cie
}
