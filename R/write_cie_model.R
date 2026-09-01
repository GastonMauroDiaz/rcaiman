#' Write and read CIE sky models
#'
#' @description
#' Persist and restore the a CIE sky model. The writer produces a human-readable
#' `.txt` model record plus CSV and GeoPackage sidecar files. The reader
#' reconstructs a list object compatible with the \pkg{rcaiman} pipeline.
#'
#' @section Functions:
#' \describe{
#'   \item{`write_cie_model`}{No return value. Writes four files to disk with the
#'   prefix `name` (see below).}
#'   \item{`read_cie_model`}{Returns a `list` similar to the output of
#'   [fit_cie_model()] or [optimize_sun_angles()].}
#' }
#'
#' @section Files written by `write_cie_model`:
#' \itemize{
#'   \item Plain text model record: `<name>_model_record.txt`
#'   \item CSV with sky radiance samples: `<name>_rr.csv`
#'   \item GeoPackage with sky sample points: `<name>_sky_points.gpkg`
#'   \item GeoPackage with the sun disk location: `<name>_sun_disk.gpkg`
#' }
#'
#' @section Text file keys:
#'
#' \describe{
#'   \item{`format_version:`}{Semantic version of the model record.}
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
#' Encoding is UTF-8. Decimal point is `.`. Unknown keys are ignored with a
#' warning. Missing required keys trigger an error. The model record begins with
#' `format_version:` which is checked for basic compatibility. It
#' reflects the on-disk model record format and differs from package
#' versioning.
#'
#' When `read_cie_model()` detects manual edits (moved sun disk or changed sky
#' points) and `refit_allowed = TRUE`, it re-fits the CIE model using the
#' current `r`, `z`, and `a`. The model record remains unchanged.
#'
#' @param name character vector of length one. File base name without extension.
#'   A path can be included, e.g., `"C:/Users/Doe/Documents/DSCN4500"`.
#' @param r numeric [terra::SpatRaster-class] of one layer. The canopy image
#'   used to obtain sky relative radiance (see [fit_cie_model()]).
#' @param refit_allowed logical vector of length one. If `TRUE`, allow automatic
#'   re-fit when manual edits are detected.
#'
#' @inheritParams optimize_sun_angles
#' @inheritParams add_directions
#'
#' @return See *Functions*.
#'
#' @name write_cie_model
#' @rdname write_cie_model
#' @aliases read_cie_model
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' plot(r)
#'
#' # Manual method following Lang et al. (2010)
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#' plot(r)
#' points(sky_points$col, nrow(r) - sky_points$r, col = 2, pch = 10)
#'
#' # x11()
#' # plot(r)
#' # sun_angles <- click(c(z, a), 1) %>% as.numeric()
#' sun_angles <- c(zen_deg = 49.5, azi_deg = 27.4) #taken with above lines then hardcoded
#'
#' sun_row_col <- row_col_from_zen_azi(z, a,
#'                                     sun_angles["zen_deg"] * pi/180,
#'                                     sun_angles["azi_deg"] * pi/180)
#' points(sun_row_col[2], nrow(r) - sun_row_col[1],
#'        col = "yellow", pch = 8, cex = 3)
#'
#' sky_points <- extract_dn(r, sky_points)
#' sky_points <- add_directions(sky_points, z, a)
#' rr <- calc_relative_radiance(sky_points, 3)
#'
#' set.seed(7)
#' model <- fit_cie_model(rr, sun_angles,
#'                        general_sky_type = "Clear")
#'
#' write_cie_model(model, z, a, "example")
#' model2 <- read_cie_model("example", r, z, a)
#' }
write_cie_model <- function(model, z, a, name) {
  .assert_cie_model(model)
  .assert_zenith_image(z)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .check_vector(name, "character", 1)

  # helpers ---------------------------------------------------------------
  .print_line <- function(...) cat(paste0(...), "\n")
  .hr <- function(ch = "-", n = 80) paste0(rep(ch, n), collapse = "")

  # modelrecord -----------------------------------------------------------
  con <- file(paste0(name, "_model_record.txt"), open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  .print_line(.hr("-"))
  .print_line("format_version:", " 1.4")
  .print_line("generated_by:", " rcaiman::write_cie_model()  # do not edit by hand")
  .print_line(.hr("-"))

  # sun and methods
  .print_line("sun_theta:", unname(model$sun_angles["zen_deg"]))
  .print_line("sun_phi:",   unname(model$sun_angles["azi_deg"]))
  if (!is.null(model$method_sun))
    .print_line("method_sun:", model$method_sun)

  .print_line(.hr("."))
  # starts
  .print_line("zenith_dn:", model$rr$zenith_dn)
  starts <- unname(model$start)
  if (length(starts) != 5) starts <- rep(NULL, 5)
  .print_line("start_a:", starts[1]); .print_line("start_b:", starts[2])
  .print_line("start_c:", starts[3]); .print_line("start_d:", starts[4])
  .print_line("start_e:", starts[5])
  .print_line(.hr("."))
  # fit
  co <- unname(model$coef[1:5])
  .print_line("fit_a:", co[1]); .print_line("fit_b:", co[2])
  .print_line("fit_c:", co[3]); .print_line("fit_d:", co[4])
  .print_line("fit_e:", co[5])
  .print_line("method:", unname(model$method))

  .print_line(.hr("-"))
  sink()

  # sidecar files ---------------------------------------------------------
  # sky points to GPKG
  cells <- terra::cellFromRowCol(z,
                                 model$rr$sky_points$row,
                                 model$rr$sky_points$col)
  xy <- terra::xyFromCell(z, cells)
  p <- terra::vect(xy, "points")
  terra::crs(p) <- terra::crs(z)
  terra::writeVector(p, paste0(name, "_sky_points.gpkg"), filetype = "GPKG")

  # sun disk to GPKG
  sun_row_col <- row_col_from_zen_azi(z, a,
                                      model$sun_angles["zen_deg"] * pi/180,
                                      model$sun_angles["azi_deg"] * pi/180)

  cells <- terra::cellFromRowCol(z,
                                 sun_row_col$row,
                                 sun_row_col$col)
  xy <- terra::xyFromCell(z, cells)
  p <- terra::vect(xy, "points")
  terra::crs(p) <- terra::crs(z)
  terra::writeVector(p, paste0(name, "_sun_disk.gpkg"), filetype = "GPKG")

  # rr CSV
  if (!is.null(model$rr$sky_points)) {
    utils::write.csv2(model$rr$sky_points,
                      paste0(name, "_rr.csv"),
                      row.names = FALSE)
  }

  invisible(NULL)
}

#' @rdname write_cie_model
#' @export
read_cie_model <- function(name, r, z, a, refit_allowed = FALSE) {
  .check_vector(name, "character", 1)
  .assert_file_exists(paste0(name, "_model_record.txt"))
  .assert_single_layer(r, allow_na = FALSE)
  .assert_zenith_image(z)
  .assert_same_geom(z, r)
  .assert_azimuth_image(a)
  .assert_same_geom(z, a)
  .check_vector(refit_allowed, "logical", 1)

  # helpers ---------------------------------------------------------------
  .read_model_record <- function(path) {
    if (!file.exists(path)) stop("Model record not found: ", path)
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
  model <- list()
  model$io <- list()
  model$rr <- list()

  # parse model record ----------------------------------------------------
  lines <- .read_model_record(paste0(name, "_model_record.txt"))

  model$io$format_version <- .get_token_after(lines, "format_version", "char", required = TRUE)
  if (!grepl("^1\\.", model$io$format_version)) {
    warning("Unknown format_version: ", model$io$format_version)
  }

  model$sun_angles <- c(
    zen_deg = .get_token_after(lines, "sun_theta", "num", TRUE),
    azi_deg = .get_token_after(lines, "sun_phi",   "num", TRUE)
  )
  msun <- .get_token_after(lines, "method_sun", "char", required = FALSE)
  if (!is.na(msun)) model$method_sun <- msun

  model$rr$zenith_dn <- .get_token_after(lines, "zenith_dn", "num", TRUE)

  model$coef <- stats::setNames(
    c(
      .get_token_after(lines, "fit_a", "num", TRUE),
      .get_token_after(lines, "fit_b", "num", TRUE),
      .get_token_after(lines, "fit_c", "num", TRUE),
      .get_token_after(lines, "fit_d", "num", TRUE),
      .get_token_after(lines, "fit_e", "num", TRUE)
    ),
    c("a","b","c","d","e")
  )

  model$start <- c(
    .get_token_after(lines, "start_a", "num", TRUE),
    .get_token_after(lines, "start_b", "num", TRUE),
    .get_token_after(lines, "start_c", "num", TRUE),
    .get_token_after(lines, "start_d", "num", TRUE),
    .get_token_after(lines, "start_e", "num", TRUE)
  )
  model$method <- .get_token_after(lines, "method", "char", TRUE)

  # sidecar files ---------------------------------------------------------

  # rr
  try(model$rr$sky_points <- utils::read.csv2(paste0(name, "_rr.csv")))

  # sky_points from GPKG → row/col
  sp <- terra::vect(paste0(name, "_sky_points.gpkg"))
  ext <- terra::extract(z, sp, cells = TRUE)
  rc  <- terra::rowColFromCell(z, ext$cell)
  sky_points <- data.frame(row = rc[,1], col = rc[,2])

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
    sza_saa <- zen_azi_from_row_col(z, a, sun_row_col$row, sun_row_col$col)
    sza_saa <- as.numeric(sza_saa)
    delta <- calc_spherical_distance(
      .degree2radian(model$sun_angles["zen_deg"]),
      .degree2radian(model$sun_angles["azi_deg"]),
      sza_saa[1],
      sza_saa[2]
    )
    sza_saa <- .radian2degree(sza_saa)
    names(sza_saa) <- c("zen_deg","azi_deg")
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
    .norm_points(model$rr$sky_points),
    .norm_points(sky_points)
  )

  if ( (delta > pi/180) || points_changed ) {
    msg <- "Detected manual changes. "
    if (isTRUE(refit_allowed)) {
      sky_points <- extract_dn(r, sky_points, use_window = TRUE)
      sky_points <- add_directions(sky_points, z, a)
      rr <- calc_relative_radiance(sky_points, names(r), zenith_dn = model$rr$zenith_dn)
      model <- fit_cie_model(rr, sza_saa, custom_sky_coef = model$start)
      message(msg, "Automatically re-fitting the CIE model.")
    } else {
      warning(msg, "Re-fitting the CIE model is recommended.")
    }
  }

  model$io$source <- name
  model
}
