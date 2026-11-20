#' HSP compatibility functions
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Read and write legacy files from HSP (HemiSPherical Project Manager) projects
#' to interoperate with existing workflows. Intended for legacy support; not
#' required when working fully within `rcaiman`.
#'
#' @section About HSP software:
#'
#' HSP was introduced in \insertCite{Lang2013}{rcaiman} and is based on the
#' method presented in \insertCite{Lang2010;textual}{rcaiman}. It runs
#' exclusively on Windows. HSP stores pre-processed images as PGM files in the
#' `manipulate` subfolder of each project (itself inside the `projects` folder).
#'
#' @section Functions:
#'
#' \describe{
#'   \item{\code{hsp_read_manual_input()}}{read sky marks and sun position
#'   defined manually within an HSP project; returns a named list with
#'   components \code{weight}, \code{max_points}, \code{angle},
#'   \code{point_radius}, \code{sun_indices}, \code{sky_indices},
#'   \code{zenith_dn}, \code{masks}, and \code{marked_areas}.
#'   The first four components and `zenith_dn` are numeric scalars read from
#'   HSP settings and statistics. `sun_indices` is a one-row `data.frame` with
#'   numeric `row` and `col` columns. `sky_indices` is a `data.frame` with
#'   numeric `col`, `row`, and `type` columns; `type` retains the HSP mark code.
#'   Both tables use cell indices numbered from 1 and can be passed to
#'   [add_directions()]. `masks` and `marked_areas` are polygon
#'   [terra::SpatVector-class] objects on the image grid from [read_caim()],
#'   or `NULL` when the corresponding HSP file contains no polygons.}
#'
#'   \item{\code{hsp_read_opt_sky_coef()}}{read optimized CIE sky coefficients
#'   from an HSP project; returns a numeric vector of length five, in CIE
#'   coefficient order `a`, `b`, `c`, `d`, `e`.}
#'
#'   \item{\code{hsp_write_sky_indices()}}{write sky cell indices to
#'   `<img_name>_points.conf` in `manipulate`, replacing an existing file.
#'   Each mark is written with HSP type code 3. Returns `NULL` invisibly.}
#'
#'   \item{\code{hsp_write_sun_indices()}}{write solar disk cell indices to
#'   `<img_name>_sun.conf` in `manipulate`, replacing an existing file.
#'   Returns `NULL` invisibly.}
#' }
#'
#' @param sky_indices `data.frame` with numeric columns `"col"` and `"row"`
#'   containing positive, integer-valued cell indices numbered from 1, without
#'   missing values. Additional columns are ignored, including `type`.
#' @param path_to_HSP_project character vector of length one. Path to the HSP
#'   project folder (e.g., `"C:/Users/johndoe/Documents/HSP/projects/my_prj/"`).
#' @param img_name character vector of length one (e.g., `"DSCN6342.pgm"` or
#'   `"DSCN6342"`). See *About HSP software*.
#' @param sun_indices one-row `data.frame` with numeric columns `"col"` and
#'   `"row"` containing positive, integer-valued cell indices of the solar
#'   disk, numbered from 1, without missing values. Additional columns are
#'   ignored. Columns are identified by name; their order does not matter.
#'
#' @return See *Functions*
#'
#' @references \insertAllCited{}
#'
#' @name hsp_compat
#' @rdname hsp_compat
#' @aliases hsp_read_manual_input hsp_read_opt_sky_coef hsp_write_sky_indices hsp_write_sun_indices
#'
#' @export
#'
#' @examples
#' \dontrun{
#' setwd("~/HSP/Projects/KS_2021_11_18")
#'
#' r <- read_caim("manipulate/DSCN6269.pgm")
#'
#' rho <- ncol(r) / 2
#' hsp_geometry <- list(
#'   model = "radial_projection",
#'   center_of_vision = c(rho, rho),
#'   parameters = equidistant_coef,
#'   normalization_radius = rho,
#'   normalization_angle = 90
#' )
#' g <- set_internal_geometry(c(rho, rho) * 2, hsp_geometry)
#' d <- direction_image(g)
#'
#' # From HSP to R in order to compare ---------------------------------------
#'
#' manual_input <- hsp_read_manual_input(".", "DSCN6269")
#' sun_angles <- manual_input$sun_indices %>% add_directions(., d = d)
#'
#' sky_points <- manual_input$sky_indices %>% add_directions(., d = d)
#' sky_points <- extract_dn(r, sky_points)
#' head(sky_points)
#'
#' rr <- calc_relative_radiance(sky_points, col_id = 3, zenith_dn = manual_input$zenith_dn)
#' head(rr)
#' model <- fit_cie_model(rr, sun_angles, parallel = TRUE)
#' sky <- cie_image(d, model$sun_angles, model$coef) * model$relative_radiance$zenith_dn
#' plot(sky)
#' names(sky) <- "above-canopy dn"
#' plot(r / sky)
#'
#' r <- read_caim("manipulate/DSCN6269.pgm")
#' sky_coef <- hsp_read_opt_sky_coef(".", "DSCN6269")
#' sky_m <- cie_image(d, sun_angles, sky_coef) * manual_input$zenith_dn
#' names(sky_m) <- "above-canopy dn"
#' plot(sky_m)
#' plot(r / sky_m)
#'
#' # From R to HSP  ----------------------------------------------------------
#'
#' w <- area_weight_image(g, hsp_geometry)
#' m <- !is.na(g$off_axis)
#'
#' bin <- binarize_with_thr(r, thr_isodata(r[m], w[m]))
#' seg <- beckers_segmentation(d, 500)
#'
#' sun_indices <- estimate_sun_angles(r, d, bin, seg, w) %>%
#'   internal_directions_from_external() %>%
#'   coord_from_internal_directions(., hsp_geometry) %>% add_indices()
#'
#' hsp_write_sun_indices(sun_indices, ".", "DSCN6269")
#'
#' sky_indices <- compute_features(distance_to_black(bin), seg, d, w) %>%
#'   internal_directions_from_external %>%
#'   coord_from_internal_directions(., hsp_geometry) %>% add_indices()
#' head(sky_indices)
#'
#' hsp_write_sky_indices(sky_indices, ".", "DSCN6269")
#'
#' }
hsp_read_manual_input <- function(path_to_HSP_project, img_name) {
  .check_vector(path_to_HSP_project, "character", 1)
  .check_vector(img_name, "character", 1)
  img_name <- basename(tools::file_path_sans_ext(img_name, compression = FALSE))

  .get_vect_from_plain <- function(file) {
    # Mait Lang, 04.11.2025.
    masks <- NULL
    # Read file content to vector of strings
    masks.vertices <- readLines(file)

    .m <- as.integer(masks.vertices[1])
    if (.m != 0) {
      while (.m >= 1 ) {

        # remove the count from vertex list
        vertices  <- strsplit( gsub("^\\S+ ", "",masks.vertices[.m+1]), split="[ \\.]") %>%
          unlist() %>%
          as.numeric() %>%
          matrix(., ncol = 2, byrow = TRUE)

        if (is.null(masks)) {
          masks <- cbind(object=.m, part=1,vertices , hole=0)
        } else {
          masks <- rbind(masks, cbind(object=.m, part=1,vertices , hole=0))
        }
        .m <- .m -1
      }
      # leave it col and row. The row has to be converted to y coordinate later when the number of image lines is known.
      colnames(masks)[3:4] <- c("col", "row" )

      #Below code is by Gastón
      r <- read_caim(file.path(path_to_HSP_project, "manipulate",
                               paste0(img_name, ".pgm")))
      x <- terra::xFromCol(r, masks[, "col"])
      y <- terra::yFromRow(r, masks[, "row"])
      masks[, c("col", "row")] <- cbind(x, y)
      colnames(masks)[3:4] <- c("x", "y")
      masks <- terra::vect(masks, "polygons", crs = crs(r))
    }
    masks
  }

  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "settings", full.names = TRUE)
  file <- files[grep(img_name, files)]
  settings <- scan(file, "character")
  settings <- settings[c(9, 11:13)]
  settings <- data.frame(
    name = Map(function(x) x[1], strsplit(settings, "=")) %>% unlist(),
    value = Map(function(x) x[2], strsplit(settings, "=")) %>% unlist()
  )

  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "sun", full.names = TRUE)
  file <- files[grep(img_name, files)]
  sun <- scan(file, "character")
  sun <- strsplit(sun, "\\.") %>% unlist() %>% as.numeric()
  sun_mark <- list()
  sun_row_col <- rev(sun)

  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "points", full.names = TRUE)
  file <- files[grep(img_name, files)]
  sky_marks <- scan(file, "character", skip = 1)
  sky_marks <- strsplit(sky_marks, "\\.") %>%
    unlist() %>%
    as.numeric() %>%
    matrix(., ncol = 3, byrow = TRUE) %>%
    as.data.frame(.)
  names(sky_marks) <- c("col", "row", "type" )

  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "statistics", full.names = TRUE)
  file <- files[grep(img_name, files)]
  content <- scan(file, "character", skip = 1, sep = "\n")
  zenith_dn <- content[grep( "Zenith", content)]
  zenith_dn <- strsplit(zenith_dn, "=")[[1]][2] %>%
    sub(",", ".", .) %>% as.numeric()

  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "masks", full.names = TRUE)
  file <- files[grep(img_name, files)]
  masks <- .get_vect_from_plain(file)

  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "polygons", full.names = TRUE)
  file <- files[grep(img_name, files)]
  marked_areas <- .get_vect_from_plain(file)

  sun_indices <- data.frame(row = sun_row_col[1], col = sun_row_col[2])

  list(weight = settings[1,2] %>% as.numeric(),
       max_points = settings[2,2] %>% as.numeric(),
       angle = settings[3,2] %>% as.numeric(),
       point_radius = settings[4,2] %>% as.numeric(),
       sun_indices = sun_indices,
       sky_indices = sky_marks,
       zenith_dn = zenith_dn,
       masks = masks,
       marked_areas = marked_areas)
}

#' @rdname hsp_compat
#' @export
hsp_read_opt_sky_coef <- function(path_to_HSP_project, img_name) {
  .check_vector(path_to_HSP_project, "character", 1)
  .check_vector(img_name, "character", 1)
  img_name <- basename(tools::file_path_sans_ext(img_name, compression = FALSE))

  files <- dir(file.path(path_to_HSP_project, "manipulate"),
               pattern = "opt-parameters", full.names = TRUE)
  file <- files[grep(img_name, files)]
  sky_coef <- scan(file, "character", skip = 1)
  sky_coef <- data.frame(
    name = Map(function(x) x[1], strsplit(sky_coef, "=")) %>% unlist(),
    value = Map(function(x) x[2], strsplit(sky_coef, "=")) %>% unlist()
  )
  sky_coef[c(2, 1, 5, 4, 3), 2] %>% sub(",", ".", .) %>% as.numeric()
}

#' @rdname hsp_compat
#' @export
hsp_write_sky_indices <- function(sky_indices, path_to_HSP_project, img_name) {
  .check_columns(sky_indices, c("row", "col"))
  for (column in c("row", "col")) {
    .check_vector(sky_indices[[column]], "integerish", sign = "positive",
                  name = paste0("sky_indices$", column))
  }
  .check_vector(path_to_HSP_project, "character", 1)
  .check_vector(img_name, "character", 1)

  no <- nrow(sky_indices)

  col.row_coordinates <- paste(sky_indices$col, sky_indices$row, "3", sep = ".")
  col.row_coordinates <- paste(col.row_coordinates, collapse = " ")

  sky_indices <- c(no, col.row_coordinates)
  sky_indices <- data.frame(sky_indices)

  img_name <- basename(tools::file_path_sans_ext(img_name, compression = FALSE))

  utils::write.table(sky_indices, file.path(path_to_HSP_project,
                                           "manipulate",
                                           paste0(img_name, "_points.conf")),
                     quote = FALSE, row.names = FALSE, col.names = FALSE,
                     fileEncoding = "UTF-8", eol = "\n")
}

#' @rdname hsp_compat
#' @export
hsp_write_sun_indices <- function(sun_indices, path_to_HSP_project, img_name) {
  .check_columns(sun_indices, c("row", "col"))
  for (column in c("row", "col")) {
    .check_vector(sun_indices[[column]], "integerish", 1, sign = "positive",
                  name = paste0("sun_indices$", column))
  }
  .check_vector(path_to_HSP_project, "character", 1)
  .check_vector(img_name, "character", 1)

  sun_col_row <- paste(sun_indices$col, sun_indices$row, sep = ".")

  img_name <- basename(tools::file_path_sans_ext(img_name, compression = FALSE))

  utils::write.table(sun_col_row, file.path(path_to_HSP_project,
                                            "manipulate",
                                            paste0(img_name, "_sun.conf")),
                     quote = FALSE, row.names = FALSE, col.names = FALSE,
                     fileEncoding = "UTF-8", eol = "\n")
}
