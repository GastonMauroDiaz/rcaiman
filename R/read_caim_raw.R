.label_cfa <- function(r) {
  r[] <- c(rep(c(0,1), ncol(r)/2), rep(c(2,3), ncol(r)/2)) %>% rep(., nrow(r)/2)
  r
}
.apply_spectral_mapping <- function(bands, spectral_mapping) {
  #bands: named list composed of single-layer rasters

  stopifnot(is.list(bands), is.list(spectral_mapping))
  stopifnot(all(names(spectral_mapping) != ""))

  out <- vector("list", length(spectral_mapping))
  names(out) <- names(spectral_mapping)

  for (nm in names(spectral_mapping)) {

    fun <- spectral_mapping[[nm]]

    # arguments expected by the function
    args <- names(formals(fun))

    if (!all(args %in% names(bands))) {
      stop("Missing bands for mapping '", nm, "'")
    }

    # stack required bands
    r <- terra::rast(bands[args])
    # apply function
    if (length(args) == 1) out[[nm]] <- fun(r[[1]])
    if (length(args) == 2) out[[nm]] <- fun(r[args[1]], r[args[2]])
    if (length(args) == 3) out[[nm]] <- fun(r[args[1]], r[args[2]], r[args[3]])
    if (length(args) == 4) out[[nm]] <- fun(r[args[1]], r[args[2]], r[args[3]], r[args[4]])
  }
  terra::rast(out)
}



#' Read a canopy image from a raw file
#'
#' @description
#' Read unprocessed sensor data from a camera RAW file and split the signal by
#' spectral band according to the in-camera color filter array (CFA). Use this
#' to obtain images with precise radiometry.
#'
#' @details
#' Uses Python [`rawpy`](https://pypi.org/project/rawpy/) through `reticulate`
#' to access sensor data and black-level metadata. Optionally extracts only the
#' blue/cyan band.
#'
#' @section Check Python Accessibility:
#'
#' To ensure that R can access a Python installation, run the following test:
#'
#' ````
#' reticulate::py_eval("1+1")
#'
#' ````
#'
#' If R can access Python successfully, you will see `2` in the console. If not,
#' you will receive instructions on how to install Python.
#'
#' **Hint**: If you are temporarily offline, try setting `Sys.setenv(UV_OFFLINE=1)`.
#' Run `reticulate::py_last_error()` for details. Try to load rawpy explicitly.
#'
#' @section Create a Virtual Environment:
#'
#' After passing the Python accessibility test, create a virtual environment
#' using the following command:
#'
#' ````
#' reticulate::virtualenv_create()
#'
#' ````
#'
#' @section Install `rawpy`:
#'
#' Install the rawpy package within the virtual environment:
#'
#' ````
#' reticulate::py_install("rawpy")
#'
#' ````
#'
#' @section For RStudio Users:
#'
#' If you are an RStudio user who works with projects, you will need a
#' _.Renviron_ file in the root of each project. To create a _.Renviron_ file,
#' follow these steps:
#'
#' * Create a "New Blank File" named ".Renviron" (without an extension) in the
#' project's root directory.
#'
#' * Run bellow code:
#'
#' ````
#' path <- file.path(reticulate::virtualenv_root(),
#' reticulate::virtualenv_list(), "Scripts", "python.exe")
#' paste("RETICULATE_PYTHON =", path)
#'
#' ````
#'
#' * Copy/paste the line from the console (the string between the quotes) into
#' the .Renviron file. This is an example `RETICULATE_PYTHON =
#' ~/.virtualenvs/r-reticulate/Scripts/python.exe`
#'
#' * Do not forget to save the changes
#'
#' By following these steps, users can easily set up their environment to access
#' raw data efficiently, but it is not the only way of doing it, you might know
#' an easier or better one.
#'
#' See the help page of [read_caim()] and [fisheye_to_equidistant()] as a
#' complement to this help page. Further details about raw files can be found in
#' \insertCite{Diaz2024;textual}{rcaiman}.
#'
#' @param path character vector of length one. Path to a file with raw data
#'   (including file extension).
#' @param only_blue logical vector of length one. If `TRUE`, return only the blue/cyan band.
#' @param offset_value numeric vector of length one. Optional black level offsets to replace
#'   [`black_level_per_channel`](https://www.libraw.org/docs/API-datastruct-eng.html#datastream_data:~:text=Per%2Dchannel%20black%20level%20correction)
#'   metadata obtained with `rawpy`.
#' @param cfa_pattern character matrix of two rows by two columns. Declares the
#'   ordered set of color filter elements used by the sensor. Values should
#'   correspond to semantic color labels (e.g. `"Red"`, `"Green"`, `"Blue"`,
#'   `"Yellow"`, `"Cyan"`, `"Magenta"`), and reflect the CFA pattern as produced
#'   by the sensor. These may intentionally override or contradict
#'   embedded metadata when such metadata is known to be incorrect.
#' @param spectral_mapping optional named list. Declares how target spectral
#'   bands should be derived from sensor channels. This argument is intended for
#'   use with `type = "interpretive_constraint"`. Each element of the list must
#'   be named according to the target band (e.g. `"Red"`, `"Green"`, `"Blue"`),
#'   and its value must be a function describing how that band is obtained from
#'   one or more spectral bands as declared in `cfa_pattern`.
#'
#' @references \insertAllCited{}
#'
#' @return Numeric [terra::SpatRaster-class]:
#' \itemize{
#'   \item single-layer if `only_blue = TRUE`.
#'   \item multi-layer if `only_blue = FALSE`, with one layer per color per CFA
#'   color (e.g., R, G, B).
#' }
#' Layers are named according to metadata in the raw file.
#'
#' @seealso [read_caim()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' file_name <- tempfile(fileext = ".NEF")
#' download.file("https://osf.io/s49py/download", file_name, mode = "wb")
#'
#' # Geometric and radiometric corrections -----------------------------------
#' zenith_colrow <- c(1290, 988)/2
#' diameter <- 756
#' z <- zenith_image(diameter, lens("Nikon_FCE9"))
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' caim <- read_caim_raw(file_name, only_blue = TRUE)
#' caim <- crop_caim(caim, zenith_colrow - diameter/2, diameter, diameter)
#' caim <- correct_vignetting(caim, z, c(0.0638, -0.101))
#' caim <- fisheye_to_equidistant(caim, z, a, m, radius = 300,
#'                                k = 1, p = 1, rmax = 100)
#' }
read_caim_raw <- function(path,
                          only_blue = FALSE,
                          offset_value = NULL,
                          cfa_pattern = NULL,
                          spectral_mapping = NULL) {
  .check_vector(path, "character", 1)
  .assert_file_exists(path)
  .check_vector(only_blue, "logical", 1)
  .check_vector(offset_value, "numeric", 1, allow_null = TRUE, sign = "positive")

  # Checks
  if (is.null(cfa_pattern) && !is.null(spectral_mapping)) {
    warning("`spectral_mapping` was ignored because `cfa_pattern` is missing")
  }
  if (!is.null(cfa_pattern) && is.null(spectral_mapping)) {
    warning("`cfa_pattern` was ignored because `spectral_mapping` is missing")
  }
  if (!is.null(cfa_pattern) && only_blue) {
    warning("`only_blue` was ignored because `cfa_pattern` was provided")
  }

  if (!is.null(cfa_pattern) && !is.null(spectral_mapping)) {
    # START code from rcaiman.registry::add_radiometry_spec()
    if (!is.matrix(cfa_pattern) || !is.character(cfa_pattern)) {
      stop("`cfa_pattern` must be a character matrix.")
    }
    if (length(unique(as.vector(cfa_pattern))) != 4) {
      stop(
        "`cfa_pattern` must define exactly four unique band identifiers.",
        call. = FALSE
      )
    }
    if (!is.null(spectral_mapping)) {
      if (!is.list(spectral_mapping) ||
          is.null(names(spectral_mapping)) ||
          any(names(spectral_mapping) == "") ||
          any(!vapply(spectral_mapping, is.function, logical(1)))) {
        stop("spectral_mapping must be a named list of functions.")
      }
      for (nm in names(spectral_mapping)) {

        fun <- spectral_mapping[[nm]]

        fmls <- formals(fun)
        args <- names(fmls)

        # no dots
        if ("..." %in% args) {
          stop(
            "Functions in `spectral_mapping` must not use `...` (found in '", nm, "').",
            call. = FALSE
          )
        }

        # number of arguments
        if (length(args) < 1 || length(args) > 4) {
          stop(
            "Function '", nm,
            "' must have between 1 and 4 arguments.",
            call. = FALSE
          )
        }

        # argument names must match cfa_pattern
        if (!all(args %in% cfa_pattern)) {
          stop(
            "Function '", nm,
            "' uses undefined bands: ",
            paste(setdiff(args, cfa_pattern), collapse = ", "),
            call. = FALSE
          )
        }
      }
    }
    # END
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(paste("Package \"reticulate\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }

  # Processing
  rawpy <- reticulate::import("rawpy")

  img <- rawpy$imread(path)
  r <- set_rcaiman_geometry(rast(img$raw_image))
  black_level <- img$black_level_per_channel
  if (!is.null(offset_value)) {
    black_level[] <- offset_value
  } else {
    message("rawpy reports black levels of ",
            paste(black_level, collapse = ", "),
            " for ", as.character(img$color_desc))
  }

  if (is.null(cfa_pattern)) {
    m <- crop_caim(rast(img$raw_colors))
    color_desc <- as.character(img$color_desc)
    color_desc <- strsplit(color_desc, "") %>% unlist()
    raw_pattern <- c(img$raw_pattern[1,1],
                     img$raw_pattern[1,2],
                     img$raw_pattern[2,2],
                     img$raw_pattern[2,1])

    .sample_r <- function(label) {
      if (length(label) == 1) {
        r[m == label] <-  r[m == label] - black_level[label+1]
        foo <- rast(r) %>% terra::aggregate(., 2)
        foo[] <- r[m == label]
        r <- set_rcaiman_geometry(foo)
      } else {
        r[m == label[1]] <- r[m == label[1]] -
          black_level[label[1]+1]
        r[m == label[2]] <- r[m == label[2]] -
          black_level[label[2]+1]
        foo1 <- foo2 <- rast(r) %>% terra::aggregate(., 2)
        foo1[] <- r[m == label[1]]
        foo2[] <- r[m == label[2]]
        r <- set_rcaiman_geometry(mean(foo1, foo2))
      }
      r
    }
    if (only_blue) {
      i <- grep("B", color_desc)
      if (length(i) == 0) {i <- grep("C", color_desc)}
      if (length(i) == 0) {stop("There is no blue or cyan filter")}
      label <- raw_pattern[i]
      r <- .sample_r(label)
    } else {
      l <- list()
      unique_color_desc <- unique(color_desc)
      for (i in seq_along(unique_color_desc)) {
        color_filter <- grep(unique_color_desc[i], color_desc)
        l[[i]] <- raw_pattern[color_filter]
      }
      r <- lapply(l, .sample_r)
      r <- terra::rast(r)
      names(r) <- unique_color_desc
    }
  } else {
    m <- .label_cfa(r)
    raw_pattern <- matrix(c(0,1,
                            2,3), byrow = TRUE, ncol = 2)
    r <- lapply(as.vector(raw_pattern), function(label) {
      r[m == label] <-  r[m == label] - black_level[label + 1]
      foo <- rast(r) %>% terra::aggregate(., 2)
      foo[] <- r[m == label]
      set_rcaiman_geometry(foo)
    })
    names(r) <- as.vector(cfa_pattern)
    r <- .apply_spectral_mapping(r, spectral_mapping)
  }
  img$close()
  r
}
