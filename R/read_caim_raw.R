#' Read a canopy image from a raw file
#'
#' Function that complements [read_caim()]
#'
#' This function facilitates the integration of the `rawpy` Python package into
#' the R environment via the `reticulate` package. This integration allows
#' `rcaiman` to access and pre-process raw data.
#'
#' Here is a step-by-step guide to assist users in setting up the environment
#' for efficient processing:
#'
#' ## Check Python Accessibility:
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
#' ## Create a Virtual Environment:
#'
#' After passing the Python accessibility test, create a virtual environment
#' using the following command:
#'
#' ````
#' reticulate::virtualenv_create()
#'
#' ````
#'
#' ## Install `rawpy`:
#'
#' Install the rawpy package within the virtual environment:
#'
#' ````
#' reticulate::py_install("rawpy")
#'
#' ````
#'
#' ## For RStudio Users:
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
#' raw data efficiently, but it is not the only way of doing it.
#'
#'
#' @param path Character vector of length one.Path to a raw file, including file
#'   extension.
#' @inheritParams read_caim
#' @inheritParams fisheye_to_equidistant
#' @inheritParams expand_noncircular
#' @param only_blue Logical vector of length one. If `TRUE`, only values from
#'   the blue or cyan wavelength will be processed.
#' @param offset_value numeric vector. This values will replace the
#' `black_level_per_channel` metadata obtained with `rawpy`.
#'
#' @return An object from class [SpatRaster-class]. Single-layer raster if
#'   `only_blue` is equal to `TRUE`. Otherwise, a raster with as many layers as
#'   there are distinct colors in the Color Filter Array. Layer names are taken
#'   from the color description metadata.
#' @export
#' @family Tool Functions
read_caim_raw <- function(path = NULL,
                          z = NULL,
                          a = NULL,
                          zenith_colrow = NULL,
                          radius = 700,
                          rmax = 100,
                          k = 1,
                          p = 1,
                          only_blue = FALSE,
                          offset_value = NULL) {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(paste("Package \"reticulate\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }

  rawpy <- reticulate::import("rawpy")

  img <- rawpy$imread(path)
  r <- crop_caim(rast(img$raw_image))
  m <- crop_caim(rast(img$raw_colors))
  black_level <- img$black_level_per_channel
  if (!is.null(offset_value)) black_level[] <- offset_value

  if (is.null(z) & is.null(a) & is.null(zenith_colrow)) {
    .fun <- function(label) {
      if (length(label) == 1) {
        r[m == label] <-  r[m == label] - black_level[label+1]
        foo <- rast(r) %>% terra::aggregate(., 2)
        foo[] <- r[m == label]
        r <- crop_caim(foo)
      } else {
        r[m == label[1]] <- r[m == label[1]] -
          black_level[label[1]+1]
        r[m == label[2]] <- r[m == label[2]] -
          black_level[label[2]+1]
        foo1 <- foo2 <- rast(r) %>% terra::aggregate(., 2)
        foo1[] <- r[m == label[1]]
        foo2[] <- r[m == label[2]]
        r <- crop_caim(mean(foo1, foo2))
      }
      warning("The output image is not in the standard equidistant projection.")
      r
    }
  } else {
    .fun <- function(label) {
      if (length(label) == 1) {
        r[m == label] <-  r[m == label] - black_level[label+1]
        m[is.na(m)] <- label
        r <- fisheye_to_equidistant(r, z, a, m == label , radius = radius,
                                    k = k, p = p, rmax = rmax)
      } else {
        r[m == label[1]] <- r[m == label[1]] -
          black_level[label[1]+1]
        r[m == label[2]] <- r[m == label[2]] -
          black_level[label[2]+1]
        m[is.na(m)] <- label[1]
        r <- fisheye_to_equidistant(r, z, a, m == label[1] |  m == label[2],
                                    radius = radius,
                                    k = k, p = p, rmax = rmax)
      }
      r[r == -1000] <- NA
      r
    }

    diameter <- terra::ncol(z)
    if (diameter > nrow(r)) {
      m <- expand_noncircular(m, z, zenith_colrow)
      r <- expand_noncircular(r, z, zenith_colrow)
      r[is.na(r)] <- -1000
    } else {
      r <- crop_caim(r,
                     upper_left = zenith_colrow - diameter/2,
                     width = diameter,
                     height = diameter)
      m <- crop_caim(m,
                     upper_left = zenith_colrow - diameter/2,
                     width = diameter,
                     height = diameter)
    }
  }

  color_desc <- as.character(img$color_desc)
  color_desc <- strsplit(color_desc, "") %>% unlist()
  raw_pattern <- c(img$raw_pattern[1,1],
                   img$raw_pattern[1,2],
                   img$raw_pattern[2,2],
                   img$raw_pattern[2,1])
  if (only_blue) {
    i <- grep("B", color_desc)
    if (length(i) == 0) {i <- grep("C", color_desc)}
    if (length(i) == 0) {stop("There is no blue or cyan filter.")}
    label <- raw_pattern[i]
    r <- .fun(label)
  } else {
    l <- list()
    unique_color_desc <- unique(color_desc)
    for (i in seq_along(unique_color_desc)) {
      color_filter <- grep(unique_color_desc[i], color_desc)
      l[[i]] <- raw_pattern[color_filter]
    }
    r <- Map(.fun, l)
    r <- terra::rast(r)
    names(r) <- unique_color_desc
  }
  img$close()
  r
}
