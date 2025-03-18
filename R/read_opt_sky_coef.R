#' Read optimized sky coefficients
#'
#' Read optimized CIE sky coefficients stored in an HSP project.
#'
#' Refer to the Details section of function
#' [write_sky_points()].
#'
#' @inheritParams write_sky_points
#'
#' @family HSP Functions
#' @seealso [cie_sky_image()]
#'
#' @return Numeric vector of length five.
#'
#' @export
read_opt_sky_coef <- function(path_to_HSP_project, img_name) {
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
