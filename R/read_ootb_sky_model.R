#' Read files writen by [write_ootb_sky_model()]
#'
#' @inheritParams write_ootb_sky_model
#' @inheritParams sky_grid_segmentation
#'
#' @returns An object of the class _list_ similar to the output of
#'   [ootb_fit_cie_sky_model()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#'
#' path <- system.file("external/ootb_sky.txt", package = "rcaiman")
#' ootb_sky <- read_ootb_sky_model(gsub(".txt", "", path), z, a)
#' }
read_ootb_sky_model <- function(name, z, a) {
  sky <- list()
  ds <- scan(paste0(name, ".txt"), "character")
  sky$model$sun_zenith_azimuth[1] <-
    ds[grep("sun_theta:", ds) + 1] %>% as.numeric()
  sky$model$sun_zenith_azimuth[2] <-
    ds[grep("sun_phi:", ds) + 1] %>% as.numeric()
  sky$model$zenith_dn <- ds[grep("zenith_dn:", ds) + 1] %>% as.numeric()
  sky$model$start[1] <- ds[grep("start_a:", ds) + 1] %>% as.numeric()
  sky$model$start[2] <- ds[grep("start_b:", ds) + 1] %>% as.numeric()
  sky$model$start[3] <- ds[grep("start_c:", ds) + 1] %>% as.numeric()
  sky$model$start[4] <- ds[grep("start_d:", ds) + 1] %>% as.numeric()
  sky$model$start[5] <- ds[grep("start_e:", ds) + 1] %>% as.numeric()
  sky$model$method <- ds[grep("method:", ds) + 1]
  sky$model$coef[1] <- ds[grep("fit_a:", ds) + 1] %>% as.numeric()
  sky$model$coef[2] <- ds[grep("fit_b:", ds) + 1] %>% as.numeric()
  sky$model$coef[3] <- ds[grep("fit_c:", ds) + 1] %>% as.numeric()
  sky$model$coef[4] <- ds[grep("fit_d:", ds) + 1] %>% as.numeric()
  sky$model$coef[5] <- ds[grep("fit_e:", ds) + 1] %>% as.numeric()
  sky$model_validation$r_squared <-
    ds[grep("r_squared:", ds) + 1] %>% as.numeric()
  sky$model_validation$rmse <-
    ds[grep("rmse:", ds) + 1] %>% as.numeric()
  sky$model_validation$mae <-
    ds[grep("mae:", ds) + 1] %>% as.numeric()
  sky$dist_to_black <- ds[grep("dist_to_black:", ds) + 1] %>% as.numeric()
  sky$use_window <- !is.null(sky$dist_to_black)
  sky$min_spherical_dist <- ds[grep("min_spherical_dist:", ds) + 1] %>%
    as.numeric()
  sky$g <- tryCatch(ds[grep("grid,", ds) + 1] %>% as.numeric(),
                    error = function(e)
                            ds[(grep("grid:", ds)+1):
                                  (grep("dist_to_black:", ds)-1)] %>%
                                                      paste(., collapse = " "))

  sky_points <- terra::vect(paste0(name, "_sky_points", ".gpkg"))
  sky_points <- terra::extract(z, sky_points, cells = TRUE)
  sky_points <- terra::rowColFromCell(z, sky_points$cell) %>% as.data.frame()
  colnames(sky_points) <- c("row", "col")
  sky$sky_points <- sky_points

  df <- utils::read.csv2(paste0(name, "_val", ".csv"))
  sky$model_validation$pred <- df$pred
  sky$model_validation$obs <- df$obs

  if (is.numeric(sky$g)) {
    sky$g <- sky_grid_segmentation(z, a, sky$g)
  }

  model <- sky$model
  sky$sky <- cie_sky_image(z, a,
                           model$sun_zenith_azimuth,
                           model$coef) * model$zenith_dn
  names(sky$sky) <- "CIE sky"

  sky
}
