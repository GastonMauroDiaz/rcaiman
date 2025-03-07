#' Read files writen by [write_ootb_sky_model()]
#'
#' @inheritParams write_ootb_sky_model
#'
#' @returns An object of the class _list_ similar to the output of
#'   [ootb_fit_cie_sky_model()].
#'
#' @family Tool Functions
#' @export
#'
#' @examples
#' path <- system.file("external/ootb_sky.txt", package = "rcaiman")
#' ootb_sky <- read_ootb_sky_model(gsub(".txt", "", path))
read_ootb_sky_model <- function(name) {

  ootb_sky <- list()
  ds <- scan(paste0(name, ".txt"), "character")
  ootb_sky$model$sun_coord$zenith_azimuth[1] <-
    ds[grep("sun_theta:", ds) + 1] %>% as.numeric()
  ootb_sky$model$sun_coord$zenith_azimuth[2] <-
    ds[grep("sun_phi:", ds) + 1] %>% as.numeric()
  ootb_sky$model$zenith_dn <- ds[grep("zenith_dn:", ds) + 1] %>% as.numeric()
  ootb_sky$model$start[1] <- ds[grep("start_a:", ds) + 1] %>% as.numeric()
  ootb_sky$model$start[2] <- ds[grep("start_b:", ds) + 1] %>% as.numeric()
  ootb_sky$model$start[3] <- ds[grep("start_c:", ds) + 1] %>% as.numeric()
  ootb_sky$model$start[4] <- ds[grep("start_d:", ds) + 1] %>% as.numeric()
  ootb_sky$model$start[5] <- ds[grep("start_e:", ds) + 1] %>% as.numeric()
  ootb_sky$model$method <- ds[grep("method:", ds) + 1]
  ootb_sky$model$coef[1] <- ds[grep("fit_a:", ds) + 1] %>% as.numeric()
  ootb_sky$model$coef[2] <- ds[grep("fit_b:", ds) + 1] %>% as.numeric()
  ootb_sky$model$coef[3] <- ds[grep("fit_c:", ds) + 1] %>% as.numeric()
  ootb_sky$model$coef[4] <- ds[grep("fit_d:", ds) + 1] %>% as.numeric()
  ootb_sky$model$coef[5] <- ds[grep("fit_e:", ds) + 1] %>% as.numeric()
  ootb_sky$model_validation$r_squared <-
    ds[grep("r_squared:", ds) + 1] %>% as.numeric()

  ootb_sky$sky_points <- utils::read.csv2(paste0(name, "_fit", ".csv"))[,-1]

  df <- data.frame(predicted = ootb_sky$model_validation$predicted,
                   observed = ootb_sky$model_validation$observed)
  df <- utils::read.csv2(paste0(name, "_val", ".csv"))
  ootb_sky$model_validation$predicted <- df$predicted
  ootb_sky$model_validation$observed <- df$observed

  ootb_sky
}
