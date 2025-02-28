#' Write report
#'
#' Write report using [base::cat()] and [base::sink()].
#'
#' @inheritParams interpolate_and_merge
#' @inheritParams base::sink
#'
#' @returns No return value. Called for side effects.
#' @export
#'
#' @family Tool Functions
write_report <- function(ootb_sky, file) {
  sink(file)
  .print <- function(x) cat(x, "\n")  # Usa cat() para evitar [1]

  .print("v0.1")

  paste0(rep("-", 80), collapse = "") %>% .print()
  ootb_sky$model$mle2_output %>% summary() %>% suppressWarnings() %>% print()
  paste0(rep("-", 80), collapse = "") %>% .print()

  ootb_sky$model$start[1] %>% paste("start_a:", .) %>% .print()
  ootb_sky$model$start[2] %>% paste("start_b:", .) %>% .print()
  ootb_sky$model$start[3] %>% paste("start_c:", .) %>% .print()
  ootb_sky$model$start[4] %>% paste("start_d:", .) %>% .print()
  ootb_sky$model$start[5] %>% paste("start_e:", .) %>% .print()
  ootb_sky$model$method %>% unname %>% paste("method:", .) %>% .print()
  ootb_sky$model$coef[1] %>% unname %>% paste("a:", .) %>% .print()
  ootb_sky$model$coef[2] %>% unname %>% paste("b:", .) %>% .print()
  ootb_sky$model$coef[3] %>% unname %>% paste("c:", .) %>% .print()
  ootb_sky$model$coef[4] %>% unname %>% paste("d:", .) %>% .print()
  ootb_sky$model$coef[5] %>% unname %>% paste("e:", .) %>% .print()
  ootb_sky$model$sun_coord$zenith_azimuth[1] %>% paste("sun_theta:", .) %>% .print()
  ootb_sky$model$sun_coord$zenith_azimuth[2] %>% paste("sun_phi:", .) %>% .print()
  ootb_sky$model$zenith_dn %>% paste("zenith_dn:", .)%>% .print()
  ootb_sky$sky_points %>% nrow() %>% paste("sky_points_no:", .) %>% .print()
  ootb_sky$sky_points$outliers %>% sum() %>% paste("outliers_no:", .) %>%  .print()
  names(ootb_sky$sky) %>% paste("sky_type:", .) %>% .print()
  ootb_sky$model_validation$rmse %>% paste("RMSE:", .) %>% .print()

  paste0(rep("-", 80), collapse = "") %>% .print()
  ootb_sky$model_validation$reg %>% summary() %>% print()
  paste0(rep("-", 80), collapse = "") %>% .print()
  sink()
}
