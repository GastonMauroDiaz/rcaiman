
# .check_*() verify structure. Classes of the base pkgs. ------------------


# .assert_*(), for semantic objects -------------------------------------

.assert_area_weight_image <- function(w) {
  .assert_single_layer(w)
  .check_vector(terra::values(w, mat = FALSE), "numeric", allow_na = TRUE,
                sign = "nonnegative", name = "w")
  invisible(TRUE)
}

.assert_cie_model <- function(model) {
  if (!is.list(model)) {
    stop("`model` must be a list returned by `fit_cie_model()`. ")
  }
  required_names <- c("relative_radiance", "opt_result", "coef", "sun_angles",
                      "method", "start" , "metric")
  if (!all(required_names %in% names(model))) {
    stop(sprintf("`model` must contain %s.",
                 paste(sprintf('"%s"', required_names), collapse = ", ")))
  }
}

.assert_neighbors <- function(x, n_samples, n_queries) {
  if (!is.list(x) || length(x) != n_queries) {
    stop("`neighbors` must be a list with one element per query_data row.",
         call. = FALSE)
  }
  for (i in seq_along(x)) {
    hit <- x[[i]]
    if (!is.list(hit) || !all(c("index", "distance") %in% names(hit)) ||
        anyDuplicated(names(hit))) {
      stop("`neighbors[[", i, "]]` must contain `index` and `distance`.",
           call. = FALSE)
    }
    .check_vector(hit$index, "integerish", sign = "positive",
                  name = paste0("neighbors[[", i, "]]$index"))
    if (any(hit$index > n_samples)) {
      stop("`neighbors[[", i, "]]$index` exceeds the number of sampling_data rows.",
           call. = FALSE)
    }
    .check_vector(hit$distance, "numeric", length(hit$index), sign = "nonnegative",
                  name = paste0("neighbors[[", i, "]]$distance"))
  }
  invisible(TRUE)
}
