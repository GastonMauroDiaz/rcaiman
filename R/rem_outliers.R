#' Remove local outliers from sampled directions
#'
#' @description
#' Remove observations whose values are outlying relative to nearby samples
#' on the sphere, optionally after removing a local sky-brightness trend.
#'
#' @details
#' Each row is evaluated against its `k` nearest other observations within
#' `search_radius`. The query is excluded by its row position, not by assuming
#' that it is the first neighbor. Repeated directions remain distinct
#' observations and can contribute to one another's neighborhoods. Distance
#' ties follow the ordering used by [nn_spherical()].
#'
#' The decision uses [is_outlier()]: compare the query value with the median and
#' scaled median absolute deviation (MAD) of its neighbors. Angular distance
#' selects the neighborhood; it is not itself the variable being tested.
#' Every neighbor contributes one observation, without distance or area weights.
#' If fewer than `k` other observations are available, retain the query.
#' A zero neighbor MAD also retains the query, following [is_outlier()].
#'
#' All decisions use the original table. Removing one observation does not
#' change the neighborhoods or decisions for others; the filter is not iterative.
#' Missing directions provide no neighborhood and are retained. They do not
#' contribute as neighbors. At a pole, azimuth may be missing for the angular
#' search, as described in [nn_spherical()]. Values selected by `col_id` must
#' be finite, including rows with missing directions.
#'
#' With `detrend = TRUE`, fit [fit_cone_shaped_model()] using
#' `method = "zenith_n_azimuth"` to the neighbors, excluding the query. Compare
#' the query's observed-minus-predicted value with the neighbor residuals using
#' the same MAD rule. The fit uses the variable selected by `col_id` and external
#' `zenith` and `azimuth` in radians. These external columns are required even
#' if internal directions are chosen for the neighbor search.
#'
#' Detrending requires `k >= 20`. If a full-rank fit cannot be obtained,
#' required external angles are missing, or fitting/prediction fails, the
#' calculation stops with the input row and cause. It never substitutes a test
#' on the original values: `cutoff_side` applies to residuals when detrending
#' is requested. Warnings from usable fits are reported with their input rows,
#' including in parallel execution. A missing azimuth at zenith permits the
#' neighbor search but not evaluation of the azimuth-dependent trend there.
#' Queries with fewer than `k` neighbors are still retained without testing.
#'
#' @param sampling_directions data frame containing the angular columns selected
#'   by `angle_columns`, in radians, and the numeric values selected by `col_id`.
#'   For detrending, also supply external `zenith` and `azimuth`. Cell indices
#'   are not required. Additional columns are preserved unchanged.
#' @param col_id character or integer-like numeric vector of length one. Name
#'   or position of the column containing the finite values to evaluate.
#' @param k integer-like numeric vector of length one, at least 3. Number of
#'   other observations required for each comparison. The query is not counted.
#'   Must be at least 20 when `detrend = TRUE`.
#' @param detrend logical vector of length one. Whether to remove a local
#'   cone-shaped trend before testing. A failed local trend stops the calculation;
#'   see Details.
#' @inheritParams nn_spherical
#' @inheritParams is_outlier
#'
#' @return Data frame containing the retained rows of `sampling_directions`,
#'   preserving their columns, row names, and original order. An empty input
#'   returns an empty table with the same columns.
#' @seealso [add_directions()], [is_outlier()], [fit_cone_shaped_model()],
#'   [apply_locally()]
#' @export
#' @examples
#' # In an image workflow, obtain these angular columns with add_directions()
#' # and the sampled band with extract_dn(). Here the data are synthetic.
#' samples <- external_directions_from_internal(fibonacci_directions(10))
#' samples$id <- seq_len(nrow(samples))
#' samples$Blue <- with(samples,
#'   100 + 40 * zenith^2 + 5 * cos(azimuth) + sin(id))
#' samples$Blue[120] <- samples$Blue[120] + 100
#' retained <- rem_outliers(samples, "Blue", k = 20, search_radius = 30,
#'                          detrend = TRUE)
#' setdiff(samples$id, retained$id)
#' plot(samples$id, samples$Blue, col = "red", pch = 16)
#' points(retained$id, retained$Blue, col = "black", pch = 16)
rem_outliers <- function(sampling_directions,
                         col_id,
                         k = 20,
                         search_radius = 20,
                         laxity = 2,
                         cutoff_side = "both",
                         detrend = FALSE,
                         parallel = FALSE,
                         cores = NULL,
                         logical = TRUE,
                         leave_free = 1,
                         angle_columns = c("zenith", "azimuth")) {
  .check_columns(sampling_directions)
  if (anyDuplicated(names(sampling_directions))) {
    stop("`sampling_directions` must have unique column names.", call. = FALSE)
  }
  if (is.numeric(col_id)) .check_vector(col_id, "integerish", 1, sign = "positive")
  col_id <- .resolve_required_column(col_id, sampling_directions)
  observed <- sampling_directions[[col_id]]
  .check_vector(observed, "numeric", name = col_id)
  .check_vector(k, "integerish", 1, sign = "positive")
  if (k < 3) stop("`k` must be at least 3.", call. = FALSE)
  .check_vector(detrend, "logical", 1)
  .check_vector(laxity, "numeric", 1, sign = "positive")
  .assert_choice(cutoff_side, c("both", "left", "right"))
  if (detrend) {
    if (k < 20) stop("`k` must be at least 20 when `detrend = TRUE`.", call. = FALSE)
    .spherical_neighbor_vectors(sampling_directions, c("zenith", "azimuth"))
  }

  # Request one extra position, then exclude the query by identity. With tied
  # directions the query may occur anywhere, or outside the returned positions.
  neighbors <- nn_spherical(sampling_directions, sampling_directions,
    k = k + 1, search_radius = search_radius, rule = "any",
    parallel = parallel, cores = cores, logical = logical, leave_free = leave_free,
    angle_columns = angle_columns)
  for (i in seq_along(neighbors)) {
    hit <- neighbors[[i]]
    take <- head(which(hit$index != i), k)
    neighbors[[i]] <- list(index = hit$index[take], distance = hit$distance[take])
  }

  if (!nrow(sampling_directions)) return(sampling_directions)
  evaluate <- function(data, distance, i) {
    warnings <- character()
    if (nrow(data) < k) return(list(keep = TRUE, warnings = warnings))
    query_value <- observed[i]
    reference_values <- data[[col_id]]
    if (detrend) {
      residuals <- tryCatch(withCallingHandlers({
        fit <- fit_cone_shaped_model(data, col_id = col_id,
          method = "zenith_n_azimuth", return = "all")
        if (is.null(fit)) {
          cause <- if (length(warnings)) paste(unique(warnings), collapse = "; ") else
            "No usable local trend model was returned."
          stop(cause, call. = FALSE)
        }
        prediction <- as.numeric(stats::predict(fit$model,
          newdata = sampling_directions[i, , drop = FALSE]))
        reference <- as.numeric(stats::residuals(fit$model))
        if (length(prediction) != 1L || !is.finite(prediction)) {
          stop("The local trend prediction is not a finite scalar; check the query's external angles.",
               call. = FALSE)
        }
        if (length(reference) != nrow(data) || any(!is.finite(reference))) {
          stop("The local trend did not return finite residuals for every neighbor.", call. = FALSE)
        }
        list(query = query_value - prediction, reference = reference)
      }, warning = function(w) {
        # Carry warnings back to the caller, including from parallel workers.
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        stop("Detrending failed: ", conditionMessage(e), call. = FALSE)
      })
      query_value <- residuals$query
      reference_values <- residuals$reference
    }
    list(keep = !is_outlier(query_value, reference_values,
                            laxity = laxity, cutoff_side = cutoff_side),
         warnings = warnings)
  }
  # A separate query table leaves all user column names, including output, free.
  queries <- data.frame(query = seq_len(nrow(sampling_directions)))
  result <- apply_locally(sampling_directions, queries, neighbors, evaluate,
    parallel = parallel, cores = cores, logical = logical, leave_free = leave_free)
  for (i in seq_along(result)) {
    if (length(result[[i]]$warnings)) {
      warning("Detrending at sampling_directions row ", i, ": ",
              paste(unique(result[[i]]$warnings), collapse = "; "), call. = FALSE)
    }
  }
  keep <- vapply(result, `[[`, logical(1), "keep")
  sampling_directions[keep, , drop = FALSE]
}
