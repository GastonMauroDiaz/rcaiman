.radian2degree <- function(x) x * 180 / pi
.degree2radian <- function(x) x * pi / 180

.nn_spherical_fast <- function(xyz, xyz0, k = 1, search_radius = NULL) {
  dot <- xyz$x * xyz0["x"] + xyz$y * xyz0["y"] + xyz$z * xyz0["z"]
  dot <- pmin(pmax(dot, -1), 1)

  d <- acos(dot)

  k <- min(k, length(d))
  idx <- order(d)[seq_len(k)]

  if (!is.null(search_radius)) {
    keep <- d[idx] <= search_radius
    idx <- idx[keep]
  }

  list(
    index = idx,
    distance = d[idx]
  )
}

.cartesian_to_directions <- function(xyz) {
  x <- xyz[,1]; y <- xyz[,2]; zc <- xyz[,3]
  zen <- acos(zc)
  azi <- atan2(y, x)
  azi[azi < 0] <- azi[azi < 0] + 2*pi
  cbind(zen = zen, azi = azi)
}

.directions_to_cartesian <- function(sampling_points) {
  zen <- sampling_points[, "zen"]
  azi <- sampling_points[, "azi"]
  data.frame(
    x = sin(zen) * cos(azi),
    y = sin(zen) * sin(azi),
    z = cos(zen)
  )
}

.get_max <- function(r){
  if (inherits(r, "SpatRaster")) {
    terra::minmax(r)[2]
  } else {
    max(r[], na.rm = TRUE)
  }
}

.get_min <- function(r) {
  if (inherits(r, "SpatRaster")){
    terra::minmax(r)[1]
  } else {
    min(r[], na.rm = TRUE)
  }
}

.decode_label <- function(label) {
  sector_ID <- trunc(label / 1000)
  ring_ID <- label - sector_ID * 1000
  ds <- data.frame(sector_ID, ring_ID)
  names(ds) <- c("sector_ID", "ring_ID")
  ds
}

.calc_rmse <- function(x) sqrt(mean(x^2))

.assert_file_extension <- function(file_name, required_extension = "tif") {
  ext <- tolower(tools::file_ext(file_name))

  if (ext == "") {
    file_name <- paste0(file_name, ".", required_extension)
  } else if (ext != required_extension) {
    stop(
      paste0(
       "Expected a ",required_extension ," extension (or no extension)."
      ),
      call. = FALSE
    )
  }
  file_name
}

.cores <- function(cores = NULL, logical = FALSE, leave_free = 0) {
  total <- parallel::detectCores(logical = logical) # physical or logical cores?
  if (is.na(total) || total < 1) total <- 1

  # explicit request
  if (!is.null(cores)) {
    cores <- min(total, cores)
  } else {
    # implicit request
    cores <- total - leave_free
  }

  # ensure meaningful parallelism
  if (cores < 2) cores <- 1

  cores
}

.with_cluster <- function(cores, expr) {
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("The 'doParallel' package is required but not available.", call. = FALSE)
  }

  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  force(expr)
}

.resolve_required_column <- function(col_id,
                                     data,
                                     name_col_id = deparse(substitute(col_id)),
                                     name_data = deparse(substitute(data))) {

  if (is.numeric(col_id)) {
    .check_vector(col_id, "numeric", 1, sign = "positive")
    if (col_id > ncol(data)) {
      stop(
        sprintf(
          "`%s` (%s) exceeds the number of columns in `%s` (%s).",
          name_col_id,
          col_id,
          name_data,
          ncol(data)
        )
      )
    }
    return(names(data)[col_id])
  } else if (is.character(col_id)) {
    .check_vector(col_id, "character", 1)
    .check_columns(data, required_columns = col_id, arg_name = name_data)
    return(col_id)
  }
  stop(sprintf("`%s` is not character or numeric.", name_col_id))
}

.noise <- function(w = 1) {
  coef_sd <- apply((rcaiman::cie_table[, 1:5]), 2, sd) * w
  Map(function(i) stats::rnorm(1, 0, coef_sd[i]), 1:5) %>% unlist()
}

.segment_values <- function(seg) {
  .assert_single_layer(seg)
  if (!terra::hasValues(seg)) stop("`seg` must contain labels.", call. = FALSE)
  values <- terra::values(seg, mat = FALSE)
  if (is.logical(values)) values <- as.integer(values)
  .check_vector(values, "integerish", allow_na = TRUE,
                sign = "nonnegative", name = "seg")
  values
}

#' Internal inverse-distance weights
#'
#' @description
#' Calculate normalized inverse-square distance weights shared by radiometric
#' reference estimators in rcaiman and rcaiman.registry. This helper is exported
#' for use across the suite and is not intended for general use.
#'
#' @param distance nonempty finite numeric vector of nonnegative distances,
#'   all in the same units. If any distance is zero, weight is distributed
#'   equally among the zero-distance observations; other weights are zero.
#'
#' @return Numeric vector of weights summing to one, in input order.
#' @keywords internal
#' @export
.inverse_distance_weights <- function(distance) {
  .check_vector(distance, "numeric", sign = "nonnegative")
  if (!length(distance)) stop("`distance` must not be empty.", call. = FALSE)
  zero <- distance == 0
  if (any(zero)) return(as.numeric(zero) / sum(zero))
  # Scaling before squaring avoids overflow close to the reference direction.
  weights <- (min(distance) / distance)^2
  weights / sum(weights)
}
