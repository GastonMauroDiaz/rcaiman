#' Internal validation helpers
#'
#' @description
#' These functions perform type and structure checks across related packages.
#' They are exported for use in Package \pkg{rcaiman.registry} but are not intended
#' for general use.
#'
#' @param x,y,m,path,allowed,type,length,allow_null,allow_na,allow_inf,multiple,case_sensitive,sign,name,x_name,y_name,mapping,cfa_pattern,min_allowed,max_allowed,open_upper,label,required_columns,reserved_columns,arg_name,parameters,projection_coef,radius,angle,max_angle,max_off_axis_angle,squared
#' Parameters used across various validation helpers.
#' @param geometry named list following the contract of
#'   [set_internal_geometry()]. Additional fields and list classes are accepted.
#'   `.check_internal_geometry()` checks the model without constructing a raster
#'   or modifying the supplied list.
#'
#' @return
#' Returns `invisible(TRUE)` on success. Raises an error otherwise.
#'
#' @name internal-validation
#' @keywords internal
NULL

#' @rdname internal-validation
#' @export
.check_vector <- function(x,
                          type = c("numeric", "character", "logical",
                                   "integerish", "even_integerish",
                                   "odd_integerish", "date"),
                          length = NULL,
                          allow_null = FALSE,
                          allow_na = FALSE,
                          allow_inf = FALSE,
                          sign = c("any", "positive", "nonnegative",
                                   "negative", "nonpositive"),
                          name = deparse(substitute(x))
) {
  .is_integerish <- function(v) {
    if (is.list(v)) return(FALSE)
    if (!is.vector(v) || !is.numeric(v)) return(FALSE)
    if (any(is.infinite(v), na.rm = TRUE)) return(FALSE)
    all(round(v[!is.na(v)]) == v[!is.na(v)])
  }

  .is_even <- function(v) {
    if (!is.vector(v) || !is.numeric(v)) return(FALSE)
    v <- v[!is.na(v)]
    all(round(v/2) == v/2)
  }

  if (is.null(x)) {
    if (allow_null) return(invisible(TRUE))
    stop(sprintf("`%s` cannot be `NULL`.", name), call. = FALSE)
  }

  if (is.list(x)) {
    stop(sprintf("`%s` must be an atomic vector, not a list.", name),
         call. = FALSE)
  }

  if (allow_na && is.vector(x) && all(is.na(x))) {
    if (!is.null(length) && length(x) != length) {
      stop(sprintf("`%s` must have length %d.", name, length), call. = FALSE)
    }
    return(invisible(TRUE))
  }

  type <- match.arg(type)
  sign <- match.arg(sign)

  ok_type <- switch(
    type,
    numeric         = is.vector(x) && is.numeric(x),
    character       = is.vector(x) && is.character(x),
    logical         = is.vector(x) && is.logical(x),
    integerish      = .is_integerish(x),
    even_integerish = .is_integerish(x) && .is_even(x),
    odd_integerish = .is_integerish(x) && !.is_even(x),
    date            = inherits(x, "Date")
  )

  type_msg <- switch(
    type,
    numeric         = "a numeric vector",
    character       = "a character vector",
    logical         = "a logical vector",
    integerish      = "an integer-like numeric vector",
    even_integerish = "an even integer-like numeric vector",
    odd_integerish  = "an odd integer-like numeric vector",
    date            = "a Date object (use as.Date())"
  )

  len_msg <- if (!is.null(length)) sprintf(" of length %d", length) else ""

  sign_msg <- switch(
    sign,
    any          = "",
    positive     = " with positive values only",
    nonnegative  = " with non-negative values only",
    negative     = " with negative values only",
    nonpositive  = " with non-positive values only"
  )

  if (!ok_type) {
    stop(sprintf("`%s` must be %s%s%s.", name, type_msg, len_msg, sign_msg),
         call. = FALSE)
  }

  if (!is.null(length) && length(x) != length) {
    stop(sprintf("`%s` must be %s%s%s.", name, type_msg, len_msg, sign_msg),
         call. = FALSE)
  }

  if (!allow_na && anyNA(x)) {
    stop(sprintf("`%s` cannot contain `NA`.", name), call. = FALSE)
  }

  if(!allow_inf && any(is.infinite(x))) {
    stop(sprintf("`%s` cannot contain `Inf`.", name), call. = FALSE)
  }

  if (sign != "any" &&
      type %in% c("numeric", "integerish", "even_integerish", "odd_integerish")) {

    vals <- as.numeric(x)
    vals <- vals[!is.na(vals)]

    ok_sign <- switch(
      sign,
      positive     = all(vals > 0),
      nonnegative  = all(vals >= 0),
      negative     = all(vals < 0),
      nonpositive  = all(vals <= 0)
    )

    if (!isTRUE(ok_sign)) {
      stop(sprintf("`%s` must be %s%s%s.", name, type_msg, len_msg, sign_msg),
           call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.assert_file_exists <- function(path, name = deparse(substitute(path))) {
  if (is.null(path))
    stop(sprintf("`%s` cannot be `NULL`.", name), call. = FALSE)
  if (!is.character(path) || length(path) == 0)
    stop(sprintf("`%s` must be a character vector of file path(s).", name), call. = FALSE)
  if (anyNA(path))
    stop(sprintf("`%s` cannot contain `NA`.", name), call. = FALSE)

  ok <- file.exists(path)
  if (!all(ok)) {
    missing <- path[!ok]
    if (length(missing) == 1) {
      stop(sprintf("File '%s' does not exist.", missing), call. = FALSE)
    } else {
      stop(sprintf("These files do not exist: %s",
                   paste(sprintf("'%s'", missing), collapse = ", ")),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.assert_choice <- function(x, allowed,
                           name = deparse(substitute(x)),
                           allow_null = FALSE,
                           multiple = FALSE,
                           case_sensitive = TRUE) {

  # sanity on allowed
  if (is.null(allowed) || !is.character(allowed) || length(allowed) < 1 || anyNA(allowed))
    stop("`allowed` must be a non-empty character vector without `NA`.", call. = FALSE)

  q <- function(v) sprintf('"%s"', v)
  collapse_or <- function(v) {
    vq <- q(v)
    n <- length(vq)
    if (n == 1) vq
    else if (n == 2) paste(vq, collapse = " or ")
    else paste0(paste(vq[-n], collapse = ", "), ", or ", vq[n])
  }

  if (is.null(x)) {
    if (allow_null) return(invisible(TRUE))
    stop(sprintf("`%s` cannot be `NULL`.", name), call. = FALSE)
  }

  if (!is.character(x))
    stop(sprintf("`%s` must be a character vector.", name), call. = FALSE)
  if (anyNA(x))
    stop(sprintf("`%s` cannot contain `NA`.", name), call. = FALSE)

  if (!multiple && length(x) != 1)
    stop(sprintf("`%s` must be a character vector of length one.", name), call. = FALSE)

  # case handling
  if (!case_sensitive) {
    x <- tolower(x)
    allowed <- tolower(allowed)
  }

  ok <- if (multiple) all(x %in% allowed) else x %in% allowed

  if (!ok) {
    if (multiple) {
      stop(sprintf("`%s` must contain only allowed values: %s.",
                   name, paste(q(allowed), collapse = ", ")), call. = FALSE)
    } else {
      msg <- if (length(allowed) == 2)
        sprintf("either %s", collapse_or(allowed))
      else
        sprintf("one of: %s", paste(q(allowed), collapse = ", "))
      stop(sprintf("`%s` must be %s.", name, msg), call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.assert_spatraster <- function(x, squared = FALSE, allow_na = TRUE, name = deparse(substitute(x))) {
  if (!inherits(x, "SpatRaster"))
    stop(sprintf("`%s` must be a SpatRaster class.", name), call. = FALSE)
  if (squared && nrow(x) != ncol(x))
    stop(sprintf("`%s` has not equal number of row and columns.", name))
  has_values <- suppressWarnings(.get_min(x))
  has_na <- is.infinite(has_values) || any(is.na(x[]))
  if (!allow_na & has_na)
    stop(sprintf("`%s` must be NA free.", name), call. = FALSE)
  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.assert_single_layer <- function(x, squared = FALSE, allow_na = TRUE, name = deparse(substitute(x))) {
  .assert_spatraster(x, squared, allow_na, name)
  if (terra::nlyr(x) != 1)
    stop(sprintf("`%s` must have a single layer.", name), call. = FALSE)
  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.assert_same_geom <- function(x, y,
                              x_name = deparse(substitute(x)),
                              y_name = deparse(substitute(y))) {
  if (!terra::compareGeom(x, y, stopOnError = FALSE, messages = FALSE))
    stop(sprintf("`%s` and `%s` must have the same geometry.", x_name, y_name), call. = FALSE)
  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.assert_logical_mask <- function(m, squared = FALSE, name = deparse(substitute(m))) {
  .assert_single_layer(m, squared = squared, name = name)
  if (!is.logical(m[]))
    stop(sprintf("`%s` must be logical.", name), call. = FALSE)
  if (anyNA(m[]))
    stop(sprintf("`%s` must be NA-free.", name), call. = FALSE)
  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.assert_raster_range <- function(x,
                                 min_allowed,
                                 max_allowed,
                                 open_upper = FALSE,
                                 label = "Values") {
  r <- terra::minmax(x)

  if (is.na(r[1]) && is.na(r[2])) {
    stop("`x` contains only NA values.", call. = FALSE)
  }

  if (open_upper) {
    invalid <- r[1] < min_allowed || r[2] >= max_allowed
  } else {
    invalid <- r[1] < min_allowed || r[2] > max_allowed
  }

  if (invalid) {
    stop(
      sprintf(
        "%s must be in [%s, %s%s].",
        label,
        min_allowed,
        max_allowed,
        if (open_upper)
          ")"
        else
          "]"
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' @rdname internal-validation
#' @export
.check_columns <- function(x,
                           required_columns = NULL,
                           reserved_columns = NULL,
                           arg_name = deparse(substitute(x))) {

  if (!is.data.frame(x)) {
    stop(sprintf("`%s` argument must be a data.frame.", arg_name), call. = FALSE)
  }

  nms <- names(x)

  if (!is.null(required_columns)) {
    missing_cols <- setdiff(required_columns, nms)
    if (length(missing_cols) > 0) {
      stop(
        sprintf(
          "`%s` argument must contain required columns: %s.",
          arg_name,
          paste(missing_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(reserved_columns)) {
    present_cols <- intersect(reserved_columns, nms)
    if (length(present_cols) > 0) {
      stop(
        sprintf(
          "`%s` argument already contains reserved columns: %s.",
          arg_name,
          paste(present_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  invisible(x)
}

#' @rdname internal-validation
#' @export
.check_cfa_pattern <- function(cfa_pattern) {
  if (is.null(cfa_pattern)) {
    stop(
      "For `type = 'interpretive_constraint'`, `cfa_pattern` must be provided.",
      call. = FALSE
    )
  }
  if (!is.matrix(cfa_pattern) || !is.character(cfa_pattern)) {
    stop("`cfa_pattern` must be a character matrix.")
  }
  if (!identical(dim(cfa_pattern), c(2L, 2L))) {
    stop("`cfa_pattern` must be a 2x2 matrix.", call. = FALSE)
  }
  if (anyNA(cfa_pattern)) {
    stop("`cfa_pattern` must not contain missing values.", call. = FALSE)
  }
  if (length(unique(as.vector(cfa_pattern))) != 4) {
    stop(
      "`cfa_pattern` must define exactly four unique band identifiers.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


#' @rdname internal-validation
#' @export
.check_spectral_mapping <- function(mapping, cfa_pattern) {

  if (!is.list(mapping) || is.null(names(mapping)) || any(names(mapping) == "")) {
    stop("`spectral_mapping` must be a named list.", call. = FALSE)
  }

  allowed_symbols <- unique(as.vector(cfa_pattern))
  allowed_ops <- c("+", "-", "*", "/", "(", ")")

  for (nm in names(mapping)) {

    expr_text <- mapping[[nm]]

    if (!is.character(expr_text) || length(expr_text) != 1) {
      stop("Each `spectral_mapping` element must be a character vector of length one.", call. = FALSE)
    }

    expr <- tryCatch(
      parse(text = expr_text)[[1]],
      error = function(e) {
        stop(sprintf("Invalid spectral expression for '%s'.", nm), call. = FALSE)
      }
    )

    tokens <- all.names(expr, functions = TRUE)

    for (tok in tokens) {

      if (tok %in% allowed_symbols) next
      if (tok %in% allowed_ops) next

      stop(
        sprintf(
          "Invalid token '%s' in spectral mapping for '%s'. Only symbols '%s', numbers, and operators '%s' are allowed.",
          tok, nm, paste0(allowed_symbols, collapse = ", "), paste0(allowed_ops, collapse = ", ")
        ),
        call. = FALSE
      )

    }

  }

  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.check_internal_geometry <- function(geometry) {
  nms <- names(geometry)
  if (!is.list(geometry) || is.null(nms) || anyNA(nms) ||
      any(nms == "") || anyDuplicated(nms)) {
    stop("`geometry` must be a list with unique, nonempty field names.",
         call. = FALSE)
  }
  .check_vector(geometry[["center_of_vision"]], "integerish", 2,
                name = "geometry$center_of_vision")
  .check_vector(geometry[["model"]], "character", 1, name = "geometry$model")
  .assert_choice(geometry[["model"]],
                 c("radial_projection", "anisotropic_projection"))
  if (geometry[["model"]] == "radial_projection") {
    max_angle <- if ("max_off_axis_angle" %in% nms)
      geometry[["max_off_axis_angle"]] else 90
    .check_internal_radial_model(
      geometry[["parameters"]], geometry[["normalization_radius"]],
      geometry[["normalization_angle"]], max_angle
    )
  } else {
    .check_internal_anisotropic_parameters(geometry[["parameters"]])
  }
  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.check_internal_radial_model <- function(parameters, radius, angle, max_angle) {
  .check_vector(parameters, "numeric")
  if (!length(parameters)) stop("`parameters` cannot be empty.", call. = FALSE)
  .check_vector(radius, "numeric", 1, sign = "positive",
                name = "normalization_radius")
  .check_vector(angle, "numeric", 1, sign = "positive",
                name = "normalization_angle")
  if (angle > 180) {
    stop("`normalization_angle` must not exceed 180 degrees.", call. = FALSE)
  }
  result <- validate_projection_function(parameters, tol = 0.01,
                                         normalization_angle = angle,
                                         max_off_axis_angle = max_angle)
  if (result$status != "valid") stop(result$reason, call. = FALSE)
  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.check_radial_inversion <- function(projection_coef, max_off_axis_angle) {
  .check_vector(projection_coef, "numeric")
  if (!length(projection_coef)) {
    stop("`projection_coef` cannot be empty.", call. = FALSE)
  }
  .check_vector(max_off_axis_angle, "numeric", 1, sign = "positive")
  if (max_off_axis_angle > 180) {
    stop("`max_off_axis_angle` must not exceed 180 degrees.", call. = FALSE)
  }
  upper <- max_off_axis_angle * pi / 180
  derivative <- projection_coef * seq_along(projection_coef)
  roots <- polyroot(derivative)
  critical <- Re(roots)[abs(Im(roots)) < 1e-7 &
                          Re(roots) > 0 & Re(roots) < upper]
  breaks <- sort(unique(c(0, critical, upper)))
  midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2
  # Evaluate the derivative on every interval between its real roots.
  slopes <- rep(0, length(midpoints))
  for (i in rev(seq_along(derivative))) {
    slopes <- slopes * midpoints + derivative[i]
  }
  if (any(!is.finite(slopes)) || any(slopes <= 0)) {
    stop("The radial polynomial must be strictly increasing on ",
         "[0, max_off_axis_angle].", call. = FALSE)
  }
  boundary <- .calc_relative_radius(upper, projection_coef)
  if (!is.finite(boundary) || boundary <= 0) {
    stop("The radial polynomial has no finite positive radius at the limit.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @rdname internal-validation
#' @export
.check_internal_anisotropic_parameters <- function(parameters) {
  required <- c("off_axis_coef", "sin_internal_azimuth_coef",
                "cos_internal_azimuth_coef", "atan_scale")
  if (!is.list(parameters) || length(parameters) != length(required) ||
      is.null(names(parameters)) || anyDuplicated(names(parameters)) ||
      !setequal(names(parameters), required)) {
    stop("`parameters` must contain exactly: ", paste(required, collapse = ", "),
         ".", call. = FALSE)
  }
  for (name in required) {
    .check_vector(parameters[[name]], "numeric",
                  length = switch(name, off_axis_coef = NULL,
                                  atan_scale = 1, 2),
                  sign = if (name == "atan_scale") "positive" else "any",
                  name = paste0("parameters$", name))
  }
  if (!length(parameters$off_axis_coef)) {
    stop("`parameters$off_axis_coef` cannot be empty.", call. = FALSE)
  }
  invisible(TRUE)
}
