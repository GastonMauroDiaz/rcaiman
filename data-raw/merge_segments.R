#' Multiresolution region merging
#'
#' @description
#' Iteratively merge adjacent image regions according to a weighted
#' radiometric heterogeneity criterion inspired by the multiresolution
#' segmentation strategy popularized by eCognition.
#'
#' @details
#' This function implements a bottom-up region-merging procedure operating on
#' an existing labeled segmentation raster.
#'
#' The algorithm starts from an initial oversegmentation and iteratively merges
#' neighboring regions that minimize the heterogeneity of the resulting merged
#' object.
#'
#' Unlike classical raster-based implementations where each pixel contributes
#' equally, this implementation supports weighted radiometric statistics through
#' an area-weight image (`w`). This enables segmentation directly in distorted
#' hemispherical projections while preserving consistency in spherical space.
#'
#' The intended workflow is:
#'
#' \enumerate{
#' \item generate an initial oversegmentation;
#' \item compute adjacency relationships among neighboring regions;
#' \item iteratively identify mutually best-fitting neighboring regions;
#' \item execute all admissible non-conflicting merges in parallel within each
#' iteration;
#' \item stop when no admissible merge satisfies the scale criterion.
#' }

#' # Multiresolution merging strategy
#'
#' The algorithm follows the general principles of classical multiresolution
#' segmentation:
#'
#' \itemize{
#' \item segmentation starts from small primitive objects;
#' \item neighboring regions are evaluated iteratively;
#' \item each region selects its best-fitting neighbor;
#' \item merges are accepted only when the relationship is mutual;
#' \item multiple independent merges may occur during the same iteration;
#' \item the process stops naturally as objects become increasingly
#' heterogeneous.
#' }
#'
#' More precisely:
#'
#' \enumerate{
#' \item for every region, all adjacent neighbors are evaluated;
#'
#' \item the neighboring region producing the minimum merged heterogeneity is
#' selected;
#'
#' \item a merge candidate is accepted only if:
#'
#' \deqn{
#' A \rightarrow B
#' \quad\text{and}\quad
#' B \rightarrow A
#' }
#'
#' meaning both regions select each other as their best-fitting neighbor;
#'
#' \item all admissible mutual pairs are sorted by increasing heterogeneity;
#'
#' \item merges are executed greedily while preventing conflicts, meaning that
#' a region can participate in at most one merge per iteration;
#'
#' \item region statistics and adjacency relationships are updated
#' incrementally after each accepted merge;
#'
#' \item the process repeats until no admissible mutual merge satisfies the
#' scale criterion.
#' }
#'
#' # Computational notes
#'
#' Region statistics are updated incrementally after each merge.
#'
#' Unlike naive implementations, the algorithm does not recompute region
#' statistics from raw pixels after every iteration.
#'
#' Adjacency relationships are also updated incrementally through a dynamic
#' neighborhood graph representation.
#'
#' Multiple non-conflicting merges are executed during each iteration,
#' substantially reducing the total number of global iterations required for
#' convergence compared to purely sequential single-merge strategies.
#'
#' Region membership is also updated incrementally using cached pixel indices,
#' avoiding repeated full-image relabeling operations.
#'
#' The implementation therefore combines conceptual transparency with
#' substantially improved computational efficiency relative to naive
#' bottom-up region-merging implementations.
#'
#' @param r A single-layer `SpatRaster` containing radiometric values used for
#' segmentation.
#'
#' @param w A single-layer `SpatRaster` containing per-pixel area weights,
#' typically generated with \code{\link{area_weight_image}}.
#'
#' @param seg A single-layer labeled segmentation raster.
#'
#' @param scale_parameter Numeric scalar controlling the maximum admissible
#' merge heterogeneity relative to the initial segmentation scale.
#'
#' @param ignore_label_0 Logical. If `TRUE`, pixels labeled `0` are ignored and
#' treated as missing.
#'
#' @param max_iterations Maximum number of merge iterations.
#'
#' @param verbose Logical. If `TRUE`, progress messages are printed during
#' execution.
#'
#' @return
#' A `SpatRaster` containing the merged segmentation.
#'
#' @seealso
#' \code{\link{compute_adjacency}}
#'
#' \code{\link{area_weight_image}}
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' # Initial oversegmentation
#' seg0 <- healpix_segmentation(z, a, 32)
#'
#' # Spherical pixel weights
#' w <- area_weight_image(z)
#'
#' # Multiresolution merging
#' seg <- merge_segments(
#'   r = r,
#'   w = w,
#'   seg = seg0,
#'   scale_parameter = 10
#' )
#'
#' plot(seg)
#' plot(extract_feature(r, seg))
#'
#' }
#'
#' @export
merge_segments <- function(r,
                           w,
                           seg,
                           scale_parameter,
                           ignore_label_0 = TRUE,
                           max_iterations = Inf,
                           verbose = FALSE) {

  .assert_single_layer(r)
  .assert_single_layer(seg)
  .assert_area_weight_image(w)
  .assert_same_geom(r, w)
  .assert_same_geom(r, seg)

  seg <- terra::deepcopy(seg)

  # ============================================================
  # Extract matrices
  # ============================================================

  x_r <- terra::as.matrix(r, wide = TRUE)
  x_seg <- terra::as.matrix(seg, wide = TRUE)
  x_w <- terra::as.matrix(w, wide = TRUE)

  if (ignore_label_0) {
    x_seg[x_seg == 0] <- NA
  }

  # ============================================================
  # Initialize region statistics
  # ============================================================

  ids <- sort(unique(as.vector(x_seg)))
  ids <- ids[!is.na(ids)]

  max_id <- max(ids)

  region_weight <- numeric(max_id)
  region_mean <- numeric(max_id)
  region_ssd <- numeric(max_id)

  active <- logical(max_id)

  region_pixels <- vector("list", length(ids))
  names(region_pixels) <- as.character(ids)

  for (id in ids) {

    idx <- which(x_seg == id)

    region_pixels[[as.character(id)]] <- idx

    values <- x_r[idx]
    weights <- x_w[idx]

    ok <- !is.na(values) & !is.na(weights)

    values <- values[ok]
    weights <- weights[ok]

    weight <- sum(weights)

    mean_value <- sum(weights * values) / weight

    ssd <- sum(weights * (values - mean_value)^2)

    region_weight[id] <- weight
    region_mean[id] <- mean_value
    region_ssd[id] <- ssd

    active[id] <- TRUE
  }

  # ============================================================
  # Initial adjacency
  # ============================================================

  adj <- compute_adjacency(
    seg,
    directed = FALSE,
    ignore_label_0 = ignore_label_0,
    return_border_length = FALSE
  )

  neighbors <- adjacency_to_list(adj)

  # ============================================================
  # Normalize scale parameter
  # ============================================================

  initial_h <- numeric(nrow(adj))

  for (i in seq_len(nrow(adj))) {

    id1 <- adj$from[i]
    id2 <- adj$to[i]

    w1 <- region_weight[id1]
    w2 <- region_weight[id2]

    mu1 <- region_mean[id1]
    mu2 <- region_mean[id2]

    ssd1 <- region_ssd[id1]
    ssd2 <- region_ssd[id2]

    initial_h[i] <- (
      ssd1 +
        ssd2 +
        ((w1 * w2) / (w1 + w2)) * (mu1 - mu2)^2
    )
  }

  scale_unit <- stats::median(initial_h, na.rm = TRUE)

  scale_parameter_internal <- scale_parameter * scale_unit

  if (verbose) {

    message(
      paste(
        "Scale normalization unit:",
        round(scale_unit, 6)
      )
    )
  }

  # ============================================================
  # Main loop
  # ============================================================

  iteration <- 1

  repeat {

    if (iteration > max_iterations) {
      break
    }

    ids <- which(active)

    if (length(ids) <= 1) {
      break
    }

    # ------------------------------------------------------------
    # Find best neighbor for each region
    # ------------------------------------------------------------

    best_to <- integer(max_id)
    best_h <- rep(Inf, max_id)

    for (id in ids) {

      key <- as.character(id)

      neighs <- neighbors[[key]]

      if (is.null(neighs) || length(neighs) == 0) {
        next
      }

      local_best_h <- Inf
      local_best_to <- NA_integer_

      w1 <- region_weight[id]
      mu1 <- region_mean[id]
      ssd1 <- region_ssd[id]

      for (nb in neighs) {

        w2 <- region_weight[nb]
        mu2 <- region_mean[nb]
        ssd2 <- region_ssd[nb]

        h <- ssd1 + ssd2 +
          ((w1 * w2) / (w1 + w2)) * (mu1 - mu2)^2

        if (!is.finite(h)) {
          next
        }

        if (h < local_best_h) {
          local_best_h <- h
          local_best_to <- nb
        }
      }

      best_to[id] <- local_best_to
      best_h[id] <- local_best_h
    }

    # ------------------------------------------------------------
    # Find mutual best-fitting pairs
    # ------------------------------------------------------------

    candidates <- list()
    k <- 1

    for (id in ids) {

      key <- as.character(id)

      nb <- best_to[id]

      if (is.na(nb) || nb == 0L) {
        next
      }

      reverse_to <- best_to[nb]

      if (is.na(reverse_to) || reverse_to != id) {
        next
      }

      h <- best_h[id]

      if (h > scale_parameter_internal) {
        next
      }

      pair <- sort(c(id, nb))

      candidates[[k]] <- data.frame(
        from = pair[1],
        to = pair[2],
        h = h
      )

      k <- k + 1
    }

    if (length(candidates) == 0) {

      if (verbose) {
        message("No admissible mutual merges remaining.")
      }

      break
    }

    candidates <- do.call(rbind, candidates)

    candidates <- unique(candidates)


    # ------------------------------------------------------------
    # Execute all non-conflicting mutual merges
    # ------------------------------------------------------------

    candidates <- candidates[order(candidates$h), ]

    used <- logical(max(ids))

    merge_from <- integer(nrow(candidates))
    merge_to <- integer(nrow(candidates))

    n_merges <- 0L

    for (i in seq_len(nrow(candidates))) {

      from <- candidates$from[i]
      to <- candidates$to[i]

      if (used[from] || used[to]) {
        next
      }

      n_merges <- n_merges + 1L

      merge_from[n_merges] <- from
      merge_to[n_merges] <- to

      used[from] <- TRUE
      used[to] <- TRUE
    }

    if (n_merges == 0L) {

      if (verbose) {
        message("No admissible non-conflicting merges.")
      }

      break
    }

    if (verbose) {

      message(
        paste(
          "Iteration:", iteration,
          "| merges:", n_merges
        )
      )
    }

    # ------------------------------------------------------------
    # Apply merges
    # ------------------------------------------------------------

    for (i in seq_len(n_merges)) {

      from <- merge_from[i]
      to <- merge_to[i]

      w1 <- region_weight[from]
      w2 <- region_weight[to]

      mu1 <- region_mean[from]
      mu2 <- region_mean[to]

      ssd1 <- region_ssd[from]
      ssd2 <- region_ssd[to]

      weight <- w1 + w2

      mean_value <- (w1 * mu1 + w2 * mu2) / weight

      ssd <- ssd1 + ssd2 +
        ((w1 * w2) / weight) * (mu1 - mu2)^2

      region_weight[from] <- weight
      region_mean[from] <- mean_value
      region_ssd[from] <- ssd

      active[to] <- FALSE

      # ----------------------------------------------------------
      # Relabel segmentation
      # ----------------------------------------------------------

      idx_to <- region_pixels[[as.character(to)]]

      x_seg[idx_to] <- from

      region_pixels[[as.character(from)]] <- c(
        region_pixels[[as.character(from)]],
        idx_to
      )

      region_pixels[[as.character(to)]] <- NULL

      # ----------------------------------------------------------
      # Update adjacency
      # ----------------------------------------------------------

      neigh_a <- neighbors[[as.character(from)]]
      neigh_b <- neighbors[[as.character(to)]]

      new_neighbors <- union(neigh_a, neigh_b)

      new_neighbors <- setdiff(
        new_neighbors,
        c(from, to)
      )

      neighbors[[as.character(from)]] <- new_neighbors

      neighbors[[as.character(to)]] <- NULL

      for (nb in new_neighbors) {

        nb_key <- as.character(nb)

        current <- neighbors[[nb_key]]

        current[current == to] <- from

        current <- unique(current)

        current <- current[current != nb]

        neighbors[[nb_key]] <- current
      }
    }

    iteration <- iteration + 1
  }

  seg[] <- as.vector(t(x_seg))

  seg
}

# compute_union_heterogeneity <- function(a, b) {
#
#   w1 <- a$weight
#   w2 <- b$weight
#
#   mu1 <- a$mean
#   mu2 <- b$mean
#
#   ssd1 <- a$ssd
#   ssd2 <- b$ssd
#
#   ssd_union <- ssd1 + ssd2 + ((w1 * w2) / (w1 + w2)) * (mu1 - mu2)^2
#
#   ssd_union
# }
#
# merge_region_stats <- function(a, b) {
#
#   w1 <- a$weight
#   w2 <- b$weight
#
#   mu1 <- a$mean
#   mu2 <- b$mean
#
#   ssd1 <- a$ssd
#   ssd2 <- b$ssd
#
#   weight <- w1 + w2
#
#   mean_value <- (w1 * mu1 + w2 * mu2) / weight
#
#   ssd <- ssd1 + ssd2 + ((w1 * w2) / weight) * (mu1 - mu2)^2
#
#   list(
#     weight = weight,
#     mean = mean_value,
#     ssd = ssd
#   )
# }

adjacency_to_list <- function(adj) {

  neighbors <- list()

  if (nrow(adj) == 0) {
    return(neighbors)
  }

  for (i in seq_len(nrow(adj))) {

    a <- adj$from[i]
    b <- adj$to[i]

    a_key <- as.character(a)
    b_key <- as.character(b)

    neighbors[[a_key]] <-
      unique(c(neighbors[[a_key]], b))

    neighbors[[b_key]] <-
      unique(c(neighbors[[b_key]], a))
  }

  neighbors
}

compute_adjacency <- function(seg,
                              directed = FALSE,
                              ignore_label_0 = TRUE,
                              return_border_length = TRUE) {

  .assert_single_layer(seg)

  x <- terra::as.matrix(seg, wide = TRUE)

  if (ignore_label_0) {
    x[x == 0] <- NA
  }

  nr <- nrow(x)
  nc <- ncol(x)

  edges <- vector("list", nr * nc * 2)

  k <- 1

  # ============================================================
  # Horizontal neighbors
  # ============================================================

  for (i in 1:nr) {

    for (j in 1:(nc - 1)) {

      a <- x[i, j]
      b <- x[i, j + 1]

      if (!is.na(a) &&
          !is.na(b) &&
          a != b) {

        if (!directed) {
          pair <- sort(c(a, b))
        } else {
          pair <- c(a, b)
        }

        edges[[k]] <- pair
        k <- k + 1
      }
    }
  }

  # ============================================================
  # Vertical neighbors
  # ============================================================

  for (i in 1:(nr - 1)) {

    for (j in 1:nc) {

      a <- x[i, j]
      b <- x[i + 1, j]

      if (!is.na(a) &&
          !is.na(b) &&
          a != b) {

        if (!directed) {
          pair <- sort(c(a, b))
        } else {
          pair <- c(a, b)
        }

        edges[[k]] <- pair
        k <- k + 1
      }
    }
  }

  edges <- edges[1:(k - 1)]

  if (length(edges) == 0) {

    if (return_border_length) {

      return(
        data.frame(
          from = integer(0),
          to = integer(0),
          border = integer(0)
        )
      )

    } else {

      return(
        data.frame(
          from = integer(0),
          to = integer(0)
        )
      )
    }
  }

  edges <- do.call(rbind, edges)

  if (is.null(dim(edges))) {

    if (return_border_length) {

      return(
        data.frame(
          from = integer(0),
          to = integer(0),
          border = integer(0)
        )
      )

    } else {

      return(
        data.frame(
          from = integer(0),
          to = integer(0)
        )
      )
    }
  }

  colnames(edges) <- c("from", "to")

  edges <- as.data.frame(edges)

  if (!return_border_length) {

    edges <- unique(edges)

    rownames(edges) <- NULL

    return(edges)
  }

  # ============================================================
  # Count shared border length
  # ============================================================

  edges$border <- 1L

  edges <- stats::aggregate(
    border ~ from + to,
    data = edges,
    FUN = sum
  )

  rownames(edges) <- NULL

  edges
}
