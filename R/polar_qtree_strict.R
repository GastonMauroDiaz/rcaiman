#' Polar Quadtree Segmentation with Optional Parallelism
#'
#' This function segments a raster `r` in polar coordinates into circular trapezoid cells
#' (defined by zenith and azimuth angle bands). Segments are recursively subdivided
#' if splitting reduces the overall standard deviation of `r`-values by more than
#' `scale_parameter`. Symmetric angular resolution is maintained (each split divides angles equally).
#'
#' @param r A `SpatRaster` of data values to segment.
#' @param z A `SpatRaster` of the same dimension giving the zenith angle for each cell.
#' @param a A `SpatRaster` of the same dimension giving the azimuth angle for each cell (in degrees).
#' @param scale_parameter Numeric threshold; a segment is split only if
#'        (sd_parent - sum(sd_children)) > `scale_parameter`.
#' @param angle_width Base angular resolution (in degrees) for zenith and azimuth bands.
#' @param parallel Logical. If TRUE, perform segmentation in parallel. Requires a `foreach` backend (e.g. via `doParallel`).
#' @param diagnose Logical. To explore a viable value for the scale parameter
#'
#' @return A `SpatRaster` of the same size as `r`, where each cell’s value is an integer ID of its segment.
#'
#' @details The raster is first converted to plain vectors (`terra::values`) of cell values and angles.
#'          Segments are represented by (z_min,z_max,azim_min,azim_max) and split into four by halving angles.
#'          Standard deviation is used as the split criterion. The final segmentation is applied
#'          via `terra::setValues` to create the output raster.
#' @export
#' @examples
#' \dontrun{
#' caim <- read_caim()
#' z <- zenith_image(ncol(caim), lens())
#' a <- azimuth_image(z)
#' m <- !is.na(z)
#' com <- compute_complementary_gradients(caim/sky)
#' chroma <- max(com$blue_yellow, com$cyan_red)
#' seg <- polar_qtree_strict(chroma, z, a, m, 0.2, 30, 50, FALSE)
#' }
polar_qtree_strict <- function(r, z, a, m,
                               scale_parameter,
                               angle_width,
                               min_size_px = 0,
                               parallel = FALSE,
                               diagnose = FALSE) {
  # Basic checks
  if (!inherits(r, "SpatRaster") || !inherits(z, "SpatRaster") || !inherits(a, "SpatRaster")) {
    stop("r, z, a must be SpatRaster objects.")
  }
  if (!terra::compareGeom(r, z, stopOnError = FALSE) ||
      !terra::compareGeom(r, a, stopOnError = FALSE)) {
    stop("r, z, a must share the same extent and resolution.")
  }

  diagnose_first_level <- function(r_vals, z_vals, a_vals, angle_width) {
    valid <- which(!is.na(r_vals))
    z_max <- max(z_vals[valid], na.rm=TRUE)
    z_breaks <- seq(0, z_max + angle_width, by=angle_width)
    a_breaks <- seq(0, 360, by=angle_width)
    if (tail(a_breaks,1) < 360) a_breaks <- c(a_breaks, 360)

    out <- list()
    id <- 1L
    for (i in seq_len(length(z_breaks)-1)) {
      for (j in seq_len(length(a_breaks)-1)) {
        zmin <- z_breaks[i]; zmax <- z_breaks[i+1]
        amin <- a_breaks[j]; amax <- a_breaks[j+1]
        sel <- valid[ z_vals[valid] >= zmin & z_vals[valid] < zmax &
                        a_vals[valid] >= amin & a_vals[valid] < amax ]
        if (length(sel) < 2) next
        sd_parent <- sd(r_vals[sel], na.rm=TRUE)
        # Child midpoints
        zmid <- (zmin+zmax)/2
        amid <- (amin+amax)/2
        children <- list(
          which(z_vals>=zmin & z_vals< zmid & a_vals>=amin & a_vals< amid),
          which(z_vals>=zmin & z_vals< zmid & a_vals>=amid & a_vals< amax),
          which(z_vals>=zmid & z_vals< zmax & a_vals>=amin & a_vals< amid),
          which(z_vals>=zmid & z_vals< zmax & a_vals>=amid & a_vals< amax)
        )
        sum_child_sd <- 0
        for (ch in children) {
          ch <- intersect(ch, sel)
          if (length(ch)>=2) sum_child_sd <- sum_child_sd + sd(r_vals[ch], na.rm=TRUE)
        }
        out[[length(out)+1]] <- data.frame(
          seg_id      = id,
          n_pixels    = length(sel),
          sd_parent   = sd_parent,
          sum_child_sd= sum_child_sd,
          delta       =  sum_child_sd - sd_parent
        )
        id <- id + 1L
      }
    }
    do.call(rbind, out)
  }

  if (diagnose) {
    r_vals <- r[m]
    z_vals <- z[m]
    a_vals <- a[m]

    df_delta <- diagnose_first_level(r_vals, z_vals, a_vals,
                                     angle_width = angle_width)
    return(df_delta)

  } else {

    # Extract vectors
    r_vals <- terra::values(r)
    z_vals <- terra::values(z)
    a_vals <- terra::values(a)
    m <- terra::values(m)

    # Normalize azimuth to [0,360)
    # a_vals[!is.na(a_vals) & a_vals >= 360] <- a_vals[!is.na(a_vals) & a_vals >= 360] - 360

    # Check parallel backend if needed
    if (parallel) {
      if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop(
        "Package \"doParallel\" needed for this function to work. Please
        install it. Additionally, this package must be manually
        loaded for this function to work in parallel, so please add
        `require(doParallel)` or `library(doParallel)` to the script."
        ,
        call. = FALSE)
      }
    }
    if (parallel & !any(grepl("doParallel", search()))) {
      stop(
        "Package \"doParallel\" must be manually loaded for this function
        to work, so please add `require(doParallel)` or `library(doParallel)`
        to the script."
        ,
        call. = FALSE)
    }
    # if (parallel && !foreach::getDoParRegistered()) {
    #   stop("No parallel backend registered. Please run doParallel::registerDoParallel().")
    # }

    # Build initial angle breaks
    z_max <- max(z_vals[m], na.rm = TRUE)
    z_breaks <- seq(0, z_max + angle_width, by = angle_width)
    if (z_breaks[length(z_breaks)] < z_max) {
      z_breaks <- c(z_breaks, z_max)
    }
    a_breaks <- seq(0, 360, by = angle_width)
    if (a_breaks[length(a_breaks)] < 360) {
      a_breaks <- c(a_breaks, 360)
    }

    # List of base segments: (zmin,zmax,amin,amax)
    segments_init <- vector("list", (length(z_breaks) - 1) * (length(a_breaks) - 1))
    idx_seg <- 1L
    for (i in seq_len(length(z_breaks) - 1)) {
      for (j in seq_len(length(a_breaks) - 1)) {
        segments_init[[idx_seg]] <- c(
          z_breaks[i], z_breaks[i+1],
          a_breaks[j], a_breaks[j+1]
        )
        idx_seg <- idx_seg + 1L
      }
    }

    # Recursive subdivision on trapezoidal segment
    subdivide_one <- function(zmin, zmax, amin, amax) {
      # select only valid cells inside this angular trapezoid
      sel <- which(
        m &
          z_vals >= zmin & z_vals < zmax &
          a_vals >= amin & a_vals < amax
      )
      if (length(sel) < 2) {
        return(list(c(zmin, zmax, amin, amax)))
      }
      vals <- r_vals[sel]
      sd_parent <- stats::sd(vals, na.rm = TRUE)
      if (is.na(sd_parent) || sd_parent == 0) {
        return(list(c(zmin, zmax, amin, amax)))
      }
      # mid angles
      zmid <- (zmin + zmax) / 2
      amid <- (amin + amax) / 2
      # child segments
      children <- list(
        c(zmin, zmid, amin, amid),
        c(zmin, zmid, amid, amax),
        c(zmid, zmax, amin, amid),
        c(zmid, zmax, amid, amax)
      )
      # sum of child SD
      sum_child_sd <- 0
      for (seg in children) {
        c_sel <- which(
          m &
            z_vals >= seg[1] & z_vals < seg[2] &
            a_vals >= seg[3] & a_vals < seg[4]
        )
        if (length(c_sel) >= 2) {
          sum_child_sd <- sum_child_sd + stats::sd(r_vals[c_sel], na.rm = TRUE)
        }
      }
      # decide
      if (all((sum_child_sd - sd_parent) > scale_parameter,
              length(vals) > min_size_px)) {
        out <- list()
        for (seg in children) {
          out <- c(out, subdivide_one(seg[1], seg[2], seg[3], seg[4]))
        }
        return(out)
      } else {
        return(list(c(zmin, zmax, amin, amax)))
      }
    }

    # Process base segments
    if (parallel) {
      # Configurar paralelismo (foreach + doParallel)
      cores <- parallel::detectCores(logical=FALSE)
      if (is.na(cores) || cores < 1) cores <- 1
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)

      final_list <- foreach::foreach(
        seg = segments_init, .combine = c, .multicombine = TRUE,
        # .export = c("r_vals","z_vals","a_vals","m","scale_parameter"),
        .packages = "stats"
      ) %dopar% {
        subdivide_one(seg[1], seg[2], seg[3], seg[4])
      }
      parallel::stopCluster(cl)
    } else {
      final_list <- list()
      for (seg in segments_init) {
        final_list <- c(final_list, subdivide_one(seg[1], seg[2], seg[3], seg[4]))
      }
    }

    # Build output ID vector
    n_cells <- length(r_vals)
    seg_ids <- integer(n_cells)
    id <- 1L
    for (seg in final_list) {
      sel <- which(
        m &
          z_vals >= seg[1] & z_vals < seg[2] &
          a_vals >= seg[3] & a_vals < seg[4]
      )
      seg_ids[sel] <- id
      id <- id + 1L
    }

    # Return SpatRaster with segment IDs
    out_r <- terra::setValues(r, seg_ids)
    names(out_r) <- "Polar quad-tree"
    return(out_r)
  }
}
