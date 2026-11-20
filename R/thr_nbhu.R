#' Nobis-Hunziker automatic thresholding
#'
#' @description
#' Compute a binary threshold using the edge-based criterion
#' proposed by \insertCite{Nobis2005;textual}{rcaiman}.
#'
#' @details
#' The method proposed by \insertCite{Nobis2005;textual}{rcaiman}
#' searches for the threshold maximizing the mean radiometric contrast across
#' binary class boundaries. For each threshold \eqn{t}, the image composed of
#' digital numbers (\eqn{\mathrm{DN}}) is binarized according to:
#'
#' \deqn{
#' f(b, t) =
#' \left\{
#' \begin{array}{ll}
#' 1 & : \mathrm{DN} > t \\
#' 0 & : \mathrm{DN} \le t
#' \end{array}
#' \right.
#' }
#'
#' A moving 2 \eqn{\times} 2 window is evaluated over the image:
#'
#' \preformatted{
#' p11  p12
#' p21  p22
#' }
#'
#' For each window, the following pixel pairs are tested:
#'
#' \itemize{
#'   \item `p11` vs `p12`
#'   \item `p11` vs `p21`
#'   \item `p11` vs `p22`
#'   \item `p12` vs `p21`
#' }
#'
#' Only pairs with opposite binary classifications contribute to the score.
#'
#' The underlying C++ implementation computes both the contrast score and the
#' number of contributing pixel pairs for all thresholds between
#' `min_thresh` and `max_thresh`.
#'
#' The original method selects the final threshold maximizing the contrast score.
#' This assumes dark objects against a relatively uniform bright background.
#' Under hemispherical diffuse-light conditions, this often favors thresholds
#' driven by the brightest sky region rather than by the overall sky-canopy
#' separation.
#'
#' This implementation optionally penalizes thresholds supported by few
#' pixel-pair comparisons:
#'
#' \deqn{
#' J(t) = S(t)
#' \cdot
#' \left(
#' \frac{n(t)}{\max(n)}
#' \right)^\alpha
#' }
#'
#' where \eqn{S(t)} is the original contrast score, \eqn{n(t)} is the number
#' of contributing pixel pairs, and \eqn{\alpha} is controlled by
#' `support_weight`.
#'
#' @param r numeric [terra::SpatRaster-class] of one layer.
#' @param min_thresh,max_thresh numeric vector of length one interpreted as
#'   positive integer. Minimum and maximum thresholds evaluated during
#'   exhaustive search.
#' @param support_weight numeric vector of length one. Must be non-negative.
#'   Controls the influence of the number of contributing pixel pairs in the
#'   modified score. Use `0` (default) to reproduce the original criterion.
#' @param smooth_curve logical vector of length one. If `TRUE`, apply a median filter to the
#'   threshold-score curve before selecting the optimum threshold. This
#'   reproduces the smoothing step used in the MATLAB implementation by
#'   [Korhonen et al](https://www.mathworks.com/matlabcentral/fileexchange/24314-matlab-codes-for-canopy-image-analysis).
#' @param smooth_width numeric vector of length one. Must be positive,
#'   interpreted as integer, and odd. Window size of the median filter applied
#'   when `smooth_curve = TRUE`.
#' @param return character vector of length one. Either `"thr"` to return the
#'   optimal threshold or `"all"` to return the complete threshold-score
#'   curve.
#'
#' @return
#' If `return = "thr"`, an integer of length one.
#'
#' If `return = "all"`, a list with:
#'
#' \describe{
#'   \item{`thr`}{
#'   Integer vector of evaluated thresholds.
#'   }
#'
#'   \item{`score`}{
#'   Mean edge-contrast score associated with each threshold.
#'   }
#'
#'   \item{`n`}{
#'   Number of pixel-pair comparisons with opposite binary classifications
#'   contributing to each score.
#'   }
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#' @examples
#' path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
#' r <- read_caim(path, fill_na = "none")$Blue
#' r <- crop_caim(r, c(600, 800), 200, 200)
#'
#' thr <- thr_nbhu(r)
#' plot(binarize_with_thr(r, thr))
#' # See apply_by_direction() for thresholds in angular windows.
thr_nbhu <- function(
    r,
    min_thresh = 1,
    max_thresh = 254,
    support_weight = 0,
    smooth_curve = FALSE,
    smooth_width = 7,
    return = "thr"
) {

  .assert_single_layer(r, allow_na = FALSE)
  .assert_choice(return, c("thr", "all"))
  .check_vector(min_thresh, "integerish", sign = "positive")
  .check_vector(max_thresh, "integerish", sign = "positive")
  .check_vector(support_weight, "numeric", sign = "nonnegative")
  .check_vector(smooth_curve, "logical")
  .check_vector(smooth_width, "odd_integerish", sign = "positive")

  img <- matrix(
    terra::values(r, mat = FALSE),
    nrow = terra::nrow(r),
    ncol = terra::ncol(r),
    byrow = TRUE
  )

  out <- thr_nbhu_cpp(
    img = img,
    min_thresh = min_thresh,
    max_thresh = max_thresh
  )

  out$score <- out$score  * (out$n/max(out$n))^support_weight


  if (smooth_curve) {
    out$score <- stats::runmed(
      x = out$score,
      k = smooth_width,
      endrule = "median"
    )
  }

  if (return == "all") {
    return(out)
  }

  out$thr[which.max(out$score)]
}
