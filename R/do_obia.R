.sRGB2LAB <- function(x) {
  lab_colors <- colorspace::sRGB(x$Red[], x$Green[], x$Blue[]) %>%
    as(., "LAB") %>%
    colorspace::coords(.)
  x$Red <- lab_colors[,1]
  x$Green <- lab_colors[,2]
  x$Red <- lab_colors[,3]
  names(x) <- c("L", "A", "B")
  x
}

.my_isodata <- function(x)
{
  if (length(x) <= 1) stop("length(x) must be greater than 1.")
  if (stats::sd(x, na.rm = TRUE) == 0) stop("sd(x) must be greater than 0.")
  thr <- mean(x, na.rm = TRUE)
  thr.back <- 0
  while (thr != thr.back) {
    thr.back <- thr
    x0 <- x[x <= thr]
    x1 <- x[x > thr]
    thr <- (mean(x0, na.rm = TRUE) + mean(x1, na.rm = TRUE)) / 2
  }
  thr
}

.automatic_selection_of_samples <- function(blueValues, index, sampleSize) {
  blueBrightness <- tapply(blueValues, index, mean)
  rows <- as.numeric(names(blueBrightness))
  rows <- rows[order(blueBrightness)]

  if (length(rows) < (sampleSize * 5)) {
    stop("You need to decrease the scaleParameter used to produce the argument seg.")
  }


  center <- round(length(rows) / 2)
  rows.mix <- rows[seq(
    from = center - sampleSize / 2,
    to = center + sampleSize / 2, by = 1
  )]
  rows.sky <- rows[seq(
    from = length(rows) - sampleSize,
    to = length(rows), by = 1
  )]
  rows.test <- rows[is.na(match(rows, c(rows.mix, rows.sky)))]
  list(mix = rows.mix, sky = rows.sky, test = rows.test)
}

.sample_based_multiTexture_classification <-
  function(segmentation, blue, red, samples, k) {
    features <- cbind(
      EBImage::computeFeatures.haralick(segmentation, blue),
      EBImage::computeFeatures.haralick(segmentation, red)
    )
    train <- rbind(features[samples$sky, ], features[samples$mix, ])
    cl <- c(rep("sky", length(samples$sky)), rep("mix", length(samples$mix)))
    test <- features[samples$test, ]
    pred.knn <- class::knn(train, test, cl, k = k)
    list(pred.knn = pred.knn, cl = cl)
  }


#' @title Do an Object-based image analysis to classify gaps.
#'
#' @description Do an Object-based image analysis with the aim of classify
#'   gaps in full-color-upward-looking hemispherical photographs.
#'
#' @param x \code{\linkS4class{CanopyPhoto}}.
#' @param bin \code{\linkS4class{BinImage}}. The standard is a call to
#'   \code{\link{enhanceHP}} followed by a call to \code{\link{autoThr}}
#' @param z \code{\linkS4class{ZenithImage}}.
#' @param seg \code{\linkS4class{PolarSegmentation}}.
#' @param g1 \code{\linkS4class{PolarSegmentation}}. The default option is a
#'   PolarSegmentation created by calling \code{makePolarGrid(z)}. To save time in
#'   batch processing of photos token with the same camera, you can compute
#'   \code{makePolarGrid(z)} only once and provide the result through this argument.
#' @param sampleSize integer. Default is \code{50}, see Details.
#' @param k integer. Default is \code{1} nearest neighbor, see Details.
#' @param zlim \code{\linkS4class{Angle}}. Defaults are \code{30} to \code{60} degrees of
#'   zenith angle, see Details.
#' @param calibration logical. Default is \code{FALSE}, see Details.
#'
#' @details This algorithm uses object-based image analysis (OBIA). The class
#'   \emph{Gap-candidate} is assigned to pixels that are white in \code{bin} and
#'   the class \emph{Plant} to the rest of the pixels. Next, the algorithm uses
#'   this result and \code{g1} to isolate hemisphere segments with \code{1}
#'   degree of resolution that are not fully cover by \emph{Plant} (i.e., Gap
#'   Fraction > 0), which are classified as \emph{Mix-OR-Gap}. Next, the
#'   algorithm get a binary mask from this result and intersect it with the
#'   argument \code{seg}. At this point, the algorithms achieve the
#'   identification of all segments of \code{seg} that could have some gaps at
#'   pixel level (i.e., \emph{Mix-OR-Gap}). Next, the algorithm classified all
#'   this segments in \emph{Gap} or \emph{Mix} in a two stage process: (1)
#'   automatic selection of samples and (2) sample-based classification. The
#'   argument \code{sampleSize} controls the sample size for both targeted
#'   classes. The algorithm uses the brightness of the blue channel to select
#'   the samples. It assumes that brighter objects belong to \emph{Gap} and
#'   objects with middle brightness belong to \emph{Mix}. The argument \code{k}
#'   is for the knn used in the second stage of the sample-based classification.
#'   Processing continues on Mix-segments in order to unmix them at pixel level
#'   (see references for more details). The arguments \code{mnZ} and \code{mxZ}
#'   can be used to delimitate the range of zenith angle in which the
#'   aforementioned process is computed. In the rest of the image the result
#'   will be the same as \code{bin}.
#'
#'   If calibrate is set as \code{TRUE}, the process stops just after the
#'   sample-based classification described in the previous paragraph and returns
#'   a classification at object level of \emph{Plan}, \emph{Mix} and \emph{Gap}.
#'   This kind of output can be used to calibrate \code{sampleSize} and
#'   \code{k}.
#'
#' @return A \code{\linkS4class{BinImage}} by default. If \code{calibrate} is
#'   set to \code{TRUE}, a \code{\linkS4class{RasterLayer}}.
#' @noRd
#' @export
do_obia <- function(x,
                    bin,
                    z,
                    seg,
                    g1,
                    sampleSize = 50,
                    k = 1,
                    zlim = c(30, 60),
                    calibration = FALSE) {
  stopifnot(compareRaster(x, bin))
  stopifnot(compareRaster(x, z))
  stopifnot(compareRaster(x, seg))
  stopifnot(compareRaster(x, g1))

  mnZ <- zlim[1]
  mxZ <- zlim[2]
  stopifnot(mnZ < mxZ)

  g1[z > mxZ | z < mnZ] <- NA

  ## step 3

  toclass <- apply_thr(extract_feature(bin, g1, mean), 0)
  rm(g1)

  ## step 5

  toclass <- seg * toclass
  rm(seg)
  rm(z)
  toclass[toclass == 0] <- NA

  samples <- .automatic_selection_of_samples(x$Blue[], toclass[], sampleSize)

  toclass[is.na(toclass)] <- 0

  mtClass <- .sample_based_multiTexture_classification(
    as.matrix(toclass),
    as.matrix(x$Blue),
    as.matrix(x$Red),
    samples,
    k
  )

  df <- data.frame(
    c(samples$sky, samples$mix, samples$test),
    c(as.factor(mtClass$cl), mtClass$pred.knn)
  )
  df <- df[!duplicated(df[, 1]), ]

  rhara <- subs(toclass, df)

  rhara[rhara > 2] <- 1
  if (calibration) {
    rhara[is.na(rhara)] <- 0
    rhara[z > mxZ | z < mnZ] <- NA
    rhara[is.na(z)] <- NA
    return(rhara)
  } else {
    ## step 6 and 7
    rL <- raster(.sRGB2LAB(x), 1)
    rm(x)
    toclass[!(bin == 1 & rhara == 1 & rL <= 95)] <- NA
    toclass[toclass == 0] <- NA

    ## step 8
    ## Feature calculation
    rL <- rL / extract_feature(rL, toclass, mean)
    rL[is.infinite(rL)] <- NA
    ## Regional thresholding with standarized Lightness feature
    thr <- .my_isodata(rL[])
    rbin <- apply_thr(rL, thr)
    rm(rL)
    rbin[rbin == 1] <- NA
    bin <- cover(rbin, bin)
    rm(rbin)
    return(bin)
  }
}
