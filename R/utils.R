radian2degree <- function(x) x * 180 / pi
degree2radian <- function(x) x * pi / 180

.get_max <- function(r) max(r[], na.rm = TRUE)
.get_min <- function(r) min(r[], na.rm = TRUE)

.is_even <- function(x) round(x/2) == x/2

.is_whole <- function(x) round(x) == x

.is_integerish <- function(x) x == round(x)

.check_if_r_was_normalized <- function(r) {
  if (max(r[], na.rm = TRUE) > 1)
    warning("Please check if \"r\" was correctly normalized")
}

.check_if_r_z_and_a_are_ok <- function(r, z, a) {
  stopifnot(class(r) == "RasterLayer")
  .check_if_r_was_normalized(r)
  stopifnot(class(z) == "RasterLayer")
  stopifnot(class(a) == "RasterLayer")
  stopifnot(.get_max(z) < 90)
  compareRaster(r, z)
  compareRaster(z, a)
}

.this_requires_EBImage <- function() {
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop(paste("Package \"EBImage\" needed for this function to work.",
               "Please install it. Find instructions here:",
               "https://bioconductor.org/packages/release/bioc/html/EBImage.html"
               ),
         call. = FALSE)
  }
}

.calc_angular_distance <- function(z1, a1, z2, a2) {
  #https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  acos(pmax(pmin(cos(z1) * cos(z2) + sin(z1) * sin(z2) * cos(abs(a2 - a1)), 1), -1))
}

.calc_rmse <- function(x) sqrt(mean(x^2))

.calc_r2 <- function(x, y) {
  model <- lm(y ~ x)
  summary(model)$r.squared
}

# .estimate_w <- function(bin) {
#   sky_border <- focal(bin, matrix(c(0, 1, 0,
#                                     1,-4, 1,
#                                     0, 1, 0), nrow = 3))
#   no_of_px_on_circle <- sum(!is.na(bin)[])
#   no_of_px_on_sky_border <- sum((sky_border != 0)[], na.rm = TRUE)
#   w <- 0.5 + no_of_px_on_sky_border / no_of_px_on_circle
#   if (w > 0.9) w <- 0.9
#   w
# }

.decode_label <- function(label) {
  sector_ID <- trunc(label / 1000)
  rings_ID <- label - sector_ID * 1000
  data.frame(sector_ID, rings_ID)
}

.make_fake_las <- function(X, Y, Z){
  data_template_names <- c("X", "Y", "Z", "gpstime","Intensity", "ReturnNumber", "NumberOfReturns",
                           "ScanDirectionFlag", "EdgeOfFlightline",  "Classification",
                           "Synthetic_flag", "Keypoint_flag", "Withheld_flag", "ScanAngleRank", "UserData",
                           "PointSourceID",  "R",  "G",  "B")
  data_template <- matrix(ncol = length(data_template_names), nrow = length(X))
  data_template <- data.frame(data_template)
  names(data_template) <- data_template_names
  data_template[] <- 0

  fake_las <- new("LAS")
  fake_las@data <- as(data_template, "data.table")
  fake_las@data$X <- X
  fake_las@data$Y <- Y
  fake_las@data$Z <- Z
  fake_las
}

.makeF8single <- function(a, ...) { # single layer output

  function(x, filename = "", ...) {
    # code from raster-package vignette
    # "Writing functions for large raster files"
    # function referred as f8
    out <- raster(x)
    big <- ! canProcessInMemory(out, 3)
    filename <- trim(filename)
    if (big & filename == '') {
      filename <- rasterTmpFile()
    }
    if (filename != '') {
      out <- writeStart(out, filename, ...)
      todisk <- TRUE
    } else {
      vv <- matrix(ncol=nrow(out), nrow=ncol(out))
      todisk <- FALSE
    }

    bs <- blockSize(x)
    pb <- pbCreate(bs$n, ...)

    if (todisk) {
      for (i in 1:bs$n) {
        v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
        v <- a(v, ...) # new code
        out <- writeValues(out, v, bs$row[i])
        pbStep(pb, i)
      }
      out <- writeStop(out)
    } else {
      for (i in 1:bs$n) {
        v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
        v <- a(v, ...) # new code
        cols <- bs$row[i]:(bs$row[i]+bs$nrows[i]-1)
        vv[,cols] <- matrix(v, nrow=out@ncols)
        pbStep(pb, i)
      }
      out <- setValues(out, as.vector(vv))
    }
    pbClose(pb)
    return(out)
  }
}

.makeF8multi <- function(a, ...) { # multi layer output

  function(x, filename = "", ...) {
    # based on code from raster-package vignette
    # "Writing functions for large raster files"
    # function referred as f8
    out <- brick(x) #
    big <- ! canProcessInMemory(out, 3)
    filename <- trim(filename)
    if (big & filename == '') {
      filename <- rasterTmpFile()
    }
    if (filename != '') {
      out <- writeStart(out, filename, ...)
      todisk <- TRUE
    } else {
      vv <- array(dim = c(ncol(out), nrow(out), nlayers(out)))
      todisk <- FALSE
    }

    bs <- blockSize(x)
    pb <- pbCreate(bs$n, ...)

    if (todisk) {
      for (i in 1:bs$n) {
        v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
        v <- a(v, ...)
        out <- writeValues(out, v, bs$row[i])
        pbStep(pb, i)
      }
      out <- writeStop(out)
    } else {
      for (i in 1:bs$n) {
        v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
        v <- a(v, ...)
        cols <- bs$row[i]:(bs$row[i]+bs$nrows[i]-1)
        vv[,cols,] <- array(v, dim=c(bs$nrows[i],nrow=out@ncols, nlayers(x)))
        pbStep(pb, i)
      }
      out <- setValues(out, as.vector(vv))
    }
    pbClose(pb)
    return(out)
  }
}


