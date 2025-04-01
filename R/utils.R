.radian2degree <- function(x) x * 180 / pi
.degree2radian <- function(x) x * pi / 180

.get_max <- function(r) max(r[], na.rm = TRUE)
.get_min <- function(r) min(r[], na.rm = TRUE)

.decode_label <- function(label) {
  sector_ID <- trunc(label / 1000)
  ring_ID <- label - sector_ID * 1000
  ds <- data.frame(sector_ID, ring_ID)
  names(ds) <- c("sector_ID", "ring_ID")
  ds
}

.calc_rmse <- function(x) sqrt(mean(x^2))

.calc_spherical_distance <- function(z1, a1, z2, a2, radians = TRUE) {
  #https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  if (radians) {
    d <- acos(pmax(pmin(cos(z1) * cos(z2) +
                          sin(z1) * sin(z2) * cos(abs(a2 - a1)), 1), -1))
  } else {
    z1 <- z1 * pi/180
    a1   <- a1 * pi/180
    z2 <- z2 * pi/180
    a2   <- a2 * pi/180
    d <- acos(pmax(pmin(cos(z1) * cos(z2) +
                          sin(z1) * sin(z2) * cos(abs(a2 - a1)), 1), -1))
  }
  d
}

.extension <- function(file_name, new_extension = "tif") {
  file_name <- filenamer::as.filename(file_name)
  file_name <- filenamer::insert(file_name, ext = new_extension, replace = TRUE)
  as.character(file_name)
}

.make_fake_las <- function(X, Y, Z){
  data_template_names <- c("X", "Y", "Z", "gpstime","Intensity", "ReturnNumber",
                           "NumberOfReturns", "ScanDirectionFlag",
                           "EdgeOfFlightline",  "Classification",
                           "Synthetic_flag", "Keypoint_flag", "Withheld_flag",
                           "ScanAngleRank", "UserData", "PointSourceID",
                           "R",  "G",  "B")
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

.show_popup <- function(message) {
  # Create the main window
  win <- tcltk::tktoplevel()
  tcltk::tkwm.title(win, "Popup")

  # Create a label with the message
  label <- tcltk::tklabel(win, text = message, padx = 10,
                          pady = 10, justify = "left")
  tcltk::tkpack(label)

  # Create a button to close the window
  button <- tcltk::tkbutton(win, text = "OK",
                            command = function() tcltk::tkdestroy(win))
  tcltk::tkpack(button, padx = 10, pady = 10)

  # Set the focus on the window
  tcltk::tkfocus(win)

  # Run the event loop
  tcltk::tkwait.window(win)
}

.get_sky_cie <- function(z, a, model) {
  sky_cie <- cie_sky_image(z, a,
                                  model$sun_coord$zenith_azimuth,
                                  model$coef) * model$zenith_dn
  names(sky_cie) <- "CIE sky"
  sky_cie
}

.noise <- function(w = 1) {
  path <- system.file("external", package = "rcaiman")
  skies <- utils::read.csv(file.path(path, "15_CIE_standard_skies.csv"))
  .sd <- apply((skies[, 1:5]), 2, sd) * w
  Map(function(i) stats::rnorm(1, 0, .sd[i]), 1:5) %>% unlist()
}

.filter <- function(ds, col_names, thr) {
  d <- as.matrix(stats::dist(ds[, col_names]))
  indices <- c()
  i <- 0
  while (i < nrow(d)) {
    i <- i + 1
    indices <- c(indices, row.names(d)[i]) #include the point itself (p)
    x <- names(d[i, d[i,] <= thr])
    if (!is.null(x)) {
      # this exclude from future search all the points near p,
      # including itself
      rows2crop <- (1:nrow(d))[match(x, rownames(d))]
      cols2crop <- (1:ncol(d))[match(x, colnames(d))]
      d <- d[-rows2crop, -cols2crop]
    }
    if (is.vector(d)) d <- matrix(d)
  }
  ds[indices,]
}
