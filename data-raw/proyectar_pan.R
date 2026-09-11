require(rcaiman.registry)
rm(list = ls())
my_entry <- readRDS("data-raw/entry_iPhone6plus.Olloclip_v3.rds")
spec <- get_geometry_spec(my_entry,  "exif_raw", "simple_method")

zenith_col_row <- c(642, 340)
diameter <- calc_diameter(spec$parameters, 1125/2, spec$max_zenith_angle)

z <- zenith_image(diameter,  spec$parameters)
a <- azimuth_image(z)

# read images -------------------------------------------------------------

# path <- "data-raw/WhatsApp Image 2026-06-10 at 3.23.32 PM.jpeg"
# caim <- read_caim(path)
#
# caim <- expand_non_circular(caim, z, zenith_col_row)
#
# r <- read_caim("data-raw/WhatsApp Image 2026-06-10 at 3.23.32 PM (1).jpeg")
# plot(r$Blue)
#
# write_caim(r, "data-raw/target", 8)

ref <- read_caim("data-raw/crosscalibration_wide_angle/img_ref/ref.tif")

r <- read_caim("data-raw/crosscalibration_wide_angle/img_target/target.tif")

# read points -------------------------------------------------------------

# pts_ref <- read.csv("data-raw/Results_ref.csv")
# pts_ref <- pts_ref[c("X", "Y")]
# colnames(pts_ref) <- c("x", "y_down")
#
# pts <- read.csv("data-raw/Results.csv")
# pts <- pts[c("X", "Y")]
# colnames(pts) <- c("x", "y_down")

pts_ref <- read_points("data-raw/crosscalibration_wide_angle/pts_ref/ref.csv")
pts <- read_points("data-raw/crosscalibration_wide_angle/pts_target/target.csv")

# Enrich ------------------------------------------------------------------

pts_ref <- add_indices(pts_ref)
pts_ref <- add_directions(pts_ref, z, a)

# modelling ---------------------------------------------------------------

calibration <- crosscalibrate_geometry(
  pts_uncal = pts,
  pts_cal = pts_ref,
  model = "anisotropic_projection",
  degree = 3,
  width_uncal = ncol(r),
  height_uncal = nrow(r),
  atan_scale = 0.7
)

summary(calibration$model$zen)
summary(calibration$model$azi_cos)
summary(calibration$model$azi_sin)

# Full image --------------------------------------------------------------

all_pixels <- data.frame(col = colFromCell(r, 1:ncell(r)), row = rowFromCell(r, 1:ncell(r)))
all_pixels <- add_coord(all_pixels)
all_pixels <- add_image_polar(all_pixels, ncol(r), nrow(r))
all_pixels <- add_directions_from_model(
  all_pixels,
  calibration$parameters
)


diameter <- calc_diameter(spec$parameters, 1125/2, spec$max_zenith_angle)

all_pixels$dn <- r$Blue[]
pts_reproj <- reproject_points(all_pixels, spec$parameters, diameter)
head(pts_reproj)
r2 <- interpolate_planar(pts_reproj, 6, z, search_radius = 5)
plot(ref$Blue)
plot(r2, add = TRUE, alpha = 0.99)

foo <- reproject_image(ref, z, a, !is.na(z), equidistant_coef, diameter)
plot(foo$Blue)
pts_reproj <- reproject_points(all_pixels, equidistant_coef, diameter)
r2 <- interpolate_planar(pts_reproj, 6, z, search_radius = 5)
plot(r2, add = TRUE, alpha = 0.99)

zfoo <- zenith_image(diameter, equidistant_coef)
m <- select_sky_region(z, 0, 100/2)
plot(m, add = TRUE, alpha = 0.2)
