require(rcaiman.registry)
path <- system.file("external/iPhone6plus.Olloclip_jpeg.yaml", package = "rcaiman")
my_entry <- read_registry_snapshot(path)
spec <- get_geometry_spec(my_entry,  "exif_jpeg", "simple_method_rescaled")

path <- "C:/Users/gdiaz/Downloads/WhatsApp Image 2026-06-10 at 3.23.32 PM.jpeg"
caim <- read_caim(path)
caim

1125
zenith_col_row <- c(642, 340)
diameter <- calc_diameter(spec$parameters, 1125/2, spec$max_zenith_angle) - 1

z <- zenith_image(diameter,  spec$parameters)
a <- azimuth_image(z)
caim <- expand_non_circular(caim, z, zenith_col_row)
plot(caim$Blue, col = seq(0, 1, 1/255) %>% grey())
m <- !is.na(caim$Red) & !is.na(z)
plot(m, add = TRUE, alpha = 0.3, legend = FALSE)

write_caim(caim, "ref", 8)

pts_ref <- read.csv("C:/Users/gdiaz/Downloads/Results_ref.csv")
pts_ref <- pts_ref[c("Y", "X")]
colnames(pts_ref) <- c("row", "col")

pts <- read.csv("C:/Users/gdiaz/Downloads/Results.csv")
pts <- pts[c("Y", "X")]
colnames(pts) <- c("row", "col")

pts_ref <- add_directions(pts_ref, z, a)

plot(caim$Blue)
points(pts_ref$col, nrow(caim$Blue) - pts_ref$row, col = 2, pch = 10)

r <- read_caim("C:/Users/gdiaz/Downloads/WhatsApp Image 2026-06-10 at 3.23.32 PM (1).jpeg")
plot(r$Blue)

# center in (0,0)
pts$y <- pts$row - nrow(r)/2
pts$x <- pts$col - ncol(r)/2
xy <- as.matrix(pts[, c("x", "y")])
pol <- geometry::cart2pol(xy)
max_r_px <- sqrt((ncol(r)/2)^2 + (nrow(r)/2)^2)

pts$zen <- pts_ref$zen
pts$azi <- pts_ref$azi
pts$r <- (pol[, 2]) / max_r_px

# Zenith
pts$atan_r <- atan(pts$r / 0.6)
modelo_zenital <- lm(zen ~ atan_r + I(atan_r^2) + I(atan_r^3), data = pts)
summary(modelo_zenital)

plot(pts$zen, predict(modelo_zenital), xlim = c(0, pi/2), ylim = c(0, pi/2))
abline(0,1)

#Azimuth
pts$azi_img <- pol[,1]

# 1. Descompones la variable respuesta (real)
pts$sin_real <- sin(pts$azi)
pts$cos_real <- cos(pts$azi)

# 2. Ajustas los modelos usando el azimut geométrico de la imagen como predictor
modelo_sin <- lm(sin_real ~ sin(azi_img) + cos(azi_img) - 1, data = pts)
modelo_cos <- lm(cos_real ~ sin(azi_img) + cos(azi_img) - 1, data = pts)

# Modelo complejo para azimut (solo si hay distorsión tangencial rara)
modelo_sin_complejo <- lm(sin_real ~ (sin(azi_img) + cos(azi_img)) * r, data = pts)

summary(modelo_cos)
summary(modelo_sin)
summary(modelo_sin_complejo)

#Calcular para toda la imagen
row_col <- rowColFromCell(r, 1:ncell(r))
row_col[, 1] <- row_col[, 1] - nrow(r)/2
row_col[, 2] <- row_col[, 2] - ncol(r)/2
pol_img <- geometry::cart2pol(row_col)
r_norm <- (pol_img[, 2]) / max_r_px
df <- data.frame(r = r_norm)
df$atan_r <- atan(df$r / 0.6)

zpan <- predict(modelo_zenital, df)
summary(zpan)

df$azi_img <- pol_img[, 1]

apan <- atan2(predict(modelo_cos, df),
              predict(modelo_sin, df))
# Pasar de [-pi, pi] a [0, 2pi]
apan <- apan %% (2 * pi)
summary(apan)

zpan_img <- apan_img <- r$Blue

zpan_img[] <- zpan
zpan_img[zpan < 0] <- 0
plot(zpan_img == 0)

apan_img[] <- apan
plot(apan_img)

#reprojectar
row_col <- rowColFromCell(r, 1:ncell(r))
row_col <- data.frame(row_col)
names(row_col) <- c("row", "col")
row_col$azi <- apan
row_col$zen <- zpan
row_col$dn <- r$Blue[]
names(row_col)

pts_reproj <- reproject_points(row_col, equal_area_coef, 1000)
head(pts_reproj)
z <- zenith_image(1000, equal_area_coef)
a <- azimuth_image(z)
r2 <- interpolate_planar(pts_reproj, "dn", z, search_radius = 5)
plot(r2)
