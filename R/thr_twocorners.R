# path <- system.file("external/iphone_olloclip.jpg", package = "rcaiman")
# caim <- read_caim(path)
# z <- zenith_image(2132/2, lens("Olloclip"))
# a <- azimuth_image(z)
# zenith_colrow <- c(1063, 771)/2
# caim <- expand_noncircular(caim, z, zenith_colrow)
# m <- !is.na(caim$Red) & !is.na(z)
# r <- gbc(caim$Blue)
# r[is.na(r)] <- 0
# # r <- normalize(caim$Blue)
#
# bin <- find_sky_pixels(r, z, a)
# sky <- ootb_sky_reconstruction(r, z, a, bin)
# ratio <- r/sky
# r <- normalize(ratio, 0,1, TRUE)
#
#
# z <- zenith_image(ncol(r), lens())
# m <- mask_hs(z, 50,60)
# r <- normalize(r, min(r[]), max(r[m]), TRUE)
# plot(r)
# m <- !is.na(z)
#
# dns <- r[m]*255 %>% round()
# breaks <- seq(min(dns), max(dns), 1) %>% length()
# dns <- graphics::hist(dns, breaks = breaks, plot = FALSE)
# plot(dns$counts, type = "l")
#
# # dn_max_left <- max()
