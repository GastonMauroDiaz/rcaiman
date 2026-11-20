# Degrade sensor data to 8-bits to save space in packaging

require(rcaiman)
rm(list = ls())

sensor_data <- read_caim("sensor_data/DSCN4248.tif")
sensor_data <- normalize_min_max(sensor_data) * 255
write_caim(sensor_data, path = "img_l0/DSCN4248.tif", 8)

