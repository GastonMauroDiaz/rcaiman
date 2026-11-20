require(rcaiman.registry.data)
summary(entry_Nikon_Coolpix5700.FCE9.CIEFAP)
export_registry_snapshot(entry_Nikon_Coolpix5700.FCE9.CIEFAP,
                         "jpeg", "rgb", "simple_method_rescaled", "wagner2001_method",
                         "inst/external/Nikon_Coolpix5700.FCE9_jpeg.yaml")
export_registry_snapshot(entry_Nikon_Coolpix5700.FCE9.CIEFAP,
                         "raw", "cygm", "simple_method", "photometric_sphere",
                         "inst/external/Nikon_Coolpix5700.FCE9_raw.yaml")
