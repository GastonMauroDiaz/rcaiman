# rcaiman 1.0.6.9000

## Breaking changes
   
* Change dependency from raster to terra
* Fix local_fuzzy_thr()
* fix expand_noncircular()
* Add examples of restricted view photography 

## New features

* New HSP functions family enable a dinamic workflow between R and the HSP software (Tartu Observatory). 
* `azimuth_image()` gains a `rotation` parameter that allows processing images whose top is oriented to any known azimuth angle.
* New `chessboard()` provides chessboard segmentation.
* New `cie_sky_model_raster()` that produces CIE sky images from parameterized CIE sky models. The CIE sky model was implemented based on Pascal code by Mait Lang.
* New `colorfulness()` provides a method to quantify image colorfulness.
* New `deffuzify()` provides a method to turn fuzzy classification into Boolean, which is an alternative to `apply_thr()`.
* New `extract_dn()` facilitates the extraction of digital number from canopy photographs. The extraction is based on points raster coordinates obtained automatically with `extract_sky_points()` or manually by working with third-party software.
* New `extract_rl()` facilitates the extraction of relative luminance from hemispherical photographs. This function uses objects from `extract_sky_points()` and returns objects essential to `fit_cie_sky_model()`.
* New `extract_sky_points()` to automatically extract sky points from canopy photos. 
* New `extract_sun_coord()` to automatically extract sky points from canopy photos. Objects returned by this function are essential to `fit_cie_sky_model()`.
* `find_sky_pixels()` now uses sample size percentage.
* New `find_sky_pixels_nonnull_criteria()` offers a method for fine-tuning working binarized images. It is based on the assumption that the threshold can be tuned as long as no new cells with zero gaps are obtained.
* New `fisheye_to_pano()` provides a method to reproject from hemispherical to cylindrical. The image resolution is substantially degraded since it is based on a sky grid produced with `sky_grid_segmentation()`.
* New `fit_cie_sky_model()` uses maximum likelihood to estimate the coefficients of the CIE sky model that best fit to data sampled from a real scene. Then, the model can be used to produce and image with `cie_sky_model_raster()`.
* New `interpolate_sky_points()` provides a method to produce raster images from point-like data, such as the objects returned by `extract_dn()` or `extract_rl()`.
* `fit_coneshaped_model()` now works with the point-like data objects returned by `extract_rl()` and returns a function that can easily produce a raster if SpatRaster objects are provided as arguments. As a result, the function is faster and more versatile.

  ```R
r <- read_caim()
z <- zenith_image(ncol(r), lens("Nikon_FCE9"))
a <- azimuth_image(z)
g <- sky_grid_segmentation(z, a, 10)
bin <- find_sky_pixels(r, z, a)
sky_points <- extract_sky_points(r, bin, g)
sky_points <- extract_rl(r, z, a, sky_points, NULL)
model <- fit_coneshaped_model(sky_points$sky_points)
model$fun(60, 10)
model$fun(z, a)
  ```

* New `Mask_sunlit_canopy()` is a wrapper function around `membership_to_color()` that facilitates masking sunlit canopy. 
* New `obia()` is a revised version of the object-based image analisys presented in <doi:10.1109/lgrs.2015.2425931>.
* New `ootb_obia()` is a revised version of the full workflow presented in <doi:10.1109/lgrs.2015.2425931>, which includes `enhance_caim()` and `obia()`.
* New `ootb_sky_reconstruction()` provides an easy to use function that will build an above canopy image from a single below canopy image, by means of `fit_cie_sky_model()` and `interpolate_sky_points()`.
* New `polar_qtree()` provides quad-tree segmentation in the polar space.
* New `qtree()` provides classical quad-tree segmentation
* `read_caim()` now is allowed to read any raster image that `terra::raster()` can read. Of course, the georeferencing is turned off by assigning a local projection and manipulating extension and resolution, as usual.

  ```R
f <- system.file("ex/elev.tif", package="terra")
read_caim(f)
terra::rast(f)
  ```

* New `thr_isodata()` alternative implementation of the IsoData method from the autothresholdr package. 


## Minor improvements and fixes
* reproject_to_equidistant() to fisheye_to_equidistant()
* calc_zenith_raster_coordinates() to calc_zenith_raster_coord()
* minor fix for fisheye_to_equidistant() 
* fix_predicted_sky to fix_reconstructed_sky
