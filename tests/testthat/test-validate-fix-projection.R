test_that("Mait's rounded coefficients require an explicit boundary correction", {
  original <- c(a=.702,b=.0434,c=-.0543)
  result <- validate_projection_function(original)
  expect_identical(result$status,"underestimation")
  expect_equal(result$val_reference,0.999329126194,tolerance=1e-11)
  expect_true(all(is.na(calc_off_axis_angle(c(.9995,1),original))))
  corrected <- fix_projection_coef(original)
  expect_identical(corrected[1:2],original[1:2])
  expect_identical(names(corrected),names(original))
  expect_gt(corrected[3],original[3])
  expect_identical(validate_projection_function(corrected)$status,"valid")
  expect_true(all(is.finite(calc_off_axis_angle(c(.9995,1),corrected))))
  geometry <- list(model="radial_projection",parameters=original,
    center_of_vision=c(693,0),normalization_radius=693,normalization_angle=90)
  expect_error(set_internal_geometry(c(1386,20),geometry),"fix_projection_coef")
  geometry$parameters <- corrected
  g <- set_internal_geometry(c(1386,20),geometry)
  xy <- coord_from_indices(rep(1:1386,20),rep(1:20,each=1386))
  inside <- (xy$x-693)^2+xy$y_down^2 <= 693^2
  expect_true(all(is.finite(terra::values(g$off_axis,mat=FALSE)[inside])))
})

test_that("normalization angle and inversion limit have separate roles", {
  original <- c(0,0,.9995/(pi/6)^3)
  corrected <- fix_projection_coef(original,normalization_angle=30,max_off_axis_angle=120)
  expect_identical(corrected[1:2],original[1:2])
  expect_gte(calc_relative_radius(30,corrected),1)
  expect_lte(calc_relative_radius(30,corrected),1.001)
  expect_identical(validate_projection_function(corrected,normalization_angle=30,
                   max_off_axis_angle=120)$status,"valid")
  expect_identical(fix_projection_coef(equidistant_coef),equidistant_coef)
  expect_identical(validate_projection_function(3/pi,normalization_angle=60,
                   max_off_axis_angle=120)$status,"valid")
})

test_that("a normalization deficit never hides a folding curve", {
  folded <- .9995*c(2,-1)
  result <- validate_projection_function(folded,normalization_angle=180/pi,
                                         max_off_axis_angle=90)
  expect_identical(result$status,"structural_invalid")
  expect_match(result$reason,"strictly increasing")
  expect_error(fix_projection_coef(folded,normalization_angle=180/pi),"strictly increasing")
  expect_identical(validate_projection_function(1)$status,"structural_invalid")
  expect_error(fix_projection_coef(1),"must equal one")
  expect_error(validate_projection_function(numeric()),"empty")
  expect_error(validate_projection_function(NA_real_),"NA")
  expect_error(validate_projection_function(1,normalization_angle=181),"180")
  expect_error(validate_projection_function(1,max_off_axis_angle=0),"positive")
})
