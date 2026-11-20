defuzzify_test_layers <- function(values, labels) {
  mem <- set_rcaiman_geometry(terra::rast(nrows=1,ncols=length(values)))
  terra::values(mem) <- values
  seg <- terra::rast(mem)
  terra::values(seg) <- labels
  list(mem=mem,seg=seg)
}

test_that("defuzzification respects excluded labels and missing memberships", {
  x <- defuzzify_test_layers(c(.9,.6,.2,1,.8,.3,NA,1),c(1,1,1,0,2,2,2,NA))
  bin <- defuzzify(x$mem,x$seg)
  values <- terra::values(bin,mat=FALSE)
  expect_identical(values,c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE))
  expect_true(terra::compareGeom(bin,x$mem))
  included <- terra::values(defuzzify(x$mem,x$seg,ignore_label_0=FALSE),mat=FALSE)
  expect_true(included[4])
  expect_false(anyNA(included))
  x$seg[] <- NA
  expect_identical(terra::values(defuzzify(x$mem,x$seg),mat=FALSE),rep(FALSE,8))
})

test_that("rounded counts and ties determine exactly which pixels are selected", {
  x <- defuzzify_test_layers(c(.5,1,.5,.5,.5,.5,0,0),c(1,1,2,2,2,2,3,3))
  # round(1.5) = 2, including the .5 pixel; equal memberships use cell order.
  expect_identical(terra::values(defuzzify(x$mem,x$seg),mat=FALSE),
                   c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE))
  x <- defuzzify_test_layers(c(NA,NA,.5),c(1,1,2))
  expect_identical(terra::values(defuzzify(x$mem,x$seg),mat=FALSE),rep(FALSE,3))
})

test_that("invalid memberships warn without breaking the logical output contract", {
  x <- defuzzify_test_layers(c(2,2,-1,-1,NA),c(1,1,2,2,2))
  expect_warning(bin <- defuzzify(x$mem,x$seg),"membership")
  expect_identical(terra::values(bin,mat=FALSE),c(TRUE,TRUE,FALSE,FALSE,FALSE))
  x$seg[] <- 0
  expect_silent(bin <- defuzzify(x$mem,x$seg))
  expect_identical(terra::values(bin,mat=FALSE),rep(FALSE,5))
})
