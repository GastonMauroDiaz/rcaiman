cie_fit_fixture <- function(sky_no = 4, reference = 1000) {
  sun <- data.frame(zenith = pi/4, azimuth = .3)
  directions <- expand.grid(zenith = seq(.1, 1.4, length.out = 6),
    azimuth = seq(0, 5.5, length.out = 8))
  coef <- as.numeric(cie_table[sky_no, 1:5])
  directions$relative_radiance <- .cie_sky_model(directions$azimuth, directions$zenith,
    sun$azimuth, sun$zenith, coef[1], coef[2], coef[3], coef[4], coef[5])
  directions$dn <- directions$relative_radiance * reference
  list(relative_radiance = list(zenith_dn = reference, sky_points = directions), sun = sun, coef = coef)
}

test_that("every initial sky is tried with every optimization method", {
  calls <- list()
  local_mocked_bindings(optim = function(par, fn, method, ...) {
    calls[[length(calls)+1L]] <<- list(par = par, method = method)
    list(par = par, value = fn(par), convergence = 0L)
  }, .package = "stats")
  f <- cie_fit_fixture()
  methods <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")
  result <- fit_cie_model(f$relative_radiance, f$sun, std_sky_no = c(3,4,5), method = methods)
  expect_length(calls, 12)
  expect_identical(vapply(calls, `[[`, character(1), "method"), rep(methods,3))
  expect_equal(nrow(result$attempts), 15)
  expect_identical(result$source, "standard")
  expect_equal(result$coef, setNames(f$coef, letters[1:5]))
  expect_null(result$opt_result)
  expect_true(is.na(result$method))
  expect_equal(result$metric, 0)
  expect_equal(sum(result$attempts$selected),1)
  expect_equal(result$relative_radiance$sky_points$dn, f$relative_radiance$sky_points$dn)
})

test_that("failed optimizations retain the best initial model and diagnostics", {
  local_mocked_bindings(optim = function(...) stop("deliberate failure"), .package = "stats")
  f <- cie_fit_fixture()
  expect_warning(result <- fit_cie_model(f$relative_radiance,f$sun,std_sky_no=c(3,4),method="BFGS"),
    "returning the best initial model")
  expect_identical(result$source,"standard")
  expect_equal(result$metric,0)
  expect_true(all(grepl("deliberate failure",
    result$attempts$message[result$attempts$source == "optimized"])))
  expect_warning(custom <- fit_cie_model(f$relative_radiance,f$sun,custom_sky_coef=f$coef,method="BFGS"),
    "returning the best initial model")
  expect_identical(custom$source,"custom")
})

test_that("zenith-reference optimization preserves DN and reports adjusted observations", {
  f <- cie_fit_fixture()
  f$relative_radiance$zenith_dn <- 1200
  f$relative_radiance$sky_points$relative_radiance <- f$relative_radiance$sky_points$dn / 1200
  original <- f$relative_radiance
  local_mocked_bindings(optim = function(par, fn, method, ...) {
    par[1:5] <- f$coef
    par[6] <- (1000 - mean(original$sky_points$dn)) / sd(original$sky_points$dn)
    list(par = par, value = fn(par), convergence = 0L)
  }, .package = "stats")
  result <- fit_cie_model(f$relative_radiance, f$sun, std_sky_no=3, method="BFGS",optim_zenith_dn=TRUE)
  expect_identical(result$source,"optimized")
  expect_equal(result$relative_radiance$zenith_dn,1000)
  expect_equal(result$relative_radiance$sky_points$relative_radiance,original$sky_points$dn/1000)
  expect_identical(result$relative_radiance$sky_points$dn,original$sky_points$dn)
  expect_identical(f$relative_radiance,original)
  expect_equal(result$metric,0,tolerance=1e-20)
  expect_length(result$start,6)
})

test_that("invalid trial references are penalized and constant DN keep the reference fixed", {
  f <- cie_fit_fixture()
  local_mocked_bindings(optim = function(par, fn, method, ...) {
    if(length(par)==6) {
      bad <- par; bad[6] <- -1e10
      expect_identical(fn(bad),Inf)
    }
    list(par = par,value=fn(par),convergence=0L)
  }, .package="stats")
  fit_cie_model(f$relative_radiance,f$sun,std_sky_no=4,method="BFGS",optim_zenith_dn=TRUE)
  f$relative_radiance$sky_points$relative_radiance <- rep(1,nrow(f$relative_radiance$sky_points))
  expect_warning(result <- fit_cie_model(f$relative_radiance,f$sun,std_sky_no=5,method="BFGS",optim_zenith_dn=TRUE),
    "keeping `zenith_dn` fixed")
  expect_equal(result$relative_radiance$zenith_dn,1000)
  expect_length(result$start,5)
})

test_that("real optimization improves an initial model and reconstructs a CIE image", {
  f <- cie_fit_fixture()
  result <- fit_cie_model(f$relative_radiance,f$sun,custom_sky_coef=f$coef+c(.2,0,0,0,0),method="BFGS")
  expect_identical(result$source,"optimized")
  expect_lt(result$metric, result$attempts$metric[1])
  expect_equal(sum(result$components),result$metric)
  d <- set_rcaiman_geometry(terra::rast(nrows=6,ncols=8,nlyrs=2))
  terra::values(d) <- as.matrix(f$relative_radiance$sky_points[c("zenith","azimuth")])*180/pi
  names(d) <- c("zenith","azimuth")
  predicted <- cie_image(d,result$sun_angles,result$coef)
  expect_equal(as.vector(predicted[]),result$relative_radiance$sky_points$pred)
})

test_that("CIE fitting validates the new table contracts and selectors", {
  f <- cie_fit_fixture()
  expect_error(fit_cie_model(f$relative_radiance,c(zen_deg=45,azi_deg=10)),"data frame")
  bad <- f$relative_radiance; bad$sky_points$zenith[1] <- 2
  expect_error(fit_cie_model(bad,f$sun),"radians")
  bad <- f$relative_radiance; bad$sky_points$relative_radiance[1] <- NA_real_
  expect_error(fit_cie_model(bad,f$sun),"NA")
  expect_error(fit_cie_model(f$relative_radiance,f$sun,std_sky_no=numeric()),"1 to 15")
  expect_error(fit_cie_model(f$relative_radiance,f$sun,method="Brent"),"multivariate")
  expect_error(fit_cie_model(f$relative_radiance,f$sun,custom_sky_coef=f$coef,std_sky_no=4),"combine")
  expect_error(fit_cie_model(f$relative_radiance,f$sun,custom_sky_coef=rep(NA_real_,5)),"finite")
})

test_that("real reference optimization improves a misnormalized sample", {
  f <- cie_fit_fixture()
  f$relative_radiance$zenith_dn <- 1200
  f$relative_radiance$sky_points$relative_radiance <- f$relative_radiance$sky_points$dn / 1200
  result <- fit_cie_model(f$relative_radiance, f$sun, custom_sky_coef=f$coef,
    method="BFGS", optim_zenith_dn=TRUE)
  expect_identical(result$source,"optimized")
  expect_lt(result$metric,result$attempts$metric[1])
  expect_gt(result$relative_radiance$zenith_dn,0)
  expect_equal(result$relative_radiance$sky_points$relative_radiance * result$relative_radiance$zenith_dn, f$relative_radiance$sky_points$dn)
})

test_that("usable nonconverged results and warning messages remain visible", {
  f <- cie_fit_fixture()
  local_mocked_bindings(optim = function(par, fn, method, ...) {
    warning("test optimizer warning")
    par[] <- f$coef
    list(par=par,value=fn(par),convergence=1L,message="iteration limit")
  }, .package="stats")
  expect_warning(result <- fit_cie_model(f$relative_radiance,f$sun,std_sky_no=3,method="BFGS"),
    "selected optimization did not converge")
  expect_identical(result$source,"optimized")
  selected <- result$attempts[result$attempts$selected,]
  expect_equal(selected$convergence,1)
  expect_match(selected$message,"test optimizer warning")
  expect_match(selected$message,"iteration limit")
})


test_that("cached geometry preserves model predictions for all standard skies", {
  sun <- data.frame(zenith=pi/4,azimuth=.3)
  z <- c(0,.1,.5,1,pi/2); a <- c(0,3,2,1,.3)
  cached <- .cie_cached_predictor(z,a,sun)
  for(i in 1:15) {
    p <- as.numeric(cie_table[i,1:5])
    expected <- .cie_sky_model(a,z,sun$azimuth,sun$zenith,p[1],p[2],p[3],p[4],p[5])
    expect_identical(cached(p),expected)
  }
})

test_that("parallel deterministic fitting matches sequential without using RNG", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 2)
  f <- cie_fit_fixture()
  set.seed(72)
  before <- .Random.seed
  sequential <- fit_cie_model(f$relative_radiance,f$sun,std_sky_no=c(3,4),method=c("Nelder-Mead","BFGS"))
  expect_identical(.Random.seed,before)
  concurrent <- fit_cie_model(f$relative_radiance,f$sun,std_sky_no=c(3,4),method=c("Nelder-Mead","BFGS"),parallel=TRUE,cores=2)
  expect_identical(.Random.seed,before)
  expect_identical(concurrent,sequential)
  expect_error(fit_cie_model(f$relative_radiance,f$sun,cores=0),"positive")
  expect_error(fit_cie_model(f$relative_radiance,f$sun,parallel=NA),"NA")
})

test_that("SANN streams depend on tasks, not scheduling or worker count", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 2)
  f <- cie_fit_fixture()
  run <- function(parallel=FALSE,cores=NULL) {
    fit_cie_model(f$relative_radiance,f$sun,std_sky_no=c(3,4),method=c("BFGS","SANN"),
      parallel=parallel,cores=cores)
  }
  set.seed(178)
  kinds <- RNGkind()
  sequential <- run()
  next_state <- .Random.seed
  expect_identical(RNGkind(),kinds)
  set.seed(178)
  concurrent <- run(TRUE,2)
  expect_identical(concurrent,sequential)
  expect_identical(.Random.seed,next_state)
  expect_identical(RNGkind(),kinds)
  set.seed(178)
  single <- run(TRUE,1)
  expect_identical(single,sequential)
  expect_identical(.Random.seed,next_state)
})

test_that("SANN preserves the caller's cached normal state", {
  skip_on_cran()
  f <- cie_fit_fixture()
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind,as.list(old_kind)),add=TRUE)
  set.seed(93,kind="Mersenne-Twister",normal.kind="Box-Muller")
  invisible(rnorm(1))
  invisible(sample.int(.Machine$integer.max,1L))
  expected <- c(rnorm(3),runif(3))
  set.seed(93,kind="Mersenne-Twister",normal.kind="Box-Muller")
  invisible(rnorm(1))
  fit_cie_model(f$relative_radiance,f$sun,std_sky_no=4,method="SANN")
  expect_identical(RNGkind()[1:2],c("Mersenne-Twister","Box-Muller"))
  expect_identical(c(rnorm(3),runif(3)),expected)
})


test_that("Fibonacci checking penalizes a negative horizon beyond positive samples", {
  f <- cie_fit_fixture()
  p <- f$relative_radiance$sky_points
  # These coefficients predict positive radiance at every observed direction,
  # but negative radiance beyond about 81.7 degrees of zenith.
  bad <- c(-2, -.1, 0, -1, 0)
  observed_prediction <- .cie_cached_predictor(p$zenith, p$azimuth, f$sun)(bad)
  expect_true(all(observed_prediction > 0))
  vault <- external_directions_from_internal(fibonacci_directions(5))
  expect_equal(nrow(vault), 826)
  expect_true(all(vault$zenith > 0 & vault$zenith < pi/2))
  checked <- .cie_cached_predictor(vault$zenith, vault$azimuth, f$sun)(bad)
  expect_true(any(checked[vault$zenith > 89*pi/180] < 0))
  expected_penalty <- 1e10 * sum(pmin(checked, 0))^2
  base_cost <- sqrt(mean((observed_prediction - p$relative_radiance)^2)) /
    mean(p$relative_radiance)
  seen <- NULL
  local_mocked_bindings(optim = function(par, fn, ...) {
    seen <<- fn(bad)
    expect_identical(fn(bad), seen)
    list(par = par, value = fn(par), convergence = 0L)
  }, .package = "stats")
  fit_cie_model(f$relative_radiance, f$sun, custom_sky_coef = f$coef, method = "BFGS")
  expect_equal(seen, base_cost + expected_penalty)
  expect_gt(seen, base_cost + 1e10)
})

test_that("twilight brightness penalty uses the solar-azimuth horizon, not sub-horizon radiance", {
  sun <- data.frame(zenith=92*pi/180,azimuth=.3)
  f <- cie_fit_fixture(12)
  p <- f$relative_radiance$sky_points
  p$relative_radiance <- .cie_cached_predictor(p$zenith,p$azimuth,sun)(f$coef)
  f$relative_radiance$sky_points <- p
  seen <- NULL
  local_mocked_bindings(optim=function(par,fn,...) {
    seen <<- fn(par)
    list(par=par,value=seen,convergence=0L)
  },.package='stats')
  result <- fit_cie_model(f$relative_radiance,sun,custom_sky_coef=f$coef,method='BFGS')
  expect_lt(.cie_cached_predictor(sun$zenith,sun$azimuth,sun)(f$coef),0)
  expect_equal(seen,0,tolerance=1e-10)
  expect_equal(result$metric,0,tolerance=1e-20)
  expect_identical(result$sun_angles,sun)
})
