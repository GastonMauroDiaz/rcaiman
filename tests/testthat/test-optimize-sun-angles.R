sun_search_fixture <- function(sky_no = 12, displaced = FALSE) {
  sun <- data.frame(zenith = 45*pi/180, azimuth = 30*pi/180)
  samples <- expand.grid(zenith = seq(.08,1.45,length.out=6),
                         azimuth = seq(0,6,length.out=9))
  coef <- as.numeric(cie_table[sky_no,1:5])
  samples$relative_radiance <- .cie_sky_model(samples$azimuth,samples$zenith,
    sun$azimuth,sun$zenith,coef[1],coef[2],coef[3],coef[4],coef[5])
  samples$dn <- samples$relative_radiance * 1000
  rr <- list(zenith_dn=1000,sky_points=samples)
  initial <- sun
  if (displaced) initial$zenith <- 50*pi/180
  model <- suppressWarnings(fit_cie_model(rr,initial,std_sky_no=sky_no,method='Nelder-Mead'))
  list(model=model, sun=sun, coef=coef)
}

test_that('caps are bounded, reproducible and work at poles and the azimuth seam', {
  for (center in list(data.frame(zenith=0,azimuth=0),
                      data.frame(zenith=pi/4,azimuth=2*pi-1e-9),
                      data.frame(zenith=95*pi/180,azimuth=1))) {
    cap <- .sun_cap_candidates(center,30,5)
    expect_identical(cap,.sun_cap_candidates(center,30,5))
    expect_equal(cap[1,],center)
    expect_lte(nrow(cap),30)
    expect_true(all(cap$zenith <= 96*pi/180))
    expect_true(all(cap$azimuth >= 0 & cap$azimuth < 2*pi))
    distance <- calc_spherical_distance(cap$zenith,cap$azimuth,center$zenith,center$azimuth)
    expect_true(all(distance <= 5*pi/180 + 1e-12))
    expect_gt(nrow(cap),1)
  }
  center <- data.frame(zenith=pi/4,azimuth=0)
  shifted <- .sun_offset_direction(center,3,4)
  expect_equal(calc_spherical_distance(center$zenith,center$azimuth,
    shifted$zenith,shifted$azimuth)*180/pi,5,tolerance=1e-10)
  expect_equal(nrow(.sun_cap_candidates(center,30,1e-6)),30)
})

test_that('no improvement is valid and does not stop at the first stationary estimate', {
  f <- sun_search_fixture(1)
  expect_true(is.na(f$model$method))
  set.seed(17); rng <- .Random.seed
  result <- optimize_sun_angles(f$model,max_iter=10,max_no_improve=3,
    n_starts=3,method_sun=NULL,std_sky_no=1,method='BFGS')
  expect_identical(.Random.seed,rng)
  expect_equal(result$metric,0,tolerance=1e-20)
  expect_identical(result$sun_angles,f$model$sun_angles)
  expect_null(result$method_sun)
  expect_identical(result$sun_search$stop_reason,'max_no_improve')
  h <- result$sun_search$history
  expect_equal(sum(h$stage=='exploration'),9)
  expect_equal(sum(h$stage=='final'),1)
  expect_identical(result$relative_radiance$sky_points$dn,
                   f$model$relative_radiance$sky_points$dn)
})

test_that('one iteration and skipped stages have finite schedules', {
  f <- sun_search_fixture(1)
  one <- optimize_sun_angles(f$model,max_iter=1,n_starts=2,
    method_sun=NULL,std_sky_no=1,method='BFGS')
  expect_equal(sum(one$sun_search$history$stage=='exploration'),2)
  expect_identical(one$sun_search$stop_reason,'max_iter')
  zero <- optimize_sun_angles(f$model,max_iter=0,method_sun=NULL,
    std_sky_no=1,method='BFGS')
  expect_identical(zero$sun_search$stop_reason,'exploration_skipped')
  expect_equal(nrow(zero$sun_search$history),2)
  expect_equal(nrow(zero$sun_search$refinements),0)
})

test_that('nested fits use fixed starts per search and selectors only at the end', {
  f <- sun_search_fixture(1)
  calls <- list()
  local_mocked_bindings(fit_cie_model=function(relative_radiance,sun_angles,...) {
    opts <- list(...)
    calls[[length(calls)+1L]] <<- c(list(sun=sun_angles,rr=relative_radiance),opts)
    x <- f$model; x$sun_angles <- sun_angles; x
  })
  result <- optimize_sun_angles(f$model,max_iter=1,n_starts=3,
    method_sun=NULL,custom_sky_coef=f$coef+1,optim_zenith_dn=TRUE,
    method_sky='BFGS',method='CG',parallel=TRUE,cores=2)
  expect_length(calls,4)
  for (i in 1:3) {
    expect_identical(calls[[i]]$custom_sky_coef,f$model$coef)
    expect_identical(calls[[i]]$method,'BFGS')
    expect_null(calls[[i]]$optim_zenith_dn)
  }
  expect_identical(calls[[4]]$custom_sky_coef,f$coef+1)
  expect_identical(calls[[4]]$optim_zenith_dn,TRUE)
  expect_identical(calls[[4]]$method,'CG')
  expect_identical(calls[[4]]$parallel,TRUE)
  expect_identical(calls[[4]]$cores,2)
})

test_that('local refinement refits coefficients, caches positions and retains best evaluated model on error', {
  f <- sun_search_fixture(12)
  f$model$metric <- 1
  calls <- list()
  local_mocked_bindings(fit_cie_model=function(relative_radiance,sun_angles,...) {
    opts <- list(...)
    calls[[length(calls)+1L]] <<- c(list(sun=sun_angles),opts)
    x <- f$model; x$sun_angles <- sun_angles
    x$metric <- if (sun_angles$zenith > f$sun$zenith + .001) .1 else 1
    x$coef <- f$model$coef + .01
    x
  })
  local_mocked_bindings(optim=function(par,fn,...) {
    expect_equal(fn(c(0,0)),1/mean(f$model$relative_radiance$sky_points$relative_radiance))
    cost <- fn(c(0,1))
    expect_equal(fn(c(0,1)),cost)
    expect_equal(fn(c(0,2)),Inf)
    stop('deliberate optimizer error')
  },.package='stats')
  result <- optimize_sun_angles(f$model,max_iter=0,max_sun_evaluations=2,
    method_sun='BFGS',method='BFGS')
  expect_length(calls,3) # two unique inner positions and one final comparison
  expect_equal(result$metric,.1)
  expect_equal(result$coef,f$model$coef+.01)
  expect_equal(result$sun_angles$zenith,46*pi/180)
  expect_identical(result$method_sun,'BFGS')
  expect_equal(result$sun_search$refinements$evaluations,2)
  expect_true(result$sun_search$refinements$budget_exhausted)
  expect_match(result$sun_search$refinements$message,'deliberate optimizer error')
  expect_identical(calls[[1]]$custom_sky_coef,calls[[2]]$custom_sky_coef)
  expect_null(calls[[3]]$custom_sky_coef) # all 15 standard skies, by default
})

test_that('the penalty center moves and twilight candidates can be selected', {
  f <- sun_search_fixture(12)
  f$model$sun_angles <- data.frame(zenith=65*pi/180,azimuth=.3)
  f$model$metric <- 1
  anchors <- numeric()
  local_mocked_bindings(fit_cie_model=function(relative_radiance,sun_angles,...) {
    x <- f$model; x$sun_angles <- sun_angles
    x$metric <- exp(-20*(sun_angles$zenith-65*pi/180))
    x
  })
  local_mocked_bindings(optim=function(par,fn,...) {
    anchors <<- c(anchors,fn(c(0,0)))
    fn(c(0,15))
    list(par=c(0,15),convergence=0L)
  },.package='stats')
  result <- optimize_sun_angles(f$model,max_iter=0,
    method_sun=c('BFGS','CG'),method='BFGS')
  expect_equal(result$sun_angles$zenith,95*pi/180)
  h <- result$sun_search$history
  hit <- which(h$stage=='refinement' & h$zenith > 94*pi/180)
  mean_r <- mean(f$model$relative_radiance$sky_points$relative_radiance)
  expect_equal(h$score[hit],sqrt(h$metric[hit])/mean_r * 1.5)
  expect_lt(anchors[2],anchors[1])
})

test_that('failed final comparison returns a usable preceding model and diagnostics', {
  f <- sun_search_fixture(1)
  local_mocked_bindings(fit_cie_model=function(...) stop('deliberate fit error'))
  expect_warning(result <- optimize_sun_angles(f$model,max_iter=0,method_sun=NULL),
    'final sky comparison failed')
  expect_identical(result$coef,f$model$coef)
  expect_match(tail(result$sun_search$history$message,1),'deliberate fit error')
})

test_that('a displaced synthetic sun improves with joint refinement', {
  f <- sun_search_fixture(12,TRUE)
  original <- f$model
  result <- optimize_sun_angles(original,max_iter=2,n_starts=5,
    max_sun_evaluations=20,std_sky_no=12,method='Nelder-Mead')
  expect_lt(result$metric,original$metric)
  expect_equal(result$relative_radiance$zenith_dn,original$relative_radiance$zenith_dn)
  expect_identical(result$relative_radiance$sky_points$dn,original$relative_radiance$sky_points$dn)
  p <- result$relative_radiance$sky_points
  expected <- .cie_cached_predictor(p$zenith,p$azimuth,result$sun_angles)(result$coef)
  expect_equal(p$pred,expected)
  expect_equal(result$metric,mean((expected-p$relative_radiance)^2))
  expect_identical(f$model,original)
})

test_that('invalid controls and old solar formats are rejected before fitting', {
  f <- sun_search_fixture(1)
  expect_error(optimize_sun_angles(f$model,start_radius=0),'positive')
  expect_error(optimize_sun_angles(f$model,tol_angle=6,start_radius=5),'tol_angle')
  expect_error(optimize_sun_angles(f$model,n_starts=0),'positive')
  expect_error(optimize_sun_angles(f$model,max_iter=-1),'non-negative')
  expect_error(optimize_sun_angles(f$model,max_sun_evaluations=0),'positive')
  expect_error(optimize_sun_angles(f$model,method_sky='SANN'),'method_sky')
  expect_error(optimize_sun_angles(f$model,method_sun='Brent'),'method_sun')
  expect_error(optimize_sun_angles(f$model,custom_sky_coef=f$coef,std_sky_no=1),'combine')
  expect_error(optimize_sun_angles(f$model,std_sky_no=16),'1 to 15')
  f$model$sun_angles <- data.frame(zenith=97*pi/180,azimuth=0)
  expect_error(suppressWarnings(optimize_sun_angles(f$model)),'96 degrees')
})


test_that('joint fitting can cross the horizon for a known twilight sky', {
  f <- sun_search_fixture(12)
  sun <- data.frame(zenith=92*pi/180,azimuth=30*pi/180)
  rr <- f$model$relative_radiance
  p <- rr$sky_points
  p$relative_radiance <- .cie_cached_predictor(p$zenith,p$azimuth,sun)(f$coef)
  p$dn <- p$relative_radiance * 1000
  rr$sky_points <- p
  initial <- data.frame(zenith=88*pi/180,azimuth=sun$azimuth)
  model <- suppressWarnings(fit_cie_model(rr,initial,std_sky_no=12,method='Nelder-Mead'))
  result <- optimize_sun_angles(model,max_iter=3,n_starts=6,
    max_sun_evaluations=40,std_sky_no=12,method='Nelder-Mead')
  expect_gt(result$sun_angles$zenith,pi/2)
  expect_lt(result$metric,model$metric/10)
  expect_lt(calc_spherical_distance(sun$zenith,sun$azimuth,
    result$sun_angles$zenith,result$sun_angles$azimuth),1*pi/180)
})


test_that('solar SANN follows the caller seed without changing RNG kind', {
  f <- sun_search_fixture(1)
  kind <- RNGkind()
  run <- function() optimize_sun_angles(f$model,max_iter=0,method_sun='SANN',
    max_sun_evaluations=3,std_sky_no=1,method='BFGS')
  set.seed(317)
  first <- run(); state <- .Random.seed
  set.seed(317)
  second <- run()
  expect_identical(first,second)
  expect_identical(.Random.seed,state)
  expect_identical(RNGkind(),kind)
})
