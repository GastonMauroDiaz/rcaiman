test_that("Gauch components partition MSD for ordinary and constant data", {
  for (pair in list(list(c(2,4,7,8), c(1,3,3,5)),
                   list(rep(2,4), 1:4), list(1:4, rep(2,4)),
                   list(rep(2,4), rep(3,4)), list(2,3), list(1:4, 1:4))) {
    predicted <- pair[[1]]; observed <- pair[[2]]
    result <- calc_msd(predicted, observed)
    expect_named(result, c("MSD", "SB", "NU", "LC"))
    expect_equal(unname(result["MSD"]), mean((predicted-observed)^2))
    expect_equal(unname(result["MSD"]), sum(result[c("SB", "NU", "LC")]))
    expect_true(all(result >= 0))
  }
  expect_equal(unname(calc_msd(c(2,4,7,8), c(1,3,3,5))),
    c(6.75, 5.0625, 1.26991758241758, .417582417582418))
  expect_equal(calc_msd(rep(2,4), 1:4), c(MSD=1.5, SB=.25, NU=0, LC=1.25))
})

test_that("components agree with observed-on-predicted regression", {
  predicted <- c(1,3,4,7,9)
  observed <- c(2,2,5,6,8)
  reg <- lm(observed ~ predicted)
  result <- calc_msd(predicted, observed)
  expect_equal(unname(result["LC"]), mean(residuals(reg)^2))
  expect_equal(unname(result["NU"]),
    unname((1-coef(reg)[2])^2 * mean((predicted-mean(predicted))^2)))
  reverse <- calc_msd(observed, predicted)
  expect_equal(result[c("MSD","SB")], reverse[c("MSD","SB")])
  expect_false(isTRUE(all.equal(result[c("NU","LC")], reverse[c("NU","LC")])))
  expect_equal(calc_msd(predicted+100, observed+100), result)
  expect_equal(calc_msd(predicted*3, observed*3), result*9)
})

test_that("MSD requires finite nonempty paired vectors", {
  for (bad in list(NA_real_, NaN, Inf, -Inf, "1", TRUE)) {
    expect_error(calc_msd(bad, 1))
    expect_error(calc_msd(1, bad))
  }
  expect_error(calc_msd(numeric(), numeric()))
  expect_error(calc_msd(1:3, 1))
  expect_error(calc_msd(matrix(1:4,2), 1:4))
  expect_equal(calc_msd(c(a=1,b=3), c(b=2,a=4)), c(MSD=1,SB=1,NU=0,LC=0))
})
