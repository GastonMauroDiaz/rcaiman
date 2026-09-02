test_that("rem_nearby_indices() filters in cell-index space", {
  sampling_indices <- data.frame(
    col = c(1, 2, 10),
    row = c(1, 1, 10),
    id = letters[1:3]
  )

  observed <- rem_nearby_indices(sampling_indices, min_dist = 2)

  expect_equal(observed, sampling_indices[c(1, 3), ])
})


test_that("rem_nearby_indices() filters in angular space", {
  sampling_indices <- data.frame(
    col = 1:3,
    row = 1:3,
    zen = rep(pi / 4, 3),
    azi = c(0, pi / 180, pi)
  )

  observed <- rem_nearby_indices(
    sampling_indices,
    min_dist = 10,
    space = "spherical"
  )

  expect_equal(observed, sampling_indices[c(1, 3), ])
})
