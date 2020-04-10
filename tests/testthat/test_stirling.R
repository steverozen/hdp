test_that("simple test of stirling numbers", {
  sfn <- make.stirling()
  expect_equal(sfn(5), c(0.48, 1.00, 0.70, 0.20, 0.02))
})

test_that("stirling numbers; have a zero", {
  sfn <- make.stirling()
  foo <- sfn(200)
  expect_equal(min(which(foo == 0)), 186)
  expect_equal(foo[15], 0.0007595248)
  foo <- sfn(220)
  expect_equal(min(which(foo == 0)), 191)
  expect_equal(foo[15], 0.000964656)
})

test_that("randnumtable", {
  if (exists("stir.closure", envir = .GlobalEnv)) {
  rm("stir.closure", envir = .GlobalEnv)
  }
  set.seed(1066)
  foo <- randnumtable(rep(1, 6) / 6,
                      c(7068, 6864, 7207, 7050, 7039, 0))
  assign("stir.closure", make.stirling(), .GlobalEnv)
  set.seed(1066)
  bar <- randnumtable(rep(1, 6) / 6,
                      c(7068, 6864, 7207, 7050, 7039, 0))
  rm("stir.closure", envir = .GlobalEnv)
  expect_equal(foo, c(2, 2, 1, 3, 3, 0))
  expect_equal(bar, c(2, 2, 1, 3, 3, 0))
})
