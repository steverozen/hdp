test_that("simple test of stirling/stirling2",
          {expect_equal(stirling2(5), c(0.48, 1.00, 0.70, 0.20, 0.02))}
)

test_that("stirling/stirling2, has a zero", {
  foo <- stirling2(200)
  expect_equal(min(which(foo == 0)), 186)
  expect_equal(foo[15], 0.0007595248)
  foo <- stirling2(220)
  expect_equal(min(which(foo == 0)), 191)
  expect_equal(foo[15], 0.000964656)
})


