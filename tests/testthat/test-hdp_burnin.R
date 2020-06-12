
test_that("hdp_burnin", {

  reg <- new.env()
  load("input.activated.hdp.state.Rdata", envir = reg)

  reg2 <- new.env()
  load("hdp.burnin.expected.Rdata", envir = reg2)

  set.seed(44)

  retvalx <- hdp_burnin(hdp = reg$hdp.state, burnin = 10, cpiter = 3)

  # save(retvalx, file = "hdp.burnin.expected.Rdata")

  expect_equal(reg2$retvalx, retvalx)
})
