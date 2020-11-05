
test_that("hdp_posterior_sample", {

  reg <- new.env()
  load("hdp.burnin.expected.Rdata", envir = reg)

  reg2 <- new.env()
  load("hdp.posterior.sample.expected.Rdata", envir = reg2)

  reg3 <- new.env()
  load("extend.hdp.posterior.sample.expected.Rdata", envir = reg3)

  set.seed(44)


  retvalx <- hdp_posterior_sample(post.input      = reg$retvalx,
                                  post.n          = 20,
                                  post.space      = 10,
                                  post.cpiter     = 3,
                                  seed            = 44,
                                  post.verbosity  = 0,
                                  checkpoint = T)

  retvalx.extend <- hdp_posterior_sample(post.input      = retvalx,
                                         post.n          = 20,
                                         post.space      = 10,
                                         post.cpiter     = 3,
                                         seed            = 44,
                                         post.verbosity  = 0,
                                         checkpoint = T)

   save(retvalx, file = "hdp.posterior.sample.expected.Rdata")
   save(retvalx.extend, file = "extend.hdp.posterior.sample.expected.Rdata")

  expect_equal(retvalx, reg2$retvalx)
  expect_equal(retvalx.extend, reg3$retvalx.extend)

})
