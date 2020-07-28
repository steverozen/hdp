
test_that("hdp_merge_and_extract_components", {

  reg <- new.env()
  load("test.MultipleSetupAndPosterior.Rdata", envir = reg)##This input is taken from mSigHdp

  reg2 <- new.env()
  load("hdp.merge.and.extract.components.expected.Rdata", envir = reg2)



  x <- hdpx::hdp_multi_chain(reg$retvalx)
  retvalx <- hdp_merge_and_extract_components(x = x,
                                              categ.CI = 0.95,
                                              exposure.CI = 0.95,
                                              cos.merge = 0.90,
                                              min.sample = 1)



  #save(retvalx, file = "hdp.merge.and.extract.components.expected.Rdata")

  expect_equal(retvalx, reg2$retvalx)


})
