
test_that("hdp_merge_and_extract_components", {

  reg <- new.env()
  load("test.MultipleSetupAndPosterior.Rdata", envir = reg)##This input is taken from mSigHdp

  reg2 <- new.env()
  load("extract_components_from_clusters.expected.Rdata", envir = reg2)



  x <- hdpx::hdp_multi_chain(reg$retvalx)
  retvalx <- extract_components_from_clusters(x = x,
                                              cos.merge = 0.90,
                                              hc.cutoff = 0.10
                                              )



  #save(retvalx, file = "extract_components_from_clusters.expected.Rdata")

  expect_equal(retvalx, reg2$retvalx)


})
