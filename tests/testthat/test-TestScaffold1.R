
test_that("TestScaffold1", {

  load("input.catalog.Rdata")
  # Get variable input.catalog

  reg <- new.env()
  load("TestScaffold1.Rdata", envir = reg)

  retvalx <- TestScaffold1(
    input.catalog = input.catalog[1:10 , 1:15],
    CPU.cores     = 1,
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = FALSE,
    verbose       = TRUE,
    post.burnin   = 50,
    num.posterior = 1
  )

  # save(retvalx, file = "TestScaffold1.Rdata")

  expect_equal(retvalx, reg$retvalx)
})
