
# input.catalog <- ICAMS::ReadCatalog("cat.csv")
load("input.catalog.Rdata")

retvalx <- hdpx::TestScaffold1(
  input.catalog = input.catalog[ , 1:15],
  CPU.cores     = 1,
  seedNumber    = 34,
  K.guess       = 5,
  multi.types   = FALSE,
  verbose       = TRUE,
  post.burnin   = 50,
  post.n        = 50,
  post.space    = 25,
  post.cpiter   = 3,
  num.posterior = 1
)
# nice valgrind --leak-check=yes ~gmssgr/debug.hdp/R-4.0.0/bin/exec/R --vanilla < s1.R >s1.out &>s1.err &
