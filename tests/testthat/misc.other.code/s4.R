
# input.catalog <- ICAMS::ReadCatalog("cat.csv")
load("input.catalog.Rdata")

colnames(input.catalog)[1:7] <- paste0("a::", colnames(input.catalog)[1:7])
colnames(input.catalog)[8:15] <- paste0("b::", colnames(input.catalog)[8:15])

input.catalog[3, ] <- 0

input.catalog[5, 1:7] <- 0
input.catalog[6, 8:15] <- 0

retvalx <- hdpx::TestScaffold1(
  input.catalog = input.catalog[ , 1:15],
  CPU.cores     = 1,
  seedNumber    = 1000,
  K.guess       = 5,
  multi.types   = TRUE,
  verbose       = TRUE,
  post.burnin   = 50,
  post.n        = 50,
  post.space    = 25,
  post.cpiter   = 3,
  num.posterior = 1
)
# nice valgrind --leak-check=yes ~/debug.hdp/R-4.0.0/bin/exec/R --vanilla < s4.R >s4.out &>s4.err &
