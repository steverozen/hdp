
stir.closure <- make.stirling()

hdp.out <- SynSigRun::Runhdp(
  input.catalog =
    devtools::package_file(
      "data-raw/syn.2.7a.7b.abst.369/sp.sp/ground.truth.syn.catalog.csv"),
  out.dir      = "data-raw/syn.2.7a.7b.abst.369/sp.sp/out",
  seedNumber   = 10,
  K.guess      = 20,
  K.range      = NULL,
  multi.types  = FALSE,
  remove.noise = FALSE,
  test.only    = FALSE,
  overwrite    = TRUE,
  verbose      = TRUE
)
