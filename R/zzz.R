.onLoad <- function(libname, pkgname) {
  # message(".onLoad ", libname, " ", pkgname)
  # print(search())
  # print(asNamespace("hdpx"))
  assign("stir.closure", xmake.s(), envir = asNamespace("hdpx"))
  # stir.closure <<- xmake.s()
}

.onAttach <- function(libname, pkgname) {
  uu <- utils::sessionInfo()$otherPkgs
  extra <- ifelse(is.null(uu), "", uu$hdpx$Version)
  packageStartupMessage("This is ", libname, " and ", pkgname,
                        " v", extra)
  # packageStartupMessage("Run citation('hdp') for citation instructions,
  #  and file.show(system.file('LICENSE', package='hdp')) for license details.")
}

