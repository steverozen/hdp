.onAttach <- function(libname, pkgname) {
  uu <- sessionInfo()$otherPkgs
  extra <- ifelse(is.null(uu), "", uu)
  packageStartupMessage("This is ", libname, " and ", pkgname,
                        " ", extra)
  # packageStartupMessage("Run citation('hdp') for citation instructions,
  #  and file.show(system.file('LICENSE', package='hdp')) for license details.")
}
