.onLoad <- function(libname, pkgname) {
  message(".onLoad ", libname, " ", pkgname)
  print(search())
  print(asNamespace("hdpx"))
  assign("stir.closure", xmake.s(), envir = asNamespace("hdpx"))
  # stir.closure <<- xmake.s()
}
