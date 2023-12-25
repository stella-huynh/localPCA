.onLoad <- function(libname = find.package("lostruct"), pkgname="lostruct") {
  repos = getOption("repos")
  repos["Github"] = "https://github.com/petrelharp/local_pca/lostruct"
  options(repos = repos)
  invisible(repos)
}

.onAttach <-  function(libname = find.package("localPCA"), pkgname = "localPCA"){
  packageStartupMessage("This R package is preliminary or provisional and is subject to revision.")
}
