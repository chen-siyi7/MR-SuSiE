mr_susie_load_native <- function(file, env, cacheDir, showOutput=FALSE, verbose=FALSE) {
  old <- Sys.getenv(c("PKG_CPPFLAGS", "PKG_LIBS"), unset=NA_character_)
  on.exit({
    for (name in names(old)) {
      if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, setNames(list(old[[name]]), name))
    }
  }, add=TRUE)
  if (grepl("RcppArmadillo", paste(readLines(file, warn=FALSE), collapse="\n"), fixed=TRUE)) {
    include <- system.file("include", package="RcppArmadillo")
    stopifnot(nzchar(include))
    flag <- paste0("-I", shQuote(include))
    Sys.setenv(PKG_CPPFLAGS=paste(flag, if (is.na(old[["PKG_CPPFLAGS"]])) "" else old[["PKG_CPPFLAGS"]]))
    libraries <- "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"
    Sys.setenv(PKG_LIBS=paste(libraries, if (is.na(old[["PKG_LIBS"]])) "" else old[["PKG_LIBS"]]))
  }
  dir.create(cacheDir, recursive=TRUE, showWarnings=FALSE)
  Rcpp::sourceCpp(file, env=env, cacheDir=cacheDir, showOutput=showOutput, verbose=verbose)
}
