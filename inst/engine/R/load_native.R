mr_susie_load_native <- function(file, env, cacheDir, showOutput=FALSE, verbose=FALSE) {
  installed <- system.file('engine', 'src', basename(file), package='MRSuSiE')
  if (nzchar(installed) && file.exists(file) &&
      identical(digest::digest(file=file, algo='sha256'), digest::digest(file=installed, algo='sha256'))) {
    spec <- switch(basename(file),
      'selector.cpp'=list(module='mr_susie_selector', functions=c('native_single_effect', 'native_update_ss', 'native_likelihood_context', 'native_cached_likelihood')),
      'profile.cpp'=list(module='mr_susie_profile', functions='native_profile_quadratic'),
      'binary_gwas.cpp'=list(module='mr_susie_binary_gwas', functions='single_snp_logistic'),
      NULL)
    if (!is.null(spec)) {
      module <- Rcpp::Module(spec$module, PACKAGE='MRSuSiE', mustStart=TRUE)
      for (name in spec$functions) assign(name, do.call("$", list(module, name)), envir=env)
      return(invisible(NULL))
    }
  }
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
