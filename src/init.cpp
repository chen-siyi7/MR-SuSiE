#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern "C" SEXP _rcpp_module_boot_mr_susie_selector();
extern "C" SEXP _rcpp_module_boot_mr_susie_profile();
extern "C" SEXP _rcpp_module_boot_mr_susie_binary_gwas();

static const R_CallMethodDef call_entries[] = {
  {"_rcpp_module_boot_mr_susie_selector", (DL_FUNC) &_rcpp_module_boot_mr_susie_selector, 0},
  {"_rcpp_module_boot_mr_susie_profile", (DL_FUNC) &_rcpp_module_boot_mr_susie_profile, 0},
  {"_rcpp_module_boot_mr_susie_binary_gwas", (DL_FUNC) &_rcpp_module_boot_mr_susie_binary_gwas, 0},
  {NULL, NULL, 0}
};

extern "C" void attribute_visible R_init_MRSuSiE(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_entries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
