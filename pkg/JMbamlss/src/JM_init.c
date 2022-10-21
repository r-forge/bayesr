#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define USE_FC_LEN_T

SEXP survint(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP survint_re(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP)

static R_CallMethodDef callMethods[] = {
  {"survint", (DL_FUNC) &survint, 5},
  {"survint_re", (DL_FUNC) &survint, 6},
  {NULL, NULL, 0}
};

void R_init_sourcetools(DllInfo* info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

