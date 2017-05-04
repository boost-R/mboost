
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP ngradientCoxPLik(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_mcumsum(SEXP);
extern SEXP R_trace_gamboost(SEXP, SEXP, SEXP);
extern SEXP R_trace_glmboost(SEXP, SEXP, SEXP);
extern SEXP R_ysum(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ngradientCoxPLik", (DL_FUNC) &ngradientCoxPLik, 4},
    {"R_mcumsum",        (DL_FUNC) &R_mcumsum,        1},
    {"R_trace_gamboost", (DL_FUNC) &R_trace_gamboost, 3},
    {"R_trace_glmboost", (DL_FUNC) &R_trace_glmboost, 3},
    {"R_ysum",           (DL_FUNC) &R_ysum,           2},
    {NULL, NULL, 0}
};

void R_init_mboost(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
