#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void csp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _sdcTable_cpp_myPaste(SEXP, SEXP, SEXP);
extern SEXP _sdcTable_cpp_mySplit(SEXP, SEXP);
extern SEXP _sdcTable_cpp_splitByIndices(SEXP, SEXP);
extern SEXP _sdcTable_greedyMultDimSuppression(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"csp", (DL_FUNC) &csp, 22},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_sdcTable_cpp_myPaste",              (DL_FUNC) &_sdcTable_cpp_myPaste,              3},
    {"_sdcTable_cpp_mySplit",              (DL_FUNC) &_sdcTable_cpp_mySplit,              2},
    {"_sdcTable_cpp_splitByIndices",       (DL_FUNC) &_sdcTable_cpp_splitByIndices,       2},
    {"_sdcTable_greedyMultDimSuppression", (DL_FUNC) &_sdcTable_greedyMultDimSuppression, 5},
    {NULL, NULL, 0}
};

void R_init_sdcTable(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
