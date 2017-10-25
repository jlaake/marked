#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP makeZtildeIdx(SEXP, SEXP);
extern SEXP sampleZ(SEXP, SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(cjs)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cjs1tlgam)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cjs1tlp)(void *, void *, void *, void *, void *);
extern void F77_NAME(cjs2tlgam)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cjs2tlp)(void *, void *, void *, void *, void *);
extern void F77_NAME(cjsgam)(void *, void *, void *, void *, void *);
extern void F77_NAME(cjsp)(void *, void *, void *, void *, void *);
extern void F77_NAME(hmmlike)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ms2gam)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msgam)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(msp)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mvmsp)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ums2p)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(umsp)(void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"makeZtildeIdx", (DL_FUNC) &makeZtildeIdx, 2},
    {"sampleZ",       (DL_FUNC) &sampleZ,       4},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"cjs",       (DL_FUNC) &F77_NAME(cjs),       16},
    {"cjs1tlgam", (DL_FUNC) &F77_NAME(cjs1tlgam),  6},
    {"cjs1tlp",   (DL_FUNC) &F77_NAME(cjs1tlp),    5},
    {"cjs2tlgam", (DL_FUNC) &F77_NAME(cjs2tlgam),  6},
    {"cjs2tlp",   (DL_FUNC) &F77_NAME(cjs2tlp),    5},
    {"cjsgam",    (DL_FUNC) &F77_NAME(cjsgam),     5},
    {"cjsp",      (DL_FUNC) &F77_NAME(cjsp),       5},
    {"hmmlike",   (DL_FUNC) &F77_NAME(hmmlike),   11},
    {"ms2gam",    (DL_FUNC) &F77_NAME(ms2gam),    10},
    {"msgam",     (DL_FUNC) &F77_NAME(msgam),      7},
    {"msp",       (DL_FUNC) &F77_NAME(msp),        6},
    {"mvmsp",     (DL_FUNC) &F77_NAME(mvmsp),     13},
    {"ums2p",     (DL_FUNC) &F77_NAME(ums2p),      8},
    {"umsp",      (DL_FUNC) &F77_NAME(umsp),       7},
    {NULL, NULL, 0}
};

void R_init_marked(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
