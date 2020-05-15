/* 
   Native symbol registration table for scuba package

*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

void HaldaneCalc(int *, int *, int *,
		 double *, double *, double *, double *, double *,
		 int *,
		 double *, double *);

static const R_CMethodDef CEntries[] = {
    {"HaldaneCalc",   (DL_FUNC) &HaldaneCalc,   11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {NULL, NULL, 0}
};

void R_init_scuba(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
  
