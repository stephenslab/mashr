#include <stdlib.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

// Declarations for .Call entry points.
SEXP mashr_calc_lik_rcpp (SEXP b_matSEXP, SEXP s_matSEXP, SEXP v_matSEXP,
			  SEXP U_3dSEXP, SEXP logdSEXP);
SEXP mashr_calc_post_rcpp(SEXP b_matSEXP, SEXP s_matSEXP, SEXP v_matSEXP,
			  SEXP U_3dSEXP, SEXP posterior_weightsSEXP);

// See "Registering native routines" in "Writing R Extensions" manual
// for an explanation of what these lines of code do.
#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

const static R_CallMethodDef R_CallDef[] = {
  CALLDEF(mashr_calc_lik_rcpp,5),
  CALLDEF(mashr_calc_post_rcpp,5),
  {NULL, NULL, 0}
};

void attribute_visible R_init_varbvs(DllInfo *dll)
{
  R_registerRoutines(dll,NULL,R_CallDef,NULL,NULL);
  R_useDynamicSymbols(dll,FALSE);
  R_forceSymbols(dll,TRUE);
}
