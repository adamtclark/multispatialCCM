#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void CCM_bootstrap(double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
extern void SSR_predict_boot(double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, int *);

static const R_CMethodDef CEntries[] = {
    {"CCM_bootstrap",    (DL_FUNC) &CCM_bootstrap,    12},
    {"SSR_predict_boot", (DL_FUNC) &SSR_predict_boot, 12},
    {NULL, NULL, 0}
};

void R_init_multispatialCCM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
