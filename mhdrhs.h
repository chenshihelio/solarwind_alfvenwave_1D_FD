#ifndef _MHDRHS_
#define _MHDRHS_

#include "macros.h"

void calcRHS(real *dudt, real *uu, real *dudx, real *d2udx2,
	real *xgrid, real *Br, real *lambda, real *A, real *dlnA,
	real *Q_AH,real Cwave, real Creflect, int nx, int nvar);

#endif