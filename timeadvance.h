#ifndef _TIMEADVANCE_
#define _TIMEADVANCE_

#include "macros.h"

extern real cc1[3], dd1[3];

void rkInit(real *cc, real *dd, real dt, real *arrRK,  size_t arrSize);

void rkExecute( int iRK,  real *cc,  real *dd, 
    real *arr, real *arrdt, real *arrdtRK,  size_t arrSize);

real calcDt(real dt, real *arr, real *Br, real *dx,
	size_t nx, int nvar, real cfl);

#endif