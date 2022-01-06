#ifndef _CALCDERIVS_
#define _CALCDERIVS_

#include "macros.h"

real* derivX(real *uu, real *xgrid, real *Br, int nx, int nvar);
real* derivXX(real *uu, real *xgrid, real *Br, int nx, int nvar);
#endif // !_CALCDERIVS_

