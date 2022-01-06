#ifndef _SUPPORT_
#define _SUPPORT_

#include "macros.h"

real* calcEigens(real *arr, real Br);

void setFirstPoint(real *uu, int nx, int nvar);

void updateBC(real timeSim, real timeScale);

#endif