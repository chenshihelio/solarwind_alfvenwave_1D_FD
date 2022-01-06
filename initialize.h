#ifndef _INITIALIZE_
#define _INITIALIZE_

#include "macros.h"


extern real x0, x1, *xgrid, *dx;
extern int nx,arrSize;

extern int nvar,typeBC,
	ifRestart, restartIdx, ifResetTime, ifUpdateBC, ifViscosity;
extern real *uu, *dudt, *dudtRK;

extern real *FD_Coeff_CenterL_1st, *FD_Coeff_CenterR_1st,
	*FD_Coeff_CenterL_2nd, *FD_Coeff_CenterR_2nd,
	*FD_Coeff_Left_1st, *FD_Coeff_Right_1st,
	*FD_Coeff_Left_2nd, *FD_Coeff_Right_2nd;


extern real *fExp,*A,*dlnAdr,*Br,*wavelength,*Q_AH;
extern real adiabaticIdx, adiabaticIdx_E, F_AH;
extern real BC_N, BC_TI, BC_TE, BC_ZOUT, BC_BR, CW, LAMBDA0, CRflect;

extern real rho0,pi0,pe0,MS;

extern real cfl,Tmax,dtOut,timeScaleUpdateBC,viscosity;

// functions-----------------------
void initialize();
#endif
