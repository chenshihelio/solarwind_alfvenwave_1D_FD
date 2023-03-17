#include <math.h>
#include <stdlib.h>
#include "initialize.h"
#include "parameter.h"
#include "macros.h"

real* calcEigens(real *arr, real Br)
{
	// Note arr is the conserved fields
	real rho, u, pi, pe, energyIon, energyElectron;
	real csound, calfven;
	real cmin, cmax;

	real *eigens;
	eigens = (real*)malloc(sizeof(real) * NEigen);

	rho = arr[0];
	u = arr[1];
	calfven = Br / sqrt(MU0 * rho);

	pi = arr[2];
	pe = arr[3];

	csound = sqrt(adiabaticIdx * (pi + pe) / rho);

	cmin = MIN(fabs(calfven), fabs(csound));
	cmax = MAX(fabs(calfven), fabs(csound));

	eigens[0] = u - cmax;
	eigens[1] = u - cmin;
	eigens[2] = u;
	eigens[3] = u + cmin;
	eigens[4] = u + cmax;

	return eigens;
}

void setFirstPoint(real *uu, int nx, int nvar)
{
	uu[IDX(0, 0, nvar)] = rho0;
	uu[IDX(0, 2, nvar)] = pi0;
	uu[IDX(0, 3, nvar)] = pe0;
	uu[IDX(0, 4, nvar)] = 0.25 * rho0 * BC_ZOUT * BC_ZOUT;
}


real updateParameter(real Y0, real Y1, real timeSim, real timeScale)
{
	real Y;
	if (timeSim >= timeScale)
	{
		Y = Y1;
	}
	else
	{
		Y = Y0 + (Y1 - Y0) * timeSim / timeScale;
	}

	return Y;
}

void updateBC(real timeSim, real timeScale)
{
	// BC_ZOUT = updateParameter(0, 5e3, timeSim, timeScale);
	
	/* adiabaticIdx = updateParameter(1.0, 1.66667, timeSim, timeScale);
	adiabaticIdx_E = updateParameter(1.0, 1.66667, timeSim, timeScale);
	if(ifInhomoGammaI==1)
	{
		for (size_t i = 0; i < nx; i++)
		{
			real s = xgrid[i] / RS;
			gammaI_inhomo[i] = 1.0 + (adiabaticIdx - 1.0) * 
				(tanh((s-4)/0.5)+1)/2.0;
		}
	} */

	/* BC_TI = updateParameter(1e6,2e6,timeSim,timeScale);
	BC_TE = updateParameter(1e6,2e6,timeSim,timeScale);
	pi0 = BC_N * KB * BC_TI;
	pe0 = BC_N * KB * BC_TE; */

	CW = updateParameter(1.0, 0.95, timeSim, timeScale); 

	// MS = updateParameter(1.5E30,1.9884E30, timeSim, timeScale);  // 1.9884E30 kg

	// BC_N = updateParameter(5e14,1e14, timeSim, timeScale);
	// // Need to update the following quantities
	// rho0 = BC_N * (MP+ME);
	// pi0 = BC_N * KB * BC_TI;
	// pe0 = BC_N * KB * BC_TE; 

	// // Update Br, remember to re-calcualte the array of Br
	// BC_BR = updateParameter(1e-4, 3e-5, timeSim, timeScale);
	// for (size_t i = 0; i < nx; i++)
	// {
	// 	Br[i] = (A[0] / A[i]) * BC_BR;
	// }

	/* // Update lambda, remember to calculate the whole array
	LAMBDA0 = updateParameter(3E7, 6E7, timeSim, timeScale);
	for (size_t i = 0; i < nx; i++)
	{
		wavelength[i] = LAMBDA0 * sqrt(A[i] / A[0]);
	} */

    // // Update Q_AH
	// F_AH = updateParameter(2e-7,5e-7,timeSim,timeScale);
	// for(int i=0;i<nx;i++)
	// {
	// 	Q_AH[i] = F_AH * (A[0]/A[i]) * exp(-(xgrid[i]/RS-1)/scaleHeight_AH);
	// }

	// C_AH_electron = updateParameter(1.0,0.4,timeSim,timeScale);


	// update CQeColless
	// CQeColless = updateParameter(0,1,timeSim, timeScale);

	// // update C_ratio_ion_colless
	// C_ratio_ion_colless = updateParameter(0.5,0,timeSim, timeScale);
}
