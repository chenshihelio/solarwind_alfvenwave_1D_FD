#include "initialize.h"
#include "macros.h"
#include "parameter.h"
#include <math.h>

void calcRHS(real *dudt, real *uu, real *dudx, real *d2udx2, 
	real *xgrid, real *Br, real *lambda, real *A, real *dlnA, 
	real *Q_AH, real Cwave, real Creflect, int nx, int nvar)
{
	real rho, ux, pi, pe, eo, ei, n, VA, r, lamb,
		divu, divVA,drho,dux,dpi,dpe,deo,dei,dlnAdr;
	real Qo, Qi, Ro, Ri, dissRateIon, dissRateElectron,
		Q_AH_i, Q_AH_e;

	for (int i = 0; i < nx; i++)
	{
		r = xgrid[i];
		lamb = lambda[i];

		rho = uu[IDX(i, 0, nvar)];
		ux = uu[IDX(i, 1, nvar)];
		pi = uu[IDX(i, 2, nvar)];
		pe = uu[IDX(i, 3, nvar)];
		eo = uu[IDX(i, 4, nvar)];
		ei = uu[IDX(i, 5, nvar)];

		drho = dudx[IDX(i, 0, nvar)];
		dux = dudx[IDX(i, 1, nvar)];
		dpi = dudx[IDX(i, 2, nvar)];
		dpe = dudx[IDX(i, 3, nvar)];
		deo = dudx[IDX(i, 4, nvar)];
		dei = dudx[IDX(i, 5, nvar)];

		dlnAdr = dlnA[i];

		n = rho / (MP + ME);

		VA = Br[i] / sqrt(MU0 * rho);

		divu = dux + ux * dlnAdr;
		divVA = -0.5*VA * (drho / rho);

		if ((ei < 0) || (eo < 0))
		{
			Qi = 0;
			Qo = 0;
		}
		else
		{
			Qo = -sqrt(ei)*eo / sqrt(rho) / lamb;
			Qi = -sqrt(eo)*ei / sqrt(rho) / lamb;
		}

		dissRateIon = Cwave * (fabs(Qo) + fabs(Qi));
		dissRateElectron = (1 - Cwave) * (fabs(Qo) + fabs(Qi));
		Ro = fabs(drho / 2. / rho * (ux - VA) * ei * Creflect);
		Ri = fabs(drho / 2. / rho * (ux + VA) * eo * Creflect);


		Q_AH_i = Q_AH[i];
		Q_AH_e = Q_AH[i];

		dudt[IDX(i, 0, nvar)] = -(ux*drho + rho * divu);
		dudt[IDX(i, 1, nvar)] = -ux * dux - (dpi + dpe + 0.5*(deo+dei)) / rho - G * MS / r / r;
		dudt[IDX(i, 2, nvar)] = -ux * dpi - adiabaticIdx * divu * pi + 
			(adiabaticIdx - 1) * (dissRateIon + Q_AH_i);
		dudt[IDX(i, 3, nvar)] = -ux * dpe - adiabaticIdx_E * divu * pe + 
			(adiabaticIdx_E - 1) * (dissRateElectron + Q_AH_e);
		dudt[IDX(i, 4, nvar)] =  -(divu + divVA)*eo - (ux + VA)*deo - 0.5*divu*eo + Qo + Ro;
		dudt[IDX(i, 5, nvar)] =  -(divu - divVA)*ei - (ux - VA)*dei - 0.5*divu*ei + Qi + Ri;

		if (ifViscosity == 1)
		{
			dudt[IDX(i, 1, nvar)] = dudt[IDX(i, 1, nvar)] + viscosity *
				(d2udx2[IDX(i, 1, nvar)]); //  + dux * dlnAdr
		}

	}

	dudt[IDX(0, 0, nvar)] = 0; // rho
	dudt[IDX(0, 2, nvar)] = 0; // Pi
	dudt[IDX(0, 3, nvar)] = 0; // Pe
	dudt[IDX(0, 4, nvar)] = 0; // Eout
}