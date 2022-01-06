#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "macros.h"
#include "parameter.h"
#include "output.h"
#include "ios_support.h"
#include <omp.h>

#define MAXBUFFSIZE 200

real x0, x1, *xgrid, *dx;
int nx, arrSize;

int nvar=NVAR,typeBC=0,ifRestart = 0, restartIdx = 0, ifResetTime,
	ifUpdateBC = 0, ifViscosity = 0;

real *uu, *dudt, *dudtRK;

real *FD_Coeff_CenterL_1st, *FD_Coeff_CenterR_1st,
	*FD_Coeff_CenterL_2nd, *FD_Coeff_CenterR_2nd,
	*FD_Coeff_Left_1st, *FD_Coeff_Right_1st,
	*FD_Coeff_Left_2nd, *FD_Coeff_Right_2nd;

real BC_N, BC_TI, BC_TE, BC_BR, BC_ZOUT, CW, LAMBDA0,CRflect;
real *A,*fExp,*dlnAdr,*Br,*wavelength;
real *Q_AH, F_AH; // ad-hoc heating term
real rho0,pi0,pe0,MS;

real cfl=0.5,Tmax=1.0,dtOut = 0.1, timeScaleUpdateBC = 1E4,
	viscosity = 1.0;
real adiabaticIdx, adiabaticIdx_E;


// functions-------------------
//real* initGrid( real x0,  real x1,
//        size_t nx,  int typeBC, real **ptr_xCI, 
//        real **ptr_xgridSize);

real* initArrays(size_t nx,  int nvar);


real stretch(int i, int N, real asymp, int offset);
int createGrid(real **ptr_xarr, real h0, real h1, real x0, real x1);

void readInput(const char *fileName);

real *calc_coeff_right_1st(real h0, real h1, real h2);
real *calc_coeff_left_1st(real h0, real h1, real h2);
real *calc_coeff_centerL_1st(real hm2, real hm1, real hp1);
real *calc_coeff_centerR_1st(real hm1, real hp1, real hp2);
real *calc_coeff_right_2nd(real h0, real h1, real h2);
real *calc_coeff_left_2nd(real h0, real h1, real h2);
real *calc_coeff_centerL_2nd(real hm2, real hm1, real hp1);
real *calc_coeff_centerR_2nd(real hm1, real hp1, real hp2);


void initialize()
{
    readInput("./input");

    printf("1D Finite-Difference nonuniform-grid code for Solar Wind + Alfven wave model\n");
    printf("---------------------------\n");
    printf("Basic configurations:\n");
    printf("    Adiabatic Index for Ion = %10.5f\n", adiabaticIdx);
	printf("    Adiabatic Index for Electron = %10.5f\n", adiabaticIdx_E);
	printf("    Mass of the Sun = %10.4e kg\n", MSUN);
    printf("---------------------------\n");
    printf("Boundary conditions at r0:\n");
    printf("    Density = %10.4e m^{-3}\n", BC_N);
    printf("    Ion Temperature = %10.4e K\n", BC_TI);
    printf("    Electron Temperature = %10.4e K\n", BC_TE);
    printf("    Magnetic field = %10.4e T\n", BC_BR);
    printf("    Zout = %10.4e m/s\n", BC_ZOUT);
    printf("    LAMBDA0 (perpendicular wavelength) = %10.4e m\n",LAMBDA0);
    printf("Other parameters:\n");
    printf("    CW (how much energy dissipation goes to ion) = %10.4f\n",CW);
	printf("    F_AH (ad-hoc heating strength) = %10.4e J/m^3/s\n", F_AH);
	printf("---------------------------\n");
	printf("    If update the boundary conditions with time?: %1d\n", ifUpdateBC);
	printf("    timescale for boundary condition update = %10.4e s\n", timeScaleUpdateBC);
	printf("---------------------------\n");
	printf("    If implement an explicity viscosity?: %1d\n", ifViscosity);
	printf("    Viscosity = %10.4e m^2/s\n", viscosity);

	MS = MSUN;

	// OMP initialize
	printf("Maximum number of processors: %d\n", omp_get_max_threads());
	// int NPROC = omp_get_max_threads();
	int NPROC = 1;
	printf("Number of processors to use: %d\n", NPROC);
	omp_set_num_threads(NPROC);
	
    // initialize grid()
	x0 = 1.1;
	x1 = 215;
	real h0 = 1e-3;
	real h1 = 0.1;
	// real h0 = 0.02;
	// real h1 = 0.02;
	nx = createGrid(&xgrid, h0, h1, x0, x1);
	printf("    x0 = %.3f Rs, x1 = %.3f Rs, nx = %8d\n", x0, xgrid[nx - 1], nx);
	for (int i = 0; i < nx; i++)
	{
		xgrid[i] = xgrid[i] * RS;
	}
	dx = initArrays(nx - 1, 1);
	for (int i = 0; i < nx - 1; i++)
	{
		dx[i] = xgrid[i + 1] - xgrid[i];
	}
	x0 = x0 * RS;
	x1 = x1 * RS;
	h0 = h0 * RS;
	h1 = h1 * RS;

    // initialize arrays
    nvar = NVAR;
    
    uu = initArrays(nx, nvar);
    dudt = initArrays(nx, nvar);
    dudtRK = initArrays(nx, nvar);

    arrSize = nx * nvar;

    // calculate expansion factor
    fExp = initArrays(nx,1);
    A = initArrays(nx,1);
    dlnAdr = initArrays(nx,1);
    wavelength = initArrays(nx,1);

    real f1,s,df;
    f1 = 1 - (FMAX - 1) * exp((1-REXP)/SIGMA_EXP);
    for(size_t i=0;i<nx;i++)
    {
        s = xgrid[i] / RS;
         fExp[i] = (FMAX + f1 * exp(-(s-REXP)/SIGMA_EXP)) / 
             (1+ exp(-(s-REXP)/SIGMA_EXP));

         df = (FMAX-f1)/(2*SIGMA_EXP)/
             (1+cosh((s-REXP)/SIGMA_EXP)) / RS;

		/*fExp[i] = 1.;
		df = 0.;*/

		A[i] = fExp[i] * s * s;
        dlnAdr[i] = 2./xgrid[i] + df/fExp[i];
		/*A[i] = 1;
		dlnAdr[i] = 0;*/

        wavelength[i] = LAMBDA0 * sqrt(A[i]/A[0]);
    }


    // calculate Br
    Br = initArrays(nx,1);
    for(size_t i=0;i<nx;i++)
    {
        Br[i] = (A[0]/A[i]) * BC_BR;
    }

	// calculate ad-hoc heating term
	Q_AH = initArrays(nx, 1);
	for (size_t i = 0; i < nx; i++)
	{
		s = xgrid[i] / RS;
		Q_AH[i] = F_AH * (A[0] / A[i]) * exp(-(s-1)/scaleHeight_AH);
	}

    // calculate some useful quantities
    rho0 = BC_N * (MP+ME);
    pi0 = BC_N * KB * BC_TI;
    pe0 = BC_N * KB * BC_TE;


	// finite-difference schemes
	FD_Coeff_CenterL_1st = initArrays(nx, NPOINTFD);
	FD_Coeff_CenterL_2nd = initArrays(nx, NPOINTFD);
	FD_Coeff_CenterR_1st = initArrays(nx, NPOINTFD);
	FD_Coeff_CenterR_2nd = initArrays(nx, NPOINTFD);
	FD_Coeff_Left_1st = initArrays(nx, NPOINTFD);
	FD_Coeff_Left_2nd = initArrays(nx, NPOINTFD);
	FD_Coeff_Right_1st = initArrays(nx, NPOINTFD);
	FD_Coeff_Right_2nd = initArrays(nx, NPOINTFD);

	real *coef_tmp;
	for (int i = 0; i < nx; i++)
	{
		// Right scheme--------------------
		if (i <= nx - NPOINTFD)
		{
			coef_tmp = calc_coeff_right_1st(dx[i], dx[i + 1], dx[i + 2]);
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_Right_1st[i*NPOINTFD + j] = coef_tmp[j];
			}
			free(coef_tmp);

			coef_tmp = calc_coeff_right_2nd(dx[i], dx[i + 1], dx[i + 2]);
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_Right_2nd[i*NPOINTFD + j] = coef_tmp[j];
			}
			free(coef_tmp);
		}
		else
		{
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_Right_1st[i*NPOINTFD + j] = NAN;
				FD_Coeff_Right_2nd[i*NPOINTFD + j] = NAN;
			}

		}

		// Left Scheme
		if (i >= (NPOINTFD - 1))
		{
			coef_tmp = calc_coeff_left_1st(dx[i - 1], dx[i - 2], dx[i - 3]);
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_Left_1st[i*NPOINTFD + j] = coef_tmp[j];
			}
			free(coef_tmp);

			coef_tmp = calc_coeff_left_2nd(dx[i - 1], dx[i - 2], dx[i - 3]);
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_Left_2nd[i*NPOINTFD + j] = coef_tmp[j];
			}
			free(coef_tmp);
		}
		else
		{
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_Left_1st[i*NPOINTFD + j] = NAN;
				FD_Coeff_Left_2nd[i*NPOINTFD + j] = NAN;
			}

		}

		// Center-Left Scheme
		if (i >= 2 && i < nx - 1)
		{
			coef_tmp = calc_coeff_centerL_1st(dx[i - 2], dx[i - 1], dx[i]);
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_CenterL_1st[i*NPOINTFD + j] = coef_tmp[j];
			}
			free(coef_tmp);

			coef_tmp = calc_coeff_centerL_2nd(dx[i - 2], dx[i - 1], dx[i]);
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_CenterL_2nd[i*NPOINTFD + j] = coef_tmp[j];
			}
			free(coef_tmp);
		}
		else
		{
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_CenterL_1st[i*NPOINTFD + j] = NAN;
				FD_Coeff_CenterL_2nd[i*NPOINTFD + j] = NAN;
			}
		}


		// Center-Right Scheme
		if (i >= 1 && i < nx - 2)
		{
			coef_tmp = calc_coeff_centerR_1st(dx[i - 1], dx[i], dx[i + 1]);
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_CenterR_1st[i*NPOINTFD + j] = coef_tmp[j];
			}
			free(coef_tmp);

			coef_tmp = calc_coeff_centerR_2nd(dx[i - 1], dx[i], dx[i + 1]);
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_CenterR_2nd[i*NPOINTFD + j] = coef_tmp[j];
			}
			free(coef_tmp);
		}
		else
		{
			for (int j = 0; j < NPOINTFD; j++)
			{
				FD_Coeff_CenterR_1st[i*NPOINTFD + j] = NAN;
				FD_Coeff_CenterR_2nd[i*NPOINTFD + j] = NAN;
			}
		}
	}

    return;
}



real stretch(int i, int N, real asymp, int offset)
{
	real den = 1 + tanh(real(offset) / real(N));
	real num = tanh(real(i - offset) / real(N))
		+ tanh(real(offset) / real(N));

	return 1 + (asymp - 1) * (num / den);
}

int createGrid(real **ptr_xarr, real h0, real h1, real x0, real x1)
{
	int NMAX = int((x1 - x0) / h0);
	real *xtmp;
	xtmp = (real*)malloc(sizeof(real) * NMAX);
	xtmp[0] = x0;
	xtmp[1] = x0 + h0;

	real x = xtmp[1];
	int i = 1;
	real s, dx;

	while (x < x1)
	{
		s = stretch(i, 200, h1 / h0, 800);
		dx = s * h0;
		x = x + dx;
		xtmp[i + 1] = x;
		i++;
	}

	int nx = i + 1;
	*ptr_xarr = (real*)malloc(sizeof(real) * nx);
	for (i = 0; i < nx; i++)
	{
		(*ptr_xarr)[i] = xtmp[i];
	}

	return nx;
}



char *strsep(char **stringp, const char *delim) {
	char *rv = *stringp;
	if (rv) {
		*stringp += strcspn(*stringp, delim);
		if (**stringp)
			*(*stringp)++ = '\0';
		else
			*stringp = 0;
	}
	return rv;
}

char *readLine(FILE *fp)
{
    char fileLine[MAXBUFFSIZE];
    char *token, *buf;

    fgets(fileLine, MAXBUFFSIZE, fp);

    buf = fileLine;
    //printf("buf=%s\n",buf);
    token = strsep(&buf, ";");

    // do not understand how the original string is changed...
    // printf("buf=%s\n",buf);
    // printf("fileLine=%s\n",fileLine); 

    return token;
}

void readInput(const char *fileName)
{
    // FILE *fpRead=fopen(fileName, "r");
	FILE *fpRead;
	errno_t err = fopen_s(&fpRead, fileName, "r");

    char *valueRead;

	valueRead = readLine(fpRead);
	ifRestart = (int)atof(valueRead);

	valueRead = readLine(fpRead);
	restartIdx = (int)atof(valueRead);

	valueRead = readLine(fpRead);
	ifResetTime = (int)atof(valueRead);

    valueRead = readLine(fpRead);
	typeBC = (int)atof(valueRead);

    valueRead = readLine(fpRead);
    cfl = atof(valueRead);

    valueRead = readLine(fpRead);
    Tmax = atof(valueRead);

    valueRead = readLine(fpRead);
    dtOut = atof(valueRead);

    valueRead = readLine(fpRead);
    adiabaticIdx = atof(valueRead);

	valueRead = readLine(fpRead);
    adiabaticIdx_E = atof(valueRead);

    valueRead = readLine(fpRead);
    BC_N = atof(valueRead);

    valueRead = readLine(fpRead);
    BC_TI = atof(valueRead);

    valueRead = readLine(fpRead);
    BC_TE = atof(valueRead);

    valueRead = readLine(fpRead);
    BC_BR = atof(valueRead);

    valueRead = readLine(fpRead);
    BC_ZOUT = atof(valueRead);

    valueRead = readLine(fpRead);
    LAMBDA0 = atof(valueRead);

    valueRead = readLine(fpRead);
    CW = atof(valueRead);

	valueRead = readLine(fpRead);
	CRflect = atof(valueRead);

	valueRead = readLine(fpRead);
	ifUpdateBC = (int)atof(valueRead);

	valueRead = readLine(fpRead);
	timeScaleUpdateBC = atof(valueRead);

	valueRead = readLine(fpRead);
	ifViscosity = (int)atof(valueRead);

	valueRead = readLine(fpRead);
	viscosity = atof(valueRead);

	valueRead = readLine(fpRead);
	F_AH = atof(valueRead);

    int closeStatus = fclose(fpRead);
}

real* initArrays( size_t nx,  int nvar)
{
    real *arr;
    arr = (real *)malloc(sizeof(real) * nx * nvar);
    memset(arr,0,sizeof(real) * nx * nvar);

    return arr;
}




// Initialize derivative schemes----------------------------
////////////////
inline real cc_right_1st(real c0, real c1, real c2)
{
	real f1 = c0 / c1;
	real f2 = c0 / c2;
	real r1 = (c0 - c2) / (c2 - c1);
	real r2 = (c1 - c0) / (c2 - c1);
	return 1 / c0 / (1 + f1 * r1 + f2 * r2);
}

real *calc_coeff_right_1st(real h0, real h1, real h2)
{
	real H0 = h0;
	real H1 = h0 + h1;
	real H2 = H1 + h2;

	real coef1 = cc_right_1st(H0, H1, H2);
	real coef2 = cc_right_1st(H1, H2, H0);
	real coef3 = cc_right_1st(H2, H0, H1);
	real coef0 = -(coef1 + coef2 + coef3);

	real *coefs;
	coefs = (real*)malloc(sizeof(real)*NPOINTFD);
	coefs[0] = coef0;
	coefs[1] = coef1;
	coefs[2] = coef2;
	coefs[3] = coef3;

	return coefs;
}

real *calc_coeff_left_1st(real h0, real h1, real h2)
{
	real *coefs = calc_coeff_right_1st(h0, h1, h2);

	real *coefs_left = (real*)malloc(sizeof(real) * NPOINTFD);
	for (int i = 0; i < NPOINTFD; i++)
	{
		coefs_left[i] = -coefs[i];
	}

	return coefs_left;
}


////////////////
inline real cc_right_2nd(real c0, real c1, real c2)
{
	real f1 = c1 / c0;
	real f2 = c2 / c0;

	real r1 = ((c0 - c2) / (c1 - c2)) * ((c0 + c2) / (c1 + c2));
	real r2 = ((c1 - c0) / (c1 - c2)) * ((c1 + c0) / (c1 + c2));

	return 2 / c0 / c0 / (1 - f1 * r1 - f2 * r2);
}

real *calc_coeff_right_2nd(real h0, real h1, real h2)
{
	real H0 = h0;
	real H1 = h0 + h1;
	real H2 = H1 + h2;

	real coef1 = cc_right_2nd(H0, H1, H2);
	real coef2 = cc_right_2nd(H1, H2, H0);
	real coef3 = cc_right_2nd(H2, H0, H1);
	real coef0 = -(coef1 + coef2 + coef3);

	real *coefs;
	coefs = (real*)malloc(sizeof(real)*NPOINTFD);
	coefs[0] = coef0;
	coefs[1] = coef1;
	coefs[2] = coef2;
	coefs[3] = coef3;

	return coefs;
}


real *calc_coeff_left_2nd(real h0, real h1, real h2)
{
	real *coefs = calc_coeff_right_2nd(h0, h1, h2);
	return coefs;
}

///////////////////

real *calc_coeff_centerL_1st(real hm2, real hm1, real hp1)
{
	real dm2 = hm2 + hm1;
	real dm1 = hm1;
	real dp1 = hp1;

	real coef_m2 = dm1 * dp1 / (dm2*(dm2 + dp1)*hm2);
	real coef_m1 = -dm2 * dp1 / (dm1*hm2*(dm1 + dp1));
	real coef_p1 = dm2 * dm1 / (dp1*(dm1 + dp1)*(dm2 + dp1));
	real coef_0 = -(coef_m2 + coef_m1 + coef_p1);

	real *coefs;
	coefs = (real*)malloc(sizeof(real)*NPOINTFD);

	coefs[0] = coef_m2;
	coefs[1] = coef_m1;
	coefs[2] = coef_0;
	coefs[3] = coef_p1;

	return coefs;
}

real *calc_coeff_centerR_1st(real hm1, real hp1, real hp2)
{
	real *coefs = calc_coeff_centerL_1st(hp2, hp1, hm1);

	real *coefs_R = (real*)malloc(sizeof(real) * NPOINTFD);
	for (int i = 0; i < NPOINTFD; i++)
	{
		coefs_R[i] = -coefs[NPOINTFD - 1 - i];
	}

	return coefs_R;
}


real *calc_coeff_centerL_2nd(real hm2, real hm1, real hp1)
{
	real dm2 = hm2 + hm1;
	real dm1 = hm1;
	real dp1 = hp1;

	real coef_m2 = 2 * (dp1 - dm1) / (dm2*(dm2 + dp1)*(dm2 - dm1));
	real coef_m1 = 2 * (dm2 - dp1) / (dm1*(dm2 - dm1)*(dm1 + dp1));
	real coef_p1 = 2 * (dm2 + dm1) / (dp1*(dm1 + dp1)*(dm2 + dp1));
	real coef_0 = -(coef_m2 + coef_m1 + coef_p1);

	real *coefs;
	coefs = (real*)malloc(sizeof(real)*NPOINTFD);

	coefs[0] = coef_m2;
	coefs[1] = coef_m1;
	coefs[2] = coef_0;
	coefs[3] = coef_p1;

	return coefs;
}

real *calc_coeff_centerR_2nd(real hm1, real hp1, real hp2)
{
	real *coefs = calc_coeff_centerL_2nd(hp2, hp1, hm1);

	real *coefs_R = (real*)malloc(sizeof(real) * NPOINTFD);
	for (int i = 0; i < NPOINTFD; i++)
	{
		coefs_R[i] = coefs[NPOINTFD - 1 - i];
	}

	return coefs_R;
}
/////////////////////////////////////
