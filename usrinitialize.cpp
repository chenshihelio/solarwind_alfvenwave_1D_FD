#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros.h"
#include "initialize.h"
#include "parameter.h"

void initFields(real *arr,  size_t nx,  int nvar)
{
    real s, amp;
    real u,rho,pi,pe,zout,zin;
    for(size_t i=0;i<nx;i++)
    {
        s = xgrid[i]/x0;

		// amp = exp(- 0.01 * (s-1) );
		amp = pow(s,-1);

        /*rho = rho0 * amp;
        u = 100;
        pi = pi0 * amp;
        pe = pe0 * amp;*/

		u = 300e3;
		// u = 150e3 * (1+tanh((s - 20) / 4));
		rho = rho0 * A[0] / A[i];
		pi = pi0 * pow(A[0] / A[i], adiabaticIdx);
		pe = pe0 * pow(A[0] / A[i], adiabaticIdx);

        arr[IDX(i,0,NVAR)] =  rho;

        arr[IDX(i,1,NVAR)] = u;

        arr[IDX(i,2,NVAR)] = pi;
        arr[IDX(i,3,NVAR)] = pe;

        zout = BC_ZOUT;
        zin = 0;
        arr[IDX(i,4,NVAR)] = 0.25 * rho * zout * zout;
        arr[IDX(i,5,NVAR)] = 0.25 * rho * zin * zin;
    }    
}
