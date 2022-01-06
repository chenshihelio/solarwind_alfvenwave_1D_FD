#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "macros.h"
#include "support.h"
#include <omp.h>

#define CC10 8.0/15.0
#define CC20 5.0/12.0
#define CC30 0.75

#define DD20 -17./60.
#define DD30 -5./12.


real cc1[3], dd1[3];


void rkInit(real *cc, real *dd, real dt, real *arrRK,  size_t arrSize)
{
    cc[0] = CC10*dt ; dd[0] = 0.0;
    cc[1] = CC20*dt ; dd[1] = DD20*dt;
    cc[2] = CC30*dt ; dd[2] = DD30*dt;
    
    memset(arrRK, 0, arrSize * sizeof(real));
}

void rkExecute( int iRK,  real *cc,  real *dd, 
    real *arr, real *arrdt, real *arrdtRK,  size_t arrSize)
{
	long i;
	#pragma omp parallel for default(none) \
        shared(arr,cc,arrSize,dd,arrdt,arrdtRK,iRK)\
		private(i)
    for(i=0; i<arrSize; i++)
    {
        arr[i] = arr[i] + cc[iRK] * arrdt[i] + dd[iRK] * arrdtRK[i];
        arrdtRK[i] = arrdt[i]; 
    }
}

real calcDt( real dt, real *arr, real *Br, real *dx, 
     size_t nx,  int nvar, real cfl)
{
    real *waveSpeeds;
    real maxSpeed, dtx, dtMin, dtMinIter;

    dtMin = -1;
    for(size_t i=0;i<nx;i++)
    {
        waveSpeeds = calcEigens(&arr[IDX(i,0,nvar)], Br[i]);
        if(waveSpeeds[2]>0)
        {
            maxSpeed = fabs(waveSpeeds[NEigen-1]); // use fabs
        }
        else 
        {
            maxSpeed = fabs(waveSpeeds[0]); // use fabs
        }
		
		if (i == nx - 1)
		{
			dtx = dx[i-1] / maxSpeed;
		}
		else
		{
			dtx = dx[i] / maxSpeed;
		}
        
        dtMinIter = dtx;

        if ( dtMin < 0 || dtMinIter < dtMin )
        {
            dtMin = dtMinIter;
        }

		free(waveSpeeds);
    }

    dtMin = cfl * dtMin;

    real dtNew;
    if (dt<0.98*dtMin || dt>1.02*dtMin)
    {
        dtNew = dtMin;
    }
    else
    {
        dtNew = dt;
    }

    return dtNew;
}
