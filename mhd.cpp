#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include "macros.h"
#include "initialize.h"
#include "usrinitialize.h"
#include "calcDerivs.h"
#include "mhdrhs.h"
#include "timeadvance.h"
#include "output.h"
#include "support.h"
#include "ios_support.h"
#include "filter.h"

using namespace std;

int checkNan(real *arr);
void writeLog(real timeSim, real dt, real diffTime, size_t iStep);

int main(int argc, char **argv)
{
    size_t iStep=0, dstepCheckNan = 20, dstepLog = 1000;
    real timeSim = 0, dt, tOut = 0, diffTime;
    time_t tStart, tRunTime;
    int iOutput = 0;
	real *dudx, *d2udx2;

    time(&tStart);

    initialize();

    write1DArray(nx, xgrid, "xgrid.dat");
	write1DArray(nx, Br, "Br.dat");

	write1DArray(nx, fExp, "fExp.dat");
    
	write1DArray(nx, A, "A.dat");
	write1DArray(nx, dlnAdr, "dlnAdr.dat");


    // define initial conditions
	if (ifRestart == 0)
	{
		printf("Calling usr-defined initialization function\n");
		initFields(uu, nx, nvar);
		timeSim = 0;
		iOutput = 0;
	}
	else
	{
		printf("Reading restart.dat, will restart from index %04d\n", 
				restartIdx);
		timeSim = readArray("restart.dat", uu, arrSize);
		iOutput = restartIdx;
		if (ifResetTime == 1) timeSim = 0;
		tOut = timeSim;
	}
    
	// // mannually add zin
	// addInitialZin(0.2, uu, nx, nvar);

    // Output the main array at t = 0
    writeOutput(iOutput,timeSim,uu,arrSize,nx*nvar);
    iOutput++;
    tOut += dtOut;

	//// debug use--------------
	//real *dudx_tmp = derivX(uu, xgrid, Br, nx, nvar);
	//writeOutput(999, timeSim, dudx_tmp, arrSize, nx*nvar);
	//free(dudx_tmp);

    // calculate dt at initial step
    dt = -1;
    dt = calcDt(dt, uu, Br, dx, nx, nvar, cfl);
    printf("dt = %e\n", dt);

    // time evolution
    while(timeSim<Tmax)
    {
        // calculate dt
        dt = calcDt(dt, uu, Br, dx, nx, nvar, cfl);
        rkInit(&cc1[0], &dd1[0], dt, dudtRK, arrSize);

        /* exit conditions */
        // check NaNs
        if ((iStep>0) && ((iStep%dstepCheckNan)== 0))
        {
            int isNanAll = checkNan(uu);
            if (isNanAll>0)
            {
                printf("Nan values encountered at t = %10.4f!!!\n", timeSim);
                break;
            }
        }

        // check dt 
        if (dt < 1e-8)
        {
            printf("****************************************************\n");
            printf("\n");
            printf("  STOP: dt < 1e-8\n");
            printf("  time = %10.4f, dt = %12.4e\n", timeSim, dt);
            printf("  iteration = %14ld\n", iStep);
            printf("\n");
            printf("****************************************************\n");
            printf("  Output at time = %10.4f, dt = %12.4e\n", timeSim, dt);

            writeOutput(iOutput,timeSim,uu,arrSize,nx*nvar);
            break;
        }


        /* Save data  */ 
        // check timeSim and save data
        // output main array
        if (timeSim >= tOut)
        {
            writeOutput(iOutput,timeSim,uu,arrSize,nx*nvar);
            iOutput++;
            tOut += dtOut;
        }

		// output log
		if (iStep % dstepLog == 0)
		{
			time(&tRunTime);
			diffTime = difftime(tRunTime, tStart);

			writeLog(timeSim, dt, diffTime, iStep);
		}

       

        // evolution module ------------
		if (ifUpdateBC == 1)
		{
			updateBC(timeSim, timeScaleUpdateBC);
		}
        for (int iRK = 0; iRK<3; iRK ++)
        {
			dudx = derivX(uu, xgrid, Br, nx, nvar);

			if (ifViscosity == 1)
			{
				d2udx2 = derivXX(uu, xgrid, Br, nx, nvar);
			}
			else
			{
				d2udx2 = (real*)malloc(sizeof(real));
			}

			calcRHS(dudt,uu,dudx,d2udx2,xgrid,Br,wavelength,
				A,dlnAdr,Q_AH,CW,CRflect,nx,nvar);

            rkExecute(iRK, cc1, dd1, uu, dudt, dudtRK, arrSize);

			filter(uu, nx, nvar);

			setFirstPoint(uu, nx, nvar);

			free(dudx);
			free(d2udx2);
        }

        timeSim+= dt;
        iStep++;
    }
    
    // after Tmax is reached or abnormal break, output one file 
    writeOutput(iOutput,timeSim,uu,arrSize,nx*nvar);

    time(&tRunTime);
    diffTime = difftime(tRunTime, tStart);

    printf("Simulation Ends: iStep = %ld, timeSim = %e, dt = %e\n",
        iStep, timeSim, dt);
    printf("Total simulation time = %12lds\n",(long)diffTime);
    
	system("pause");
    return 0;
}


int checkNan(real *arr)
{
    int isNan = isnan(arr[0]);
    int isNanAll = 0;

	//isNanAll = isNan;

	for (int i = 0; i <arrSize; i++)
	{
		if(isnan(arr[i])) isNanAll++;
	}

    
    // MPI_Allreduce(&isNan, &isNanAll, 1, MPI_INT, 
    //     MPI_MAX, MPI_COMM_WORLD);
    
    return isNanAll;
}



void writeLog(real timeSim, real dt, real diffTime, size_t iStep)
{
	long seconds = (long)diffTime;
	int hour, minute, second;

	hour = seconds / 3600;
	minute = (seconds % 3600) / 60;
	second = seconds % 60;

	//FILE *fpWrite = fopen("log", "w");
	FILE *fpWrite;
	errno_t err = fopen_s(&fpWrite, "log", "w");

	fprintf(fpWrite, "  Simulation time = %f\n", timeSim);
	fprintf(fpWrite, "  dt = %12.4e\n", dt);
	fprintf(fpWrite, "  Real time used (sec) = %12ld\n", seconds);
	fprintf(fpWrite, "  Real time used (hh:mm:ss) = (%03d:%02d:%02d)\n",
		hour, minute, second);
	fprintf(fpWrite, "  Iteration = %14ld\n", iStep);
	fclose(fpWrite);
}