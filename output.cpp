#include <stdio.h>
#include <stdlib.h>
#include "macros.h"
#include "ios_support.h"

void write1DArray(int nx, real *xgrid, const char *filename)
{
    FILE *fpWrite;

	real nxf = (real)nx;

    // fpWrite = fopen(filename,"wb");
	errno_t err = fopen_s(&fpWrite, filename, "wb");

    fwrite(&nxf, sizeof(real), 1, fpWrite);
    fwrite(xgrid, sizeof(real), nx, fpWrite);
    fclose(fpWrite);
}

int readGrid(real **ptrGrid, const char *filename)
{
	FILE *fpRead;

	real nxf;

	// fpRead = fopen(filename, "rb");
	errno_t err = fopen_s(&fpRead, filename, "rb");

	if (fpRead == NULL)
	{
		printf("Grid File: %s not found!\n", filename);
		return -1;
	}

	fread(&nxf, sizeof(real), 1, fpRead);

	int nx = (int)nxf;

	*ptrGrid = (real *)malloc(sizeof(real) * nx);
	fread(*ptrGrid, sizeof(real), nx, fpRead);

	fclose(fpRead);
	
	return nx;
}

void writeArray(int iOutput, real timeSim, real *arr, 
    size_t arrSize, int ifPrim)
{
    int bufsz,snprintfReturn;
    char *fileName;

    if(ifPrim==0)
    {
        bufsz = snprintf(NULL, 0, "out%04d.dat", iOutput);
        fileName = (char *)malloc(sizeof(char) * (bufsz+1));
        snprintfReturn = snprintf(fileName,bufsz + 1, "out%04d.dat",iOutput);
        printf("Output file out%04d.dat at time %12.4e\n", iOutput, timeSim);
    }
    else
    {
        bufsz = snprintf(NULL, 0, "prout%04d.dat", iOutput);
        fileName = (char *)malloc(sizeof(char) * (bufsz+1));
        int snprintfReturn = snprintf(fileName,bufsz + 1, "prout%04d.dat",iOutput);
        printf("Output file prout%04d.dat at time %12.4e\n", iOutput, timeSim);
    }

    FILE *fpWrite;
    // fpWrite = fopen(fileName,"wb");
	errno_t err = fopen_s(&fpWrite, fileName, "wb");
    fwrite(&timeSim, sizeof(real), 1, fpWrite);
    fwrite(arr, sizeof(real), arrSize, fpWrite);
    fclose(fpWrite);

	free(fileName);
}


void writeOutput(int iOutput, real timeSim, real *uu, 
    size_t arrSize, size_t arrPrSize)
{
    writeArray(iOutput,timeSim,uu,arrSize,0);
}

real readArray(const char *fileName, real *arr, size_t arrSize)
{
	FILE *fpRead;
	real timeSim;

	// fpRead = fopen(fileName, "rb");
	errno_t err = fopen_s(&fpRead, fileName, "rb");
	fread(&timeSim, sizeof(real), 1, fpRead);
	fread(arr, sizeof(real), arrSize, fpRead);
	fclose(fpRead);

	return timeSim;
}