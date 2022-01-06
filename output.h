#ifndef _OUTPUT_
#define _OUTPUT_

void write1DArray(int nx, real *xgrid, const char *filename);

int readGrid(real **ptrGrid, const char *filename);

void writeArray(int iOutput, real timeSim, real *arr, 
    size_t arrSize, int ifPrim);

void writeOutput(int iOutput, real timeSim, real *uu, 
    size_t arrSize, size_t arrPrSize);

real readArray(const char *fileName, real *arr, size_t arrSize);
#endif