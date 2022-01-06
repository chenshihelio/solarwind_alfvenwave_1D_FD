#ifndef _MACROS_
#define _MACROS_

typedef double real;

#define IDX(i,iv,nv) (i) * (nv) + (iv)
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

#define NEigen 5
#define NVAR 6

#define NPOINTFD 4 // finite-difference scheme order

#endif