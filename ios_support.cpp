#include <stdlib.h>
#include <stdio.h>
#include "ios_support.h"

// In order to make the code portable between Windows and IOS
errno_t fopen_s(FILE **fp, const char *filename, const char *mode)
{

    errno_t stat = 0;

    *fp = fopen(filename, mode);

    return stat;
}