#ifndef _IOS_SUPPORT_
#define _IOS_SUPPORT_

typedef int errno_t;

errno_t fopen_s(FILE **fp, const char *filename, const char *mode);

#endif