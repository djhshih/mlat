#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <stddef.h>

/* Let's pretend C has a boolean type. */
#define TRUE 1
#define FALSE 0
#define boolean int
#ifndef __cplusplus
#ifndef bool
#define bool char
#endif
#endif

#endif
