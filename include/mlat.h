#ifdef __cplusplus
extern "C"{
#endif 

#include "mlat/gfResult.h"
#include "mlat/aliType.h"
#include "mlat/mlatParams.h"
#include "mlat/gfDb.h"

/* minimal BLAT fast sequence alignment tool */
void mlat(char *dbFile, char *queryFile, char *outName, struct mlatParams *p);

#ifdef __cplusplus
}
#endif 
