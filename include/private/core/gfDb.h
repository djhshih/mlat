#include "mlatParams.h"

#ifndef _GFDB_H_
#define _GFDB_H_

/* Private definition */
struct gfDb;

/* Open and read a database sequence, mask, and index it */
/* Free gfDb* with freeGfDb */
struct gfDb *newGfDb(char *dbFile, struct mlatParams *p);

/* Free gfDb using reference to pointer */
void freeGfDb(struct gfDb **pDb);

/* Search a query sequence against a target database index */
/* Free gfResult* with freeGfResult */
struct gfResult *searchSeq(struct gfDb *db, char *querySeq,
                           struct mlatParams *p);

#endif /* _GFDB_H_ */
