#include "mlat.h"
#include "gfResult.h"

#ifndef _MLAT_ALT_H_
#define _MLAT_ALT_H_

struct gfDb {
	struct genoFind *gf;
	bioSeq *seqList;
	int seqCount;
	struct hash *maskHash;
};

/* Open and read a database sequence, mask, and index it */
/* Free gfDb* with freeGfDb */
struct gfDb *newGfDb(char *dbFile, struct mlatParams *p);
void freeGfDb(struct gfDb **pDb);

/* Search for seq in index and align it */
void searchDnaStrand(struct gfDb *db, bioSeq *seq, boolean isRc, Bits *qMaskBits, struct mlatParams *p, struct gfOutput *gvo);

/* Search a query sequence against a target database index */
/* Free gfOutput* with freeGfOutputResult */
struct gfOutput *searchSeq(struct gfDb *db, char *querySeq, struct mlatParams *p);

#endif /* _MLAT_ALT_H_ */
