//#include "mlatPriv.h"
#include "gfResult.h"

#include "aliType.h"

#ifndef _MLAT_ALT_H_
#define _MLAT_ALT_H_

#define boolean int

struct mlatParams {
  int tileSize;
  int stepSize;
  int minMatch;
  int minScore;
  int maxGap;
  int repMatch;
  boolean oneOff;
  boolean noHead;
  boolean trimA;
  boolean trimHardA;
  boolean trimT;
  boolean fastMap;
  boolean fine;
  char *makeOoc;
  char *ooc;
  enum gfType qType;
  enum gfType tType;
  char *mask;
  char *repeats;
  char *qMask;
  double minRepDivergence;
  double minIdentity;
  char *outputFormat;
};

/* Initialize default mlat parameters */
/* Free mlatParams* with freeMlatParams */
struct mlatParams *newMlatParams();

/* Free mlatParams using reference to pointer */
void freeMlatParams(struct mlatParams **pp);

/* Container for sequence database */
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

#endif /* _MLAT_ALT_H_ */
