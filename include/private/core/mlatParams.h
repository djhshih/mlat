/** Parameters for mlat **/

#ifndef _MLAT_PARAMS_H_
#define _MLAT_PARAMS_H_

#include "aliType.h"

struct mlatParams {
  /* Integer parameters */
  int tileSize;
  int stepSize;
  int minMatch;
  int minScore;
  int maxGap;
  int repMatch;
  /* Boolean flags */
  int oneOff;
  int noHead;
  int trimA;
  int trimHardA;
  int trimT;
  int fastMap;
  int fine;
  /* Other parameters */
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

#endif
