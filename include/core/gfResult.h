/* gfResult - In-memory results from genoFind */

#ifndef _GFRESULT_H_
#define _GFRESULT_H_

#include "common.h"
#include "dnaseq.h"
#include "fuzzyfind.h"
#include "hash.h"
#include "aliType.h"
#include "localmem.h"
#include "bits.h"
#include "genoFind.h"
#include "trans3.h"

struct gfAlignBlock {
  int size;
  int qStart;
  int tStart;
};

struct gfAlign {
  int qStart, qEnd;
  int tStart, tEnd;
  int qInsertBaseCount;
  int qInsertCount;
  int tInsertBaseCount;
  int tInsertCount;
  int matchCount;
  int mismatchCount;
  int repMatchCount;
  int nCount;
  char qStrand;
  char tStrand;
  int score;
  struct gfAlignBlock *blocks;
  size_t blockCount;
};

struct gfResult {
  char *tName;
  char *qName;

  /* externally read-only members */

  /* alignment result */
  struct gfAlign *aligns;
  /* number of alignments */
  size_t size;
  /* currently allocated space */
  size_t capacity;
};

struct gfResult *newGfResult();

void freeGfResult(struct gfResult **pp);

/* free with freeGfOutputResult */
struct gfOutput *gfOutputResult(int minGood, boolean qIsProt, boolean tIsProt);

/* Free gfOutputResult */
void freeGfOutputResult(struct gfOutput **pp);

/* Free gfOutputResult but return the contained gfResult */
struct gfResult *unpackGfOutputResult(struct gfOutput **pp);

#endif /* GFRESULT_H */
