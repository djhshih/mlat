/* gfResult - In-memory results from genoFind */

#ifndef _GFRESULT_H_
#define _GFRESULT_H_

#include <stddef.h>

struct gfAlignBlock {
  int size;
  int qStart;
  int tStart;
};

struct gfAlign {
  char *tName;
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
  /* externally read-only members */

  /* alignment result */
  struct gfAlign *aligns;
  /* number of alignments */
  size_t size;
  /* currently allocated space */
  size_t capacity;
};

struct gfResult *newGfResult();

struct gfResult *cloneGfResult(struct gfResult *other);


void freeGfResult(struct gfResult **pp);

#endif /* GFRESULT_H */
