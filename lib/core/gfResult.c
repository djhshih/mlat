/* gfResult - In-memory results from genoFind */

#include "gfResult.h"
#include "common.h"

struct gfResult *newGfResult() {
  struct gfResult *r;
  AllocVar(r);
  ZeroVar(r);

  // set an initial capacity
  r->capacity = 10;
  r->aligns = AllocN(struct gfAlign, r->capacity);

  return r;
}

void freeGfResult(struct gfResult **pp) {
  struct gfResult *p = *pp;
  if (p != NULL) {
    freeMem(p->tName);
    freeMem(p->qName);
    size_t i;
    for (i = 0; i < p->size; ++i) {
      freeMem(p->aligns[i].blocks);
    }
    freeMem(p->aligns);
    freez(pp);
  }
}

