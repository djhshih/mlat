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

static struct gfAlign *cloneGfAlign(struct gfAlign *other) {
  struct gfAlign *x;
  AllocVar(x);
  *x = *other;

  x->blocks = AllocN(struct gfAlignBlock, x->blockCount);
  for (size_t i = 0; i < x->blockCount; ++i) {
    x->blocks[i] = other->blocks[i];
  }

  return x;
}

struct gfResult *cloneGfResult(struct gfResult *other) {
  struct gfResult *x = newGfResult();
  x->tName = cloneString(other->tName);
  x->qName = cloneString(other->qName);

  x->size = other->size;
  x->capacity = other->capacity;

  x->aligns = AllocN(struct gfAlign, x->capacity);
  for (size_t i = 0; i < x->size; ++i) {
    const struct gfAlign *align = &other->aligns[i];
    x->aligns[i] = *align;

    size_t blockCount = x->aligns[i].blockCount;
    x->aligns[i].blocks = AllocN(struct gfAlignBlock, blockCount);
    for (size_t j = 0; j < blockCount; ++j) {
      x->aligns[i].blocks[j] = align->blocks[j];
    }
  }
  return x;
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
