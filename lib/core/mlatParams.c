#include "mlatParams.h"
#include "aliType.h"
#include "common.h"

struct mlatParams *newMlatParams() {
  struct mlatParams *p;
  AllocVar(p);

  p->tileSize = 11;
  p->stepSize = 5;
  p->minMatch = 2;
  p->minScore = 30;
  p->maxGap = 2;
  p->repMatch = 1024 * p->tileSize / p->stepSize;
  p->oneOff = FALSE;
  p->noHead = FALSE;
  p->trimA = FALSE;
  p->trimHardA = FALSE;
  p->trimT = FALSE;
  p->fastMap = FALSE;
  p->makeOoc = NULL;
  p->ooc = NULL;
  p->qType = gftDna;
  p->tType = gftDna;
  p->mask = NULL;
  p->repeats = NULL;
  p->qMask = NULL;
  p->minRepDivergence = 15;
  p->minIdentity = 90;
  p->outputFormat = NULL;

  return p;
}

void freeMlatParams(struct mlatParams **pp) { freez(pp); }
