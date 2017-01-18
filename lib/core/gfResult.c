/* gfResult - In-memory results from genoFind */

#include "gfResult.h"

static void gfResultOut(char *qName, int qSize, int qOffset, struct ffAli *align,
		bioSeq *tSeq, struct hash *t3Hash, bioSeq *qSeq,
		boolean qIsRc, boolean tIsRc, enum ffStringency stringency, int score,
		struct gfOutput *out)
{
	// Basd on gfOut.c savePslx(.)
	// NB  this function needs to have the same signature as gfOutput::out
	// NB  minMatch is not used, because upstream has already done filtering
	//     minMatch is changed to score (i.e. saveAlignments)

	struct hash *maskHash = out->maskHash;
	int minIdentity = out->minGood * 10;
	boolean qIsProt = out->qIsProt;
	boolean tIsProt = out->tIsProt;

  struct ffAli *ff, *nextFf;
  struct ffAli *right = ffRightmost(align);

	struct gfResult *r = out->data;
	// clone the strings because they may be deallocated during r's lifetime
	r->tName = cloneString(tSeq->name);
	r->qName = cloneString(qSeq->name);

  DNA *needle = qSeq->dna;
  DNA *hay = tSeq->dna;

  int nStart = align->nStart - needle;
  int nEnd = right->nEnd - needle;
  int hStart, hEnd;
  int nInsertBaseCount = 0;
  int nInsertCount = 0;
  int hInsertBaseCount = 0;
  int hInsertCount = 0;
  int matchCount = 0;
  int mismatchCount = 0;
  int repMatch = 0;
  int countNs = 0;

  DNA *np, *hp, n, h;
  int blockSize;
  int i;
  struct trans3 *t3List = NULL;
  Bits *maskBits = NULL;

  if (maskHash != NULL)
    maskBits = hashMustFindVal(maskHash, tSeq->name);

  if (t3Hash != NULL)
    t3List = hashMustFindVal(t3Hash, tSeq->name);

  hStart = trans3GenoPos(align->hStart, tSeq, t3List, FALSE) + qOffset;
  hEnd = trans3GenoPos(right->hEnd, tSeq, t3List, TRUE) + qOffset;

  /* Count up matches, mismatches, inserts, etc. */
  for (ff = align; ff != NULL; ff = nextFf) {
    nextFf = ff->right;
    blockSize = ff->nEnd - ff->nStart;
    np = ff->nStart;
    hp = ff->hStart;
    for (i = 0; i < blockSize; ++i) {
      n = np[i];
      h = hp[i];
      if (n == 'n' || h == 'n')
        ++countNs;
      else {
        if (n == h) {
          if (maskBits != NULL) {
            int seqOff = hp + i - hay;
            if (bitReadOne(maskBits, seqOff))
              ++repMatch;
            else
              ++matchCount;
          } else
            ++matchCount;
        } else
          ++mismatchCount;
      }
    }
    if (nextFf != NULL) {
      int nhStart =
          trans3GenoPos(nextFf->hStart, tSeq, t3List, FALSE) + qOffset;
      int ohEnd = trans3GenoPos(ff->hEnd, tSeq, t3List, TRUE) + qOffset;
      int hGap = nhStart - ohEnd;
      int nGap = nextFf->nStart - ff->nEnd;

      if (nGap != 0) {
        ++nInsertCount;
        nInsertBaseCount += nGap;
      }
      if (hGap != 0) {
        ++hInsertCount;
        hInsertBaseCount += hGap;
      }
    }
  }

	int gaps = nInsertCount + (stringency == ffCdna ? 0 : hInsertCount);
	if ((matchCount + repMatch + mismatchCount) > 0) {
		// `id` is scaled to be 10 * percentage
		// so `minIdentity` needs to be 10 * percentage (done above)
		int id = roundingScale(1000, matchCount + repMatch - 2 * gaps,
													 matchCount + repMatch + mismatchCount);
		if (id >= minIdentity) {
			if (qIsRc) {
				int temp;
				int oSize = qSeq->size;
				temp = nStart;
				nStart = oSize - nEnd;
				nEnd = oSize - temp;
			}
			if (tIsRc) {
				int temp;
				temp = hStart;
				hStart = qSize - hEnd;
				hEnd = qSize - temp;
			}

			// get pointer for new hit
			struct gfAlign *hit = gfResultNewHit(r);
			hit->qStart = nStart;
			hit->qEnd = nEnd;
			hit->tStart = hStart;
			hit->tEnd = hEnd;
			hit->qInsertBaseCount = nInsertBaseCount;
			hit->qInsertCount = nInsertCount;
			hit->tInsertBaseCount = hInsertBaseCount;
			hit->tInsertCount = hInsertCount;
			hit->matchCount = matchCount;
			hit->mismatchCount = mismatchCount;
			hit->repMatchCount = repMatch;
			hit->nCount = countNs;
			hit->qStrand = (qIsRc ? '-' : '+');
			hit->tStrand = (tIsRc ? '-' : '+');
			hit->score = score;

			size_t blockCount = ffAliCount(align);

			/* Ideally, this could be enclosed within a function.
			 * However, trans3GenoPos with its many passed parameters makes this
			 * task unwieldy without closures or functors ... */
			struct gfAlignBlock* blocks = AllocN(struct gfAlignBlock, blockCount);
			for (ff = align, i = 0; ff != NULL; ff = ff->right) {
				struct gfAlignBlock *block = &blocks[i++];
				block->size = ff->nEnd - ff->nStart;
				block->qStart = ff->nStart - needle;
				block->tStart = trans3GenoPos(ff->hStart, tSeq, t3List, FALSE) + qOffset;
			}
			hit->blocks = blocks;
			
		} /* if (id >= minIdentity) */
	} /* if ((matchCount + repMatch + mismatchCount) > 0) */
}

struct gfResult *newGfResult()
{
	struct gfResult *r;
	AllocVar(r);
	ZeroVar(r);

	// set an initial capacity
	r->capacity = 10;
	r->aligns = AllocN(struct gfAlign, r->capacity);

	return r;
}

struct gfAlign *gfResultNewHit(struct gfResult* r) {
	++(r->size);
	if (r->size > r->capacity) {
		size_t oldCapacity = r->capacity;
		r->size *= 2;
		ExpandArray(r->aligns, oldCapacity, r->capacity);
	}
	return &r->aligns[r->size - 1];
}

void freeGfResult(struct gfResult **pp)
{
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

// free with freeGfOutputResult
struct gfOutput *gfOutputResult(int minGood, boolean qIsProt, boolean tIsProt)
{
	struct gfOutput *out = gfOutputInit(minGood, qIsProt, tIsProt);
	out->out = gfResultOut;
	out->data = newGfResult();

	return out;
}

void freeGfOutputResult(struct gfOutput **pp)
{
	struct gfOutput *p = *pp;
	if (p != NULL) {
		freeGfResult((struct gfResult **)&p->data);
		freez(pp);
	}
}

