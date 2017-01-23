#ifndef _MLAT_H_
#define _MLAT_H_

#include "common.h"
#include "memalloc.h"
#include "linefile.h"
#include "bits.h"
#include "hash.h"
#include "dnautil.h"
#include "dnaseq.h"
#include "fa.h"
#include "twoBit.h"
#include "psl.h"
#include "sig.h"
#include "options.h"
#include "obscure.h"
#include "genoFind.h"
#include "trans3.h"
#include "gfClientLib.h"

enum constants {
  qWarnSize = 5000000, /* Warn if more than this many bases in one query. */
};

/* Parameters of mlat passed down the function call chain
 *
 * Do not dynamically allocate char* members (e.g. using malloc, or Kent's
 * functions cloneString, AllocVar, AllocN, needMem, etc.).
 *
 * Simply assign char* members to point to statically allocated or dynamically
 * allocated by the function that owns the memory.
 *
 * With this constraint, mlatParams* p can be deallocated by freez(&p).
 */
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

struct mlatParams *newMlatParams();

/* Search for seq in index, align it, and write results to psl. */
void searchOneStrand(struct dnaSeq *seq, struct genoFind *gf, FILE *psl,
                     boolean isRc, struct hash *maskHash, Bits *qMaskBits,
                     struct mlatParams *p, struct gfOutput *gvo);

/* Search for protein seq in index and write results to psl. */
void searchOneProt(aaSeq *seq, struct genoFind *gf, FILE *f, int minScore,
                   struct gfOutput *gvo);

/* Search for seq on either strand in index. */
void searchOne(bioSeq *seq, struct genoFind *gf, FILE *f, struct hash *maskHash,
               Bits *qMaskBits, struct mlatParams *p, struct gfOutput *gvo);

/* Copy seq to trimmed (shallow copy) and optionally trim
 * off polyA tail or polyT head. */
void trimSeq(struct dnaSeq *seq, struct dnaSeq *trimmed, struct mlatParams *p);

/* Massage query sequence a bit, converting it to correct
 * case (upper for protein/lower for DNA) and optionally
 * returning upper/lower case info , and trimming poly A. */
Bits *maskQuerySeq(struct dnaSeq *seq, boolean tIsProt, boolean maskQuery,
                   boolean lcMask);

/* Search a single sequence against a single genoFind index. */
void searchOneMaskTrim(struct dnaSeq *seq, struct genoFind *gf, FILE *outFile,
                       struct hash *maskHash, struct mlatParams *p,
                       struct gfOutput *gvo, long long *retTotalSize,
                       int *retCount);

/* Search all sequences in all files against single genoFind index. */
void searchOneIndex(int fileCount, char *files[], struct genoFind *gf,
                    char *outName, struct hash *maskHash, FILE *outFile,
                    boolean showStatus, struct mlatParams *p,
                    struct gfOutput *gvo);

/* Convert sequence list to a trans3 list and lists for each of three frames. */
struct trans3 *seqListToTrans3List(struct dnaSeq *seqList, aaSeq *transLists[3],
                                   struct hash **retHash);

/* Look for qSeq in indices for three frames.  Then do rest of alignment. */
void tripleSearch(aaSeq *qSeq, struct genoFind *gfs[3], struct hash *t3Hash,
                  boolean dbIsRc, FILE *f, int minScore, struct gfOutput *gvo);

/* Translate qSeq three ways and look for each in three frames of index. */
void transTripleSearch(struct dnaSeq *qSeq, struct genoFind *gfs[3],
                       struct hash *t3Hash, boolean dbIsRc, boolean qIsDna,
                       FILE *f, int minScore, struct gfOutput *gvo);

/* Run query against translated DNA database (3 frames on each strand). */
void bigMlat(struct dnaSeq *untransList, int queryCount, char *queryFiles[],
             char *outFile, boolean transQuery, boolean qIsDna, FILE *out,
             boolean showStatus, struct mlatParams *p, struct gfOutput *gvo);

/* mlat - Minimal BLAT fast sequence search command line tool. */
void mlat(char *dbFile, char *queryFile, char *outName, struct mlatParams *p);

#endif /* _MLAT_H_ */
