#include "mlat.h"

/* TODO const-correctness on:
 * const char* querySeq -> const char* bioSeq.dna -> worms
 * const struct mlatParams* p
 */

/* Open and read a database sequence, mask, and index it */
/* Free gfDb* with freeGfDb */
/* TODO support other tType and qType with calls to bigMlat */
struct gfDb *newGfDb(char *dbFile, struct mlatParams *p) {
  struct gfDb *db;
  AllocVar(db);
  db->maskHash = NULL;

  const boolean showStatus = FALSE;

  bool bothSimpleNuc =
      (p->tType == gftDna && (p->qType == gftDna || p->qType == gftRna));
  bool bothSimpleProt = (p->tType == gftProt && p->qType == gftProt);
  // only valid for simple query and target types
  if (!(bothSimpleNuc || bothSimpleProt))
    return NULL;

  struct dnaSeq *seq;
  char **dbFiles;
  int dbCount;

  gfClientFileArray(dbFile, &dbFiles, &dbCount);

  db->seqList = gfClientSeqList(dbCount, dbFiles, p->tType == gftProt,
                                p->tType == gftDnaX, p->repeats,
                                p->minRepDivergence, showStatus);
  db->seqCount = slCount(db->seqList);

  // save masking info
  if (p->repeats != NULL) {
    db->maskHash = newHash(0);
    for (seq = db->seqList; seq != NULL; seq = seq->next) {
      hashAdd(db->maskHash, seq->name, maskFromUpperCaseSeq(seq));
    }
  }

  if (p->mask == NULL) {
    gfClientUnmask(db->seqList);
  }

  // create indexer
  db->gf =
      gfIndexSeq(db->seqList, p->minMatch, p->maxGap, p->tileSize, p->repMatch,
                 p->ooc, p->tType == gftProt, p->oneOff, FALSE, p->stepSize);

  // turn off masking, because sequences should always be unmaksed for the
  // extension phase
  if (p->mask != NULL) {
    gfClientUnmask(db->seqList);
  }

  freeArrays((void **)dbFiles, dbCount);

  return db;
}

void freeGfDb(struct gfDb **pDb) {
  struct gfDb *db = *pDb;
  if (db != NULL) {
    genoFindFree(&db->gf);
    freeHash(&db->maskHash);
    freeDnaSeqList(&db->seqList);
    freez(pDb);
  }
}

/* Search for a query sequence in the database index on one DNA strand */
/* FIXME reconsile with searchOneStrand */
static void searchDnaStrand(struct gfDb *db, bioSeq *seq, boolean isRc,
                     Bits *qMaskBits, struct mlatParams *p,
                     struct gfOutput *gvo) {
  if (p->fastMap && (seq->size > MAXSINGLEPIECESIZE)) {
    errAbort(
        "Maximum single piece size (%d) exceeded by query %s of size (%d). "
        "Larger pieces will have to be split up until no larger than this "
        "limit when the -fastMap option is used.",
        MAXSINGLEPIECESIZE, seq->name, seq->size);
  }

  gfLongDnaInMem(seq, db->gf, isRc, p->minScore, qMaskBits, gvo, p->fastMap,
                 p->fine);
}

/* Search a query sequence against a target database index */
/* Free gfOutput* with freeGfOutputResult */
struct gfOutput *searchSeq(struct gfDb *db, char *querySeq,
                           struct mlatParams *p) {
  // create temporary shallow copy
  bioSeq seq0;
  ZeroVar(&seq0);
  seq0.dna = querySeq;
  seq0.size = strlen(querySeq);

  // apply mask to query sequence
  boolean lcMask = (p->qMask != NULL && sameWord(p->qMask, "lower"));
  Bits *qMaskBits =
      maskQuerySeq(&seq0, p->tType == gftProt, p->qMask != NULL, lcMask);

  // trim query sequence, creating a deep copy
  bioSeq seq;
  ZeroVar(&seq);
  trimSeq(&seq0, &seq, p);

  // convert u to t for RNA sequences
  if (p->qType == gftRna || p->qType == gftRnaX)
    memSwapChar(seq.dna, seq.size, 'u', 't');

  struct gfOutput *gvo =
      gfOutputResult(p->minIdentity, p->qType == gftProt, p->tType == gftProt);

  if (p->tType == gftProt) {
    searchOneProt(&seq, db->gf, p->minScore, gvo);
  } else {
    gvo->maskHash = db->maskHash;
    // search forward strand
    searchDnaStrand(db, &seq, FALSE, qMaskBits, p, gvo);

    // search reverse strand
    reverseComplement(seq.dna, seq.size);
    searchDnaStrand(db, &seq, TRUE, qMaskBits, p, gvo);
    reverseComplement(seq.dna, seq.size);
  }

  bitFree(&qMaskBits);

  return gvo;
}
