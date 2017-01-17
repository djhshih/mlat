/* Copyright 2001-2004 Jim Kent.  All rights reserved. */
#include "mlat.h"

struct mlatParams *makeMlatParams()
{
  struct mlatParams *p;
  AllocVar(p);

  p->tileSize = 11;
  p->stepSize = 0;  /* Default (same as tileSize) */
  p->minMatch = 2;
  p->minScore = 30;
  p->maxGap = 2;
  p->repMatch = 1024 * 4;
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

void searchOneStrand(struct dnaSeq *seq, struct genoFind *gf, FILE *psl,
                     boolean isRc, struct hash *maskHash, Bits *qMaskBits,
                     struct mlatParams *p, struct gfOutput *gvo)
/* Search for seq in index, align it, and write results to psl. */
{
  if (p->fastMap && (seq->size > MAXSINGLEPIECESIZE))
    errAbort(
        "Maximum single piece size (%d) exceeded by query %s of size (%d). "
        "Larger pieces will have to be split up until no larger than this "
        "limit "
        "when the -fastMap option is used.",
        MAXSINGLEPIECESIZE, seq->name, seq->size);
  gfLongDnaInMem(seq, gf, isRc, p->minScore, qMaskBits, gvo, p->fastMap,
                 optionExists("fine"));
}

void searchOneProt(aaSeq *seq, struct genoFind *gf, FILE *f, int minScore, struct gfOutput *gvo)
/* Search for protein seq in index and write results to psl. */
{
  int hitCount;
  struct lm *lm = lmInit(0);
  struct gfClump *clumpList = gfFindClumps(gf, seq, lm, &hitCount);
  gfAlignAaClumps(gf, clumpList, seq, FALSE, minScore, gvo);
  gfClumpFreeList(&clumpList);
  lmCleanup(&lm);
}

void searchOne(bioSeq *seq, struct genoFind *gf, FILE *f,
               struct hash *maskHash, Bits *qMaskBits, struct mlatParams* p,
               struct gfOutput *gvo)
/* Search for seq on either strand in index. */
{
  if (p->tType == gftProt) {
    searchOneProt(seq, gf, f, p->minScore, gvo);
  } else {
    gvo->maskHash = maskHash;
    searchOneStrand(seq, gf, f, FALSE, maskHash, qMaskBits, p, gvo);
    reverseComplement(seq->dna, seq->size);
    searchOneStrand(seq, gf, f, TRUE, maskHash, qMaskBits, p, gvo);
    reverseComplement(seq->dna, seq->size);
  }
  gfOutputQuery(gvo, f);
}

void trimSeq(struct dnaSeq *seq, struct dnaSeq *trimmed, struct mlatParams* p)
/* Copy seq to trimmed (shallow copy) and optionally trim
 * off polyA tail or polyT head. */
{
  DNA *dna = seq->dna;
  int size = seq->size;
  *trimmed = *seq;
  if (p->trimT)
    maskHeadPolyT(dna, size);
  if (p->trimA || p->trimHardA) {
    int trimSize = maskTailPolyA(dna, size);
    if (p->trimHardA) {
      trimmed->size -= trimSize;
      dna[size - trimSize] = 0;
    }
  }
}

Bits *maskQuerySeq(struct dnaSeq *seq, boolean tIsProt, boolean maskQuery,
                   boolean lcMask)
/* Massage query sequence a bit, converting it to correct
 * case (upper for protein/lower for DNA) and optionally
 * returning upper/lower case info , and trimming poly A. */
{
  Bits *qMaskBits = NULL;
  verbose(2, "%s\n", seq->name);
  if (tIsProt)
    faToProtein(seq->dna, seq->size);
  else {
    if (maskQuery) {
      if (lcMask)
        toggleCase(seq->dna, seq->size);
      qMaskBits = maskFromUpperCaseSeq(seq);
    }
    faToDna(seq->dna, seq->size);
  }
  if (seq->size > qWarnSize) {
    warn("Query sequence %s has size %d, it might take a while.", seq->name,
         seq->size);
  }
  return qMaskBits;
}

void searchOneMaskTrim(struct dnaSeq *seq, struct genoFind *gf, FILE *outFile,
                       struct hash *maskHash, struct mlatParams* p,
                       struct gfOutput *gvo,
                       long long *retTotalSize, int *retCount)
/* Search a single sequence against a single genoFind index. */
{
	boolean tIsProt = (p->tType == gftProt);
  boolean maskQuery = (p->qMask != NULL);
  boolean lcMask = (p->qMask != NULL && sameWord(p->qMask, "lower"));
  Bits *qMaskBits = maskQuerySeq(seq, tIsProt, maskQuery, lcMask);
  struct dnaSeq trimmedSeq;
  ZeroVar(&trimmedSeq);
  trimSeq(seq, &trimmedSeq, p);
  if (p->qType == gftRna || p->qType == gftRnaX)
    memSwapChar(trimmedSeq.dna, trimmedSeq.size, 'u', 't');
  searchOne(&trimmedSeq, gf, outFile, maskHash, qMaskBits, p, gvo);
  *retTotalSize += seq->size;
  *retCount += 1;
  bitFree(&qMaskBits);
}

void searchOneIndex(int fileCount, char *files[], struct genoFind *gf,
                    char *outName, struct hash *maskHash, FILE *outFile, 
                    boolean showStatus, struct mlatParams* p, struct gfOutput *gvo)
/* Search all sequences in all files against single genoFind index. */
{
  int i;
  char *fileName;
  int count = 0;
  long long totalSize = 0;

  gfOutputHead(gvo, outFile);
  for (i = 0; i < fileCount; ++i) {
    fileName = files[i];
    if (twoBitIsSpec(fileName)) {
      struct twoBitSpec *tbs = twoBitSpecNew(fileName);
      struct twoBitFile *tbf = twoBitOpen(tbs->fileName);
      if (p->tType == gftProt)
        errAbort("%s is a two bit file, which doesn't work for proteins.",
                 fileName);
      if (tbs->seqs != NULL) {
        struct twoBitSeqSpec *ss = NULL;
        for (ss = tbs->seqs; ss != NULL; ss = ss->next) {
          struct dnaSeq *seq =
              twoBitReadSeqFrag(tbf, ss->name, ss->start, ss->end);
          searchOneMaskTrim(seq, gf, outFile, maskHash, p, gvo,
                            &totalSize, &count);
          dnaSeqFree(&seq);
        }
      } else {
        struct twoBitIndex *index = NULL;
        for (index = tbf->indexList; index != NULL; index = index->next) {
          struct dnaSeq *seq = twoBitReadSeqFrag(tbf, index->name, 0, 0);
          searchOneMaskTrim(seq, gf, outFile, maskHash, p, gvo,
                            &totalSize, &count);
          dnaSeqFree(&seq);
        }
      }
      twoBitClose(&tbf);
    } else {
      static struct dnaSeq seq;
      struct lineFile *lf = lineFileOpen(fileName, TRUE);
      while (faMixedSpeedReadNext(lf, &seq.dna, &seq.size, &seq.name)) {
        searchOneMaskTrim(&seq, gf, outFile, maskHash, p, gvo,
                          &totalSize, &count);
      }
      lineFileClose(&lf);
    }
  }
  carefulClose(&outFile);
  if (showStatus)
    printf("Searched %lld bases in %d sequences\n", totalSize, count);
}

struct trans3 *seqListToTrans3List(struct dnaSeq *seqList, aaSeq *transLists[3],
                                   struct hash **retHash)
/* Convert sequence list to a trans3 list and lists for each of three frames. */
{
  int frame;
  struct dnaSeq *seq;
  struct trans3 *t3List = NULL, *t3;
  struct hash *hash = newHash(0);

  for (seq = seqList; seq != NULL; seq = seq->next) {
    t3 = trans3New(seq);
    hashAddUnique(hash, t3->name, t3);
    slAddHead(&t3List, t3);
    for (frame = 0; frame < 3; ++frame) {
      slAddHead(&transLists[frame], t3->trans[frame]);
    }
  }
  slReverse(&t3List);
  for (frame = 0; frame < 3; ++frame) {
    slReverse(&transLists[frame]);
  }
  *retHash = hash;
  return t3List;
}

void tripleSearch(aaSeq *qSeq, struct genoFind *gfs[3], struct hash *t3Hash,
                  boolean dbIsRc, FILE *f, int minScore, struct gfOutput *gvo)
/* Look for qSeq in indices for three frames.  Then do rest of alignment. */
{
  gvo->reportTargetStrand = TRUE;
  gfFindAlignAaTrans(gfs, qSeq, t3Hash, dbIsRc, minScore, gvo);
}

void transTripleSearch(struct dnaSeq *qSeq, struct genoFind *gfs[3],
                       struct hash *t3Hash, boolean dbIsRc, boolean qIsDna,
                       FILE *f, int minScore, struct gfOutput *gvo)
/* Translate qSeq three ways and look for each in three frames of index. */
{
  int qIsRc;
  gvo->reportTargetStrand = TRUE;
  for (qIsRc = 0; qIsRc <= qIsDna; qIsRc += 1) {
    gfLongTransTransInMem(qSeq, gfs, t3Hash, qIsRc, dbIsRc, !qIsDna, minScore,
                          gvo);
    if (qIsDna)
      reverseComplement(qSeq->dna, qSeq->size);
  }
}

void bigMlat(struct dnaSeq *untransList, int queryCount, char *queryFiles[],
             char *outFile, boolean transQuery, boolean qIsDna, FILE *out,
             boolean showStatus, struct mlatParams *p, struct gfOutput *gvo)
/* Run query against translated DNA database (3 frames on each strand). */
{
  int frame, i;
  struct dnaSeq *seq, trimmedSeq;
  struct genoFind *gfs[3];
  aaSeq *dbSeqLists[3];
  struct trans3 *t3List = NULL;
  int isRc;
  struct lineFile *lf = NULL;
  struct hash *t3Hash = NULL;
  boolean forceUpper = FALSE;
  boolean forceLower = FALSE;
  boolean toggle = FALSE;
  boolean maskUpper = FALSE;

  ZeroVar(&trimmedSeq);
  if (showStatus)
    printf("Blatx %d sequences in database, %d files in query\n",
           slCount(untransList), queryCount);

  /* Figure out how to manage query case.  Proteins want to be in
   * upper case, generally, nucleotides in lower case.  But there
   * may be repeatMasking based on case as well. */
  if (transQuery) {
    if (p->qMask == NULL)
      forceLower = TRUE;
    else {
      maskUpper = TRUE;
      toggle = !sameString(p->qMask, "upper");
    }
  } else {
    forceUpper = TRUE;
  }

  if (gvo->fileHead != NULL)
    gvo->fileHead(gvo, out);

  for (isRc = FALSE; isRc <= 1; ++isRc) {
    /* Initialize local pointer arrays to NULL to prevent surprises. */
    for (frame = 0; frame < 3; ++frame) {
      gfs[frame] = NULL;
      dbSeqLists[frame] = NULL;
    }

    t3List = seqListToTrans3List(untransList, dbSeqLists, &t3Hash);
    for (frame = 0; frame < 3; ++frame) {
      gfs[frame] = gfIndexSeq(dbSeqLists[frame], p->minMatch, p->maxGap, p->tileSize,
                              p->repMatch, p->ooc, TRUE, p->oneOff, FALSE, p->stepSize);
    }

    for (i = 0; i < queryCount; ++i) {
      aaSeq qSeq;

      lf = lineFileOpen(queryFiles[i], TRUE);
      while (faMixedSpeedReadNext(lf, &qSeq.dna, &qSeq.size, &qSeq.name)) {
        /* Put it into right case and optionally mask on case. */
        if (forceLower)
          toLowerN(qSeq.dna, qSeq.size);
        else if (forceUpper)
          toUpperN(qSeq.dna, qSeq.size);
        else if (maskUpper) {
          if (toggle)
            toggleCase(qSeq.dna, qSeq.size);
          upperToN(qSeq.dna, qSeq.size);
        }
        if (qSeq.size > qWarnSize) {
          warn("Query sequence %s has size %d, it might take a while.",
               qSeq.name, qSeq.size);
        }
        trimSeq(&qSeq, &trimmedSeq, p);
        if (transQuery)
          transTripleSearch(&trimmedSeq, gfs, t3Hash, isRc, qIsDna, out, p->minScore, gvo);
        else
          tripleSearch(&trimmedSeq, gfs, t3Hash, isRc, out, p->minScore, gvo);
        gfOutputQuery(gvo, out);
      }
      lineFileClose(&lf);
    }

    /* Clean up time. */
    trans3FreeList(&t3List);
    freeHash(&t3Hash);
    for (frame = 0; frame < 3; ++frame) {
      genoFindFree(&gfs[frame]);
    }

    for (seq = untransList; seq != NULL; seq = seq->next) {
      reverseComplement(seq->dna, seq->size);
    }
  }
  carefulClose(&out);
}

void mlat(char *dbFile, char *queryFile, char *outName, struct mlatParams *p)
/* mlat - Minimal BLAT fast sequence search command line tool. */
{
  char **dbFiles, **queryFiles;
  int dbCount, queryCount;
  struct dnaSeq *dbSeqList, *seq;
  struct genoFind *gf;
  boolean tIsProt = (p->tType == gftProt);
  boolean qIsProt = (p->qType == gftProt);
  boolean bothSimpleNuc =
      (p->tType == gftDna && (p->qType == gftDna || p->qType == gftRna));
  boolean bothSimpleProt = (tIsProt && qIsProt);
  FILE *f = mustOpen(outName, "w");
  boolean showStatus = (f != stdout);
  int databaseSeqCount = 0;          /* Number of sequences in database. */
  unsigned long databaseLetters = 0; /* Number of bases in database. */
	/* Stuff to support various output formats. */
	struct gfOutput *gvo; /* Overall output controller */

  gfClientFileArray(dbFile, &dbFiles, &dbCount);
  if (p->makeOoc != NULL) {
    gfMakeOoc(p->makeOoc, dbFiles, dbCount, p->tileSize, p->repMatch, p->tType);
    if (showStatus)
      printf("Done making %s\n", p->makeOoc);
    exit(0);
  }
  gfClientFileArray(queryFile, &queryFiles, &queryCount);
  dbSeqList = gfClientSeqList(dbCount, dbFiles, tIsProt, p->tType == gftDnaX,
                              p->repeats, p->minRepDivergence, showStatus);
  databaseSeqCount = slCount(dbSeqList);
  for (seq = dbSeqList; seq != NULL; seq = seq->next)
    databaseLetters += seq->size;

  gvo = gfOutputAny(p->outputFormat, p->minIdentity * 10, qIsProt, tIsProt, p->noHead,
                    dbFile, databaseSeqCount, databaseLetters,
                    p->minIdentity, f);

  if (bothSimpleNuc || bothSimpleProt) {
    struct hash *maskHash = NULL;

    /* Save away masking info for output. */
    if (p->repeats != NULL) {
      maskHash = newHash(0);
      for (seq = dbSeqList; seq != NULL; seq = seq->next) {
        Bits *maskedBits = maskFromUpperCaseSeq(seq);
        hashAdd(maskHash, seq->name, maskedBits);
      }
    }

    /* Handle masking and indexing.  If masking is off, we want the indexer
     * to see unmasked sequence, otherwise we want it to see masked.  However
     * after indexing we always want it unmasked, because things are always
     * unmasked for the extension phase. */
    if (p->mask == NULL && !bothSimpleProt)
      gfClientUnmask(dbSeqList);
    gf = gfIndexSeq(dbSeqList, p->minMatch, p->maxGap, p->tileSize, p->repMatch, p->ooc,
                    tIsProt, p->oneOff, FALSE, p->stepSize);
    if (p->mask != NULL)
      gfClientUnmask(dbSeqList);

    searchOneIndex(queryCount, queryFiles, gf, outName, 
                   maskHash, f, showStatus, p, gvo);
    freeHash(&maskHash);
  } else if (p->tType == gftDnaX && p->qType == gftProt) {
    bigMlat(dbSeqList, queryCount, queryFiles, outName, FALSE, TRUE, f,
            showStatus, p, gvo);
  } else if (p->tType == gftDnaX && (p->qType == gftDnaX || p->qType == gftRnaX)) {
    bigMlat(dbSeqList, queryCount, queryFiles, outName, TRUE, p->qType == gftDnaX,
            f, showStatus, p, gvo);
  } else {
    errAbort("Unrecognized combination of target and query types\n");
  }
  freeDnaSeqList(&dbSeqList);
}

