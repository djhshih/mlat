/* mlat - Minimal BLAT fast sequence search command line tool. */
/* Copyright 2001-2004 Jim Kent.  All rights reserved. */
#include "mlat.h"

void usage()
/* Explain usage and exit. */
{
  printf(
      "blat - Standalone BLAT v. %s fast sequence search command line tool\n"
      "usage:\n"
      "   blat database query [-ooc=11.ooc] output.psl\n"
      "where:\n"
      "   database and query are each either a .fa or .2bit file,\n"
      "   or a list these files one file name per line.\n"
      "   -ooc=11.ooc tells the program to load over-occurring 11-mers from\n"
      "               and external file.  This will increase the speed\n"
      "               by a factor of 40 in many cases, but is not required\n"
      "   output.psl is where to put the output.\n"
      "   Subranges of .2bit files may specified using the syntax:\n"
      "      /path/file.2bit:seqid:start-end\n"
      "   With the second form, a sequence id of file:start-end will be used.\n"
      "options:\n"
      "   -t=type     Database type.  Type is one of:\n"
      "                 dna - DNA sequence\n"
      "                 prot - protein sequence\n"
      "                 dnax - DNA sequence translated in six frames to "
      "protein\n"
      "               The default is dna\n"
      "   -q=type     Query type.  Type is one of:\n"
      "                 dna - DNA sequence\n"
      "                 rna - RNA sequence\n"
      "                 prot - protein sequence\n"
      "                 dnax - DNA sequence translated in six frames to "
      "protein\n"
      "                 rnax - DNA sequence translated in three frames to "
      "protein\n"
      "               The default is dna\n"
      "   -prot       Synonymous with -t=prot -q=prot\n"
      "   -ooc=N.ooc  Use overused tile file N.ooc.  N should correspond to \n"
      "               the tileSize\n"
      "   -tileSize=N sets the size of match that triggers an alignment.  \n"
      "               Usually between 8 and 12\n"
      "               Default is 11 for DNA and 5 for protein.\n"
      "   -stepSize=N spacing between tiles. Default is tileSize.\n"
      "   -oneOff=N   If set to 1 this allows one mismatch in tile and still\n"
      "               triggers an alignments.  Default is 0.\n"
      "   -minMatch=N sets the number of tile matches.  Usually set from 2 to "
      "4\n"
      "               Default is 2 for nucleotide, 1 for protein.\n"
      "   -minScore=N sets minimum score.  This is the matches minus the \n"
      "               mismatches minus some sort of gap penalty.  Default is "
      "30\n"
      "   -minIdentity=N Sets minimum sequence identity (in percent).  Default "
      "is\n"
      "               90 for nucleotide searches, 25 for protein or "
      "translated\n"
      "               protein searches.\n"
      "   -maxGap=N   sets the size of maximum gap between tiles in a clump.  "
      "Usually\n"
      "               set from 0 to 3.  Default is 2. Only relevent for "
      "minMatch > 1.\n"
      "   -noHead     suppress .psl header (so it's just a tab-separated "
      "file)\n"
      "   -makeOoc=N.ooc Make overused tile file. Target needs to be complete "
      "genome.\n"
      "   -repMatch=N sets the number of repetitions of a tile allowed before\n"
      "               it is marked as overused.  Typically this is 256 for "
      "tileSize\n"
      "               12, 1024 for tile size 11, 4096 for tile size 10.\n"
      "               Default is 1024.  Typically only comes into play with "
      "makeOoc.\n"
      "               Also affected by stepSize. When stepSize is halved "
      "repMatch is\n"
      "               doubled to compensate.\n"
      "   -mask=type  Mask out repeats.  Alignments won't be started in masked "
      "region\n"
      "               but may extend through it in nucleotide searches.  "
      "Masked areas\n"
      "               are ignored entirely in protein or translated searches. "
      "Types are\n"
      "                 lower - mask out lower cased sequence\n"
      "                 upper - mask out upper cased sequence\n"
      "                 out   - mask according to database.out RepeatMasker "
      ".out file\n"
      "                 file.out - mask database according to RepeatMasker "
      "file.out\n"
      "   -qMask=type Mask out repeats in query sequence.  Similar to -mask "
      "above but\n"
      "               for query rather than target sequence.\n"
      "   -repeats=type Type is same as mask types above.  Repeat bases will "
      "not be\n"
      "               masked in any way, but matches in repeat areas will be "
      "reported\n"
      "               separately from matches in other areas in the psl "
      "output.\n"
      "   -minRepDivergence=NN - minimum percent divergence of repeats to "
      "allow \n"
      "               them to be unmasked.  Default is 15.  Only relevant for "
      "\n"
      "               masking using RepeatMasker .out files.\n"
      "   -trimT      Trim leading poly-T\n"
      "   -noTrimA    Don't trim trailing poly-A\n"
      "   -trimHardA  Remove poly-A tail from qSize as well as alignments in \n"
      "               psl output\n"
      "   -fastMap    Run for fast DNA/DNA remapping - not allowing introns, \n"
      "               requiring high %%ID. Query sizes must not exceed %d.\n"
      "   -out=type   Controls output file format.  Type is one of:\n"
      "                   psl - Default.  Tab separated format, no sequence\n"
      "                   pslx - Tab separated format with sequence\n"
      "                   axt - blastz-associated axt format\n"
      "                   maf - multiz-associated maf format\n"
      "                   sim4 - similar to sim4 format\n"
      "                   wublast - similar to wublast format\n"
      "                   blast - similar to NCBI blast format\n"
      "                   blast8- NCBI blast tabular format\n"
      "                   blast9 - NCBI blast tabular format with comments\n"
      "   -fine       For high quality mRNAs look harder for small initial "
      "and\n"
      "               terminal exons.  Not recommended for ESTs\n"
      "   -maxIntron=N  Sets maximum intron size. Default is %d\n"
      "   -extendThroughN - Allows extension of alignment through large blocks "
      "of N's\n",
      gfVersion, MAXSINGLEPIECESIZE, ffIntronMaxDefault);
  exit(-1);
}

struct optionSpec options[] = {
    {"t", OPTION_STRING},
    {"q", OPTION_STRING},
    {"prot", OPTION_BOOLEAN},
    {"ooc", OPTION_STRING},
    {"tileSize", OPTION_INT},
    {"stepSize", OPTION_INT},
    {"oneOff", OPTION_INT},
    {"minMatch", OPTION_INT},
    {"minScore", OPTION_INT},
    {"minIdentity", OPTION_FLOAT},
    {"maxGap", OPTION_INT},
    {"noHead", OPTION_BOOLEAN},
    {"makeOoc", OPTION_STRING},
    {"repMatch", OPTION_INT},
    {"mask", OPTION_STRING},
    {"qMask", OPTION_STRING},
    {"repeats", OPTION_STRING},
    {"minRepDivergence", OPTION_FLOAT},
    {"trimT", OPTION_BOOLEAN},
    {"noTrimA", OPTION_BOOLEAN},
    {"trimHardA", OPTION_BOOLEAN},
    {"fastMap", OPTION_BOOLEAN},
    {"out", OPTION_STRING},
    {"fine", OPTION_BOOLEAN},
    {"maxIntron", OPTION_INT},
    {"extendThroughN", OPTION_BOOLEAN},
    {NULL, 0},
};

int main(int argc, char *argv[])
/* Process command line into global variables and call blat. */
{
  boolean tIsProtLike, qIsProtLike;
  struct mlatParams *p = newMlatParams();

#ifdef DEBUG
  {
    char *cmd = "blat hCrea.geno hCrea.mrna foo.psl -t=dnax -q=rnax";
    char *words[16];

    printf("Debugging parameters\n");
    cmd = cloneString(cmd);
    argc = chopLine(cmd, words);
    argv = words;
  }
#endif /* DEBUG */

  optionInit(&argc, argv, options);
  if (argc != 4)
    usage();

  /* Get database and query sequence types and make sure they are
   * legal and compatable. */
  if (optionExists("prot"))
    p->qType = p->tType = gftProt;
  if (optionExists("t"))
    p->tType = gfTypeFromName(optionVal("t", NULL));
  p->trimA = optionExists("trimA") || optionExists("trima");
  p->trimT = optionExists("trimT") || optionExists("trimt");
  p->trimHardA = optionExists("trimHardA");
  switch (p->tType) {
  case gftProt:
  case gftDnaX:
    tIsProtLike = TRUE;
    break;
  case gftDna:
    tIsProtLike = FALSE;
    break;
  default:
    tIsProtLike = FALSE;
    errAbort("Illegal value for 't' parameter");
    break;
  }
  if (optionExists("q"))
    p->qType = gfTypeFromName(optionVal("q", NULL));
  if (p->qType == gftRnaX || p->qType == gftRna)
    p->trimA = TRUE;
  if (optionExists("noTrimA"))
    p->trimA = FALSE;
  switch (p->qType) {
  case gftProt:
  case gftDnaX:
  case gftRnaX:
    p->minIdentity = 25;
    qIsProtLike = TRUE;
    break;
  default:
    qIsProtLike = FALSE;
    break;
  }
  if ((tIsProtLike ^ qIsProtLike) != 0)
    errAbort("t and q must both be either protein or dna");

  /* Set default tile size for protein-based comparisons. */
  if (tIsProtLike) {
    p->tileSize = 5;
    p->minMatch = 1;
    p->oneOff = FALSE;
    p->maxGap = 0;
  }

  /* Get tile size and related parameters from user and make sure
   * they are within range. */
  p->tileSize = optionInt("tileSize", p->tileSize);
  p->stepSize = optionInt("stepSize", p->tileSize);
  p->minMatch = optionInt("minMatch", p->minMatch);
  p->oneOff = optionExists("oneOff");
  p->fastMap = optionExists("fastMap");
  p->fine = optionExists("fine");
  p->minScore = optionInt("minScore", p->minScore);
  p->maxGap = optionInt("maxGap", p->maxGap);
  p->minRepDivergence = optionFloat("minRepDivergence", p->minRepDivergence);
  p->minIdentity = optionFloat("minIdentity", p->minIdentity);
  gfCheckTileSize(p->tileSize, tIsProtLike);
  if (p->minMatch < 0)
    errAbort("minMatch must be at least 1");
  if (p->maxGap > 100)
    errAbort("maxGap must be less than 100");

  /* Set repMatch parameter from command line, or
   * to reasonable value that depends on tile size. */
  if (optionExists("repMatch"))
    p->repMatch = optionInt("repMatch", p->repMatch);
  else
    p->repMatch = gfDefaultRepMatch(p->tileSize, p->stepSize, tIsProtLike);

  /* Gather last few command line options. */
  p->noHead = optionExists("noHead");
  p->ooc = optionVal("ooc", NULL);
  p->makeOoc = optionVal("makeOoc", NULL);
  p->mask = optionVal("mask", NULL);
  p->qMask = optionVal("qMask", NULL);
  p->repeats = optionVal("repeats", NULL);
  if (p->repeats != NULL && p->mask != NULL &&
      differentString(p->repeats, p->mask))
    errAbort("The -mask and -repeat settings disagree.  "
             "You can just omit -repeat if -mask is on");
  if (p->mask != NULL) /* Mask setting will also set repeats. */
    p->repeats = p->mask;
  p->outputFormat = optionVal("out", p->outputFormat);
  /* set global for fuzzy find functions */
  setFfIntronMax(optionInt("maxIntron", ffIntronMaxDefault));
  setFfExtendThroughN(optionExists("extendThroughN"));

  /* Call routine that does the work. */
  mlat(argv[1], argv[2], argv[3], p);

  /* No members of p are dynamically allocated */
  freez(&p);

  return 0;
}
