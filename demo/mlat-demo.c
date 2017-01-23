#include <stdio.h>
#include "mlat.h"

int main(int argc, char *argv[]) {
  struct mlatParams *p = newMlatParams();
  p->tileSize = 6;
  p->stepSize = 2;
  p->minScore = 0;

  if (argc < 2) {
    fprintf(stderr, "usage: %s <database> <sequence> [<sequence> ...]\n",
            argv[0]);
    return 1;
  }

  char *dbFile = argv[1];
  struct gfDb *db = newGfDb(dbFile, p);

  fprintf(stderr, "INFO: allocated gfDb in %s\n", __func__);
  fprintf(stderr, "\n");

  size_t queryCount = argc - 2;
  for (size_t i = 0; i < queryCount; ++i) {
    fprintf(stderr, "INFO: query %zu\n", i);

    char *querySeq = argv[i + 2];
    struct gfResult *res = searchSeq(db, querySeq, p);

    fprintf(stderr, "INFO: got gfResult *res == %p in %s\n", res, __func__);

    if (res->size == 0) {
      // no match: move onto next query
      continue;
    }

    struct gfAlign *align = &res->aligns[0];

    fprintf(stderr, "INFO: got gfAlign *align == %p in %s\n", align, __func__);

    fprintf(stderr, "INFO: got res->size == %zu in %s\n", res->size, __func__);
    fprintf(stderr, "INFO: got res->capacity == %zu in %s\n", res->capacity,
            __func__);
    fprintf(stderr, "INFO: got res->tName == %p in %s\n", res->tName, __func__);
    fprintf(stderr, "INFO: got res->qName == %p in %s\n", res->qName, __func__);

    fprintf(stdout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "matchCount",
            "mismatchCount", "repMatchCount", "nCount", "qInsertCount",
            "qInsertBaseCount", "tInsertCount", "tInsertBaseCount",
            "qStrand/tStrand");
    fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c%c\n", align->matchCount,
            align->mismatchCount, align->repMatchCount, align->nCount,
            align->qInsertCount, align->qInsertBaseCount, align->tInsertCount,
            align->tInsertBaseCount, align->qStrand, align->tStrand);

    fprintf(stderr, "INFO: freeing res in %s\n", __func__);
    freeGfResult(&res);
    fprintf(stderr, "INFO: freed res in %s\n", __func__);

    fprintf(stderr, "\n");
  }

  fprintf(stderr, "INFO: finished queries in %s\n", __func__);

  freeGfDb(&db);

  fprintf(stderr, "INFO: freed db in %s\n", __func__);

  freeMlatParams(&p);

  fprintf(stderr, "INFO: freed p in %s\n", __func__);

  return 0;
}
