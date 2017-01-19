#include "mlat-alt.h"

int main(int argc, char *argv[]) {
  struct mlatParams *p = newMlatParams();
  p->tileSize = 6;
  p->stepSize = 2;
  p->minScore = 0;

  char *dbFile = argv[1];
  struct gfDb *db = newGfDb(dbFile, p);

  fprintf(stderr, "INFO: allocated gfDb in %s\n", __func__);
  fprintf(stderr, "\n");

  size_t queryCount = argc - 2;
  for (size_t i = 0; i < queryCount; ++i) {
    fprintf(stderr, "INFO: query %zu\n", i);

    char *querySeq = argv[i + 2];
    struct gfOutput *out = searchSeq(db, querySeq, p);

    fprintf(stderr, "INFO: got gfOutput *out == %p in %s\n", out, __func__);

    struct gfResult *res = out->data;

    if (res->size == 0) {
      // no match: move onto next query
      continue;
    }

    struct gfAlign *hit = &res->aligns[0];

    fprintf(stderr, "INFO: got gfAlign *hit == %p in %s\n", hit, __func__);

    fprintf(stderr, "INFO: got res->size == %zu in %s\n", res->size, __func__);
    fprintf(stderr, "INFO: got res->capacity == %zu in %s\n", res->capacity,
            __func__);
    fprintf(stderr, "INFO: got res->tName == %p in %s\n", res->tName, __func__);
    fprintf(stderr, "INFO: got res->qName == %p in %s\n", res->qName, __func__);

    fprintf(stdout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "matchCount",
            "mismatchCount", "repMatchCount", "nCount", "qInsertCount",
            "qInsertBaseCount", "tInsertCount", "tInsertBaseCount",
            "qStrand/tStrand");
    fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c%c\n", hit->matchCount,
            hit->mismatchCount, hit->repMatchCount, hit->nCount,
            hit->qInsertCount, hit->qInsertBaseCount, hit->tInsertCount,
            hit->tInsertBaseCount, hit->qStrand, hit->tStrand);

    fprintf(stderr, "INFO: freeing out in %s\n", __func__);
    freeGfOutputResult(&out);
    fprintf(stderr, "INFO: freed out in %s\n", __func__);
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "INFO: finished queries in %s\n", __func__);

  freeGfDb(&db);

  fprintf(stderr, "INFO: freed db in %s\n", __func__);

  freez(&p);

  fprintf(stderr, "INFO: freed p in %s\n", __func__);

  return 0;
}
