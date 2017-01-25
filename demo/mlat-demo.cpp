#include <stdio.h>
#include "mlat.hpp"

int main(int argc, char *argv[]) {

  if (argc < 2) {
    fprintf(stderr, "usage: %s <database> <sequence> [<sequence> ...]\n",
            argv[0]);
    return 1;
  }

  char *dbFile = argv[1];
  mlat::Database m(dbFile);

  mlatParams& p = m.params();
  p.tileSize = 6;
  p.stepSize = 2;
  p.minScore = 1;

  fprintf(stderr, "INFO: allocated gfDb in %s\n", __func__);
  fprintf(stderr, "\n");

  size_t queryCount = argc - 2;
  for (size_t i = 0; i < queryCount; ++i) {
    fprintf(stderr, "INFO: query %zu\n", i);

    char *querySeq = argv[i + 2];

    mlat::Result res = m.search(querySeq);

    fprintf(stderr, "INFO: got gfResult *res == %p in %s\n", &res, __func__);

    fprintf(stderr, "INFO: got res->size == %zu in %s\n", res.size(), __func__);
    fprintf(stderr, "INFO: got res->capacity == %zu in %s\n", res.capacity(),
            __func__);

    if (res.size() == 0) {
      // no match: move onto next query
      continue;
    }

    const gfAlign *align = &res[0];

    fprintf(stderr, "INFO: got gfAlign *align == %p in %s\n", align, __func__);

    fprintf(stderr, "INFO: got align->tName == %s in %s\n", align->tName, __func__);
    fprintf(stderr, "INFO: got align->qName == %s in %s\n", align->qName, __func__);

    fprintf(stdout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "matchCount",
            "mismatchCount", "repMatchCount", "nCount", "qInsertCount",
            "qInsertBaseCount", "tInsertCount", "tInsertBaseCount",
            "qStrand/tStrand");
    fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c%c\n", align->matchCount,
            align->mismatchCount, align->repMatchCount, align->nCount,
            align->qInsertCount, align->qInsertBaseCount, align->tInsertCount,
            align->tInsertBaseCount, align->qStrand, align->tStrand);

    fprintf(stderr, "\n");
  }

  return 0;
}
