#include "mlat.h"
#include <algorithm>
#include <stdexcept>

#ifndef _MLAT_HPP_
#define _MLAT_HPP_

namespace mlat {

namespace detail {

bool gfAlignGreater(gfAlign &lhs, gfAlign &rhs) {
  return lhs.score > rhs.score;
}

}

class Result {

  gfResult *rep;

public:
  Result(gfResult *r) : rep(r) {
    // sort aligns by descending score
    std::sort(rep->aligns, rep->aligns + rep->size, detail::gfAlignGreater);
  }

  Result(const Result &that) : rep(cloneGfResult(that.rep)) {}

  Result & operator=(const Result &that) {
    if (this != &other) {
      freeGfResult(rep);
      rep = cloneGfResult(that.rep);
    }
    return *this;
  }

  const char *target() const { return rep->tName; }

  const char *query() const { return rep->qName; }

  const size_t size() const { return rep->size; }

  const size_t capacity() const { return rep->capacity; }

  // read-only access
  const gfAlign &operator[](size_t i) const {
    // no bounds checking, as per operator[] convention
    return rep->aligns[i];
  }

  // read-only access
  const gfAlign& at(size_t i) const {
    if (i < rep->size) {
      return rep->aligns[i];
    }
    throw std::out_of_range(__func__);
  }

  ~Result() { freeGfResult(&rep); }
};

class Database {

  mlatParams *p;
  gfDb *db;

public:
  // TODO const char
  Database(char *dbFname)
      : p(newMlatParams()), db(newGfDb(dbFname, p)) {}

  ~Database() {
    freeMlatParams(&p);
    freeGfDb(&db);
  }

  const mlatParams &params() const { return *p; }

  mlatParams &params() { return *p; }

  Result search(char *querySeq) {
    return Result(searchSeq(db, querySeq, p));
  }
};
}

#endif // _MLAT_HPP_
