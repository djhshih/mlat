#include "mlat.h"

class Mlat
{

	mlatParams* p;
	gfDb* db;
	gfResult* res;

public:

	Mlat(char *dbFname)
		: p(newMlatParams()), db(newGfDb(dbFname, p)), res(NULL)
	{}
	
	~Mlat() {
		freeMlatParams(&p);
		freeGfDb(&db);
		freeGfResult(&res);
	}

	const mlatParams& params() const {
		return *p;
	}

	mlatParams& params(){
		return *p;
	}

	gfResult* search(char* querySeq) {
		freeGfResult(&res);
		res = searchSeq(db, querySeq, p);
		return res;
	}

};
