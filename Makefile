ROOT=.

include $(ROOT)/common.mk

all: mlat
	

lib/libmlat.a:
	cd lib && make libmlat.a

mlat: lib/libmlat.a
	$(CBUILD) mlat.c lib/libmlat.a -o mlat

mlat-shared:
	$(CBUILD) mlat.c -Llib -lmlat -o mlat-shared

2bit:
	$(CBUILD) 2bit.c lib/core/*.o lib/aux/*.o -o 2bit

blatc:
	$(CBUILD) blatc.c lib/core/*.o lib/aux/*.o lib/net/*.o -o blatc

blatd:
	$(CBUILD) blatd.c lib/core/*.o lib/aux/*.o lib/net/*.o -o blatd

clean:
	rm -f *.{o,gcda,gcno,gcov} lib/{core,aux,net}/*.{o,gcda,gcno,gcov}

