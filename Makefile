ROOT=.

include $(ROOT)/common.mk

bin:
	mkdir bin

mlat: 
	$(CBUILD) mlat.c lib/core/*.o lib/aux/*.o -o mlat

2bit:
	$(CBUILD) 2bit.c lib/core/*.o lib/aux/*.o -o 2bit

blatc:
	$(CBUILD) blatc.c lib/core/*.o lib/aux/*.o lib/net/*.o -o blatc

blatd:
	$(CBUILD) blatd.c lib/core/*.o lib/aux/*.o lib/net/*.o -o blatd

clean:
	rm -f *.{o,gcda,gcno,gcov} lib/{core,aux,net}/*.{o,gcda,gcno,gcov}

