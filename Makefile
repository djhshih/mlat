ROOT=.

include $(ROOT)/common.mk

targets=mlat mlat-alt mlat-shared 2bit blatc blatd

all: $(targets)
	mkdir -p bin
	mv $(targets) bin

lib/libmlat.a:
	cd lib && make libmlat.a

mlat: lib/libmlat.a
	$(CBUILD) mlat.c lib/libmlat.a -o mlat

mlat-alt:
	$(CBUILD) mlat-alt.c lib/core/*.o lib/aux/*.o lib/mlat.o lib/mlat-alt.o \
		-o mlat-alt

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

