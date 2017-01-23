ROOT=.

include $(ROOT)/common.mk

targets=mlat mlat-alt mlat-shared 2bit blatc blatd

all: $(targets)
	mkdir -p bin
	mv $(targets) bin

lib/libmlat.a:
	cd lib && make libmlat.a

lib/libmlatnet.a:
	cd lib && make libmlatnet.a

mlat: lib/libmlat.a
	$(CBUILD) mlat.c lib/libmlat.a -o mlat

mlat-alt:
	$(CBUILD) mlat-alt.c lib/core/*.o lib/aux/*.o lib/mlat-alt.o \
		-o mlat-alt

mlat-shared:
	$(CBUILD) mlat.c -Llib -lmlat -o mlat-shared

2bit:
	$(CBUILD) 2bit.c lib/core/*.o lib/aux/*.o -o 2bit

blatc: lib/libmlatnet.a
	$(CBUILD) blatc.c lib/libmlatnet.a -o blatc

blatd: lib/libmlatnet.a
	$(CBUILD) blatd.c lib/libmlatnet.a -o blatd

clean:
	rm -f *.{o,a,gcda,gcno,gcov} lib/{.,core,aux,net}/*.{o,a,gcda,gcno,gcov} bin/*

