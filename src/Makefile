ROOT=..

include $(ROOT)/common.mk

all: mlat 2bit
	

server: blatd blatc 
	

mlat: $(ROOT)/lib/libmlat.a
	$(CBUILD) mlat.c -L$(ROOT)/lib -lmlat -lm -o mlat

mlat-shared: $(ROOT)/lib/libmlat.so
	$(CBUILD) mlat.c -L$(ROOT)/lib -lmlat -lm -o mlat-shared

mlat-static: $(ROOT)/lib/libmlat.a
	$(CBUILD) mlat.c -L$(ROOT)/lib -lmlat -lm -lpthread -static -static-libgcc -D_STATIC -D_FORTIFY_SOURCE=0

2bit: $(ROOT)/lib/libmlat.a
	$(CBUILD) 2bit.c -L$(ROOT)/lib -lmlat -lm -o 2bit

blatc: $(ROOT)/lib/libmlatnet.a
	$(CBUILD) blatc.c -L$(ROOT)/lib -lmlatnet -lm -lpthread -o blatc

blatd: $(ROOT)/lib/libmlatnet.a
	$(CBUILD) blatd.c -L$(ROOT)/lib -lmlatnet -lm -lpthread -o blatd

clean:
	rm -f *.{gcno,gcov,gcda,log}
	rm -f mlat mlat-shared mlat-static 2bit blatc blatd

