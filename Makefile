ROOT=.

include $(ROOT)/common.mk

bin:
	mkdir bin

mlat: 
	$(CC) $(COPTS) $(CFLAGS) mlat.c lib/core/*.o lib/aux/*.o -o mlat

2bit:
	$(CC) $(COPTS) $(CFLAGS) 2bit.c lib/core/*.o lib/aux/*.o -o 2bit

blatc:
	$(CC) $(COPTS) $(CFLAGS) blatc.c lib/core/*.o lib/aux/*.o lib/net/*.o -o blatc

blatd:
	$(CC) $(COPTS) $(CFLAGS) blatd.c lib/core/*.o lib/aux/*.o lib/net/*.o -o blatd

clean:
	rm -f *.{o,gcda,gcno,gcov} lib/{core,aux,net}/*.{o,gcda,gcno,gcov}

