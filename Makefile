ROOT=.

include $(ROOT)/common.mk

bin:
	mkdir bin

mlat: 
	$(CC) $(COPTS) $(CFLAGS) mlat.c lib/core/*.o lib/aux/*.o -o mlat

2bit:
	$(CC) $(COPTS) $(CFLAGS) 2bit.c lib/core/*.o lib/aux/*.o -o 2bit

mlatc:
	$(CC) $(COPTS) $(CFLAGS) mlatc.c lib/core/*.o lib/aux/*.o lib/net/*.o -o mlatc

mlatd:
	$(CC) $(COPTS) $(CFLAGS) mlatd.c lib/core/*.o lib/aux/*.o lib/net/*.o -o mlatd

clean:
	rm -f *.{o,gcda,gcno,gcov} lib/{core,aux,net}/*.{o,gcda,gcno,gcov}

