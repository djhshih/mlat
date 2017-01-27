MACHTYPE?=$(shell uname -m)

COPTS=-I$(ROOT)/include -I$(ROOT)/include/private/core -I$(ROOT)/include/private/aux -I$(ROOT)/include/private/net -DMACHTYPE_${MACHTYPE}

CFLAGS_RELEASE=-O3 -Os
CFLAGS_DEBUG=-g
CFLAGS_COV=-coverage -O0
# disable fortify to circumvent missing __memchk in glibc
CFLAGS_STATIC=-O2 -Os -static -static-libgcc -D_STATIC -D_FORTIFY_SOURCE=0

CFLAGS=$(CFLAGS_RELEASE)

CBUILD=$(CC) $(COPTS) $(CFLAGS)

%.o: %.c
	$(CBUILD) -o $@ -c $<

