MACHTYPE?=$(shell uname -m)

COPTS=-I$(ROOT)/include -I$(ROOT)/include/private/core -I$(ROOT)/include/private/aux -I$(ROOT)/include/private/net -DMACHTYPE_${MACHTYPE}

# NB  code is not compatible with -O2 and -O3 with gcc-4.4.7
# NB  use -O2 or -O3 only with gcc-5.4.0 (tested) or later
CFLAGS_RELEASE=-O3 -Os

CFLAGS_DEBUG=-g

CFLAGS_COV=-coverage -O0

# disable fortify to circumvent missing __memchk in glibc
CFLAGS_STATIC=-O2 -Os -static -static-libgcc -D_STATIC -D_FORTIFY_SOURCE=0

CFLAGS=-O1

CBUILD=$(CC) $(COPTS) $(CFLAGS)

%.o: %.c
	$(CBUILD) -o $@ -c $<

