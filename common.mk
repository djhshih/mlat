COPTS=-I$(ROOT)/include -I$(ROOT)/include/private/core -I$(ROOT)/include/private/aux -I$(ROOT)/include/private/net
CFLAGS=-O3 -Os
#CFLAGS=-g
#CFLAGS=-coverage -O0

CBUILD=$(CC) $(COPTS) $(CFLAGS)

%.o: %.c
	$(CBUILD) -o $@ -c $<

