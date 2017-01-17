COPTS=-I$(ROOT)/include -I$(ROOT)/include/core -I$(ROOT)/include/aux -I$(ROOT)/include/net
CFLAGS=-O3 -Os
#CFLAGS=-coverage -O0

CBUILD=$(CC) $(COPTS) $(CFLAGS)

%.o: %.c
	$(CBUILD) -o $@ -c $<

