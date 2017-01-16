COPTS=-I$(ROOT)/include/core -I$(ROOT)/include/aux -I$(ROOT)/include/net
CFLAGS=-O3 -Os
#CFLAGS=-coverage -O0

%.o: %.c
	$(CC) $(COPTS) $(CFLAGS) -o $@ -c $<

