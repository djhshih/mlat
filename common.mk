COPTS=-I$(ROOT)/include/core -I$(ROOT)/include/aux -I$(ROOT)/include/net
CFLAGS=-coverage -O0

%.o: %.c
	$(CC) $(COPTS) $(CFLAGS) -o $@ -c $<

