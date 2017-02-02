ROOT=..

include $(ROOT)/common.mk

all: mlat-demo mlat-demo-cpp
	

mlat-demo: $(ROOT)/lib/libmlat.a
	$(CC) $(CFLAGS) -I$(ROOT)/include -L$(ROOT)/lib mlat-demo.c -lmlat -lm -o mlat-demo

mlat-demo-cpp: $(ROOT)/lib/libmlat.a
	$(CXX) $(CFLAGS) -I$(ROOT)/include -L$(ROOT)/lib mlat-demo-cpp.cpp -lmlat -lm -o mlat-demo-cpp

clean:
	rm -f *.{gcno,gcov,gcda,log}
	rm -f mlat-demo mlat-demo-cpp

