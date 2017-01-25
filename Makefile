ROOT=.

include $(ROOT)/common.mk

targets=mlat 2bit blatc blatd
DESTDIR?=build

all: $(targets) demo/mlat-demo
	mkdir -p $(DESTDIR)/{bin,include,lib}
	cp lib/libmlat.{a,so,dylib} $(DESTDIR)/lib || true
	cp include/*.{h,hpp} $(DESTDIR)/include
	cp -r include/mlat $(DESTDIR)/include
	mv $(targets) $(DESTDIR)/bin

shared:
	cd lib && make libmlat.so
	cp lib/*.so build/lib

lib/libmlat.a:
	cd lib && make libmlat.a

lib/libmlatnet.a:
	cd lib && make libmlatnet.a

lib/libmlat.so:
	cd lib && make libmlat.so

mlat: lib/libmlat.a
	$(CBUILD) src/mlat.c lib/libmlat.a -o mlat

example: demo/mlat-demo demo/mlat-demo-cpp
	cp lib/libmlat.so demo

demo/mlat-demo: lib/libmlat.so
	$(CC) $(CFLAGS) -Iinclude -Llib -lmlat demo/mlat-demo.c -o demo/mlat-demo

demo/mlat-demo-cpp: lib/libmlat.so
	$(CXX) $(CFLAGS) -Iinclude -Llib -lmlat demo/mlat-demo.cpp -o demo/mlat-demo-cpp

mlat-shared:
	$(CBUILD) src/mlat.c -Llib -lmlat -o mlat-shared

2bit:
	$(CBUILD) src/2bit.c lib/core/*.o lib/aux/*.o -o 2bit

blatc: lib/libmlatnet.a
	$(CBUILD) src/blatc.c lib/libmlatnet.a -o blatc

blatd: lib/libmlatnet.a
	$(CBUILD) src/blatd.c lib/libmlatnet.a -o blatd

test: build/bin/mlat demo/mlat-demo
	./test.sh

clean:
	rm -f *.{o,a,gcda,gcno,gcov} lib/{.,core,aux,net}/*.{o,a,gcda,gcno,gcov,so,dylib}
	rm -rf build
	rm -f demo/mlat-demo demo/mlat-demo-cpp

