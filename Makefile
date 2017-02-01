ROOT=.

include $(ROOT)/common.mk

targets=mlat 2bit blatc blatd
DESTDIR?=build

all: $(targets)
	mkdir -p $(DESTDIR)/{bin,include,lib}
	cp lib/libmlat.{a,so,dylib} $(DESTDIR)/lib || true
	cp include/*.{h,hpp} $(DESTDIR)/include
	cp -rL include/mlat $(DESTDIR)/include
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
	$(CBUILD) src/mlat.c lib/libmlat.a -lm -o mlat

mlat-shared: lib/libmlat.so
	$(CBUILD) src/mlat.c -Llib -lmlat -lm -o mlat-shared

mlat-static: lib/libmlat.a
	$(CBUILD) src/mlat.c -Llib -lmlat -lm -lpthread -static -static-libgcc -D_STATIC -D_FORTIFY_SOURCE=0

2bit:
	$(CBUILD) src/2bit.c lib/core/*.o lib/aux/*.o -lm -o 2bit

blatc: lib/libmlatnet.a
	$(CBUILD) src/blatc.c lib/libmlatnet.a -lm -lpthread -o blatc

blatd: lib/libmlatnet.a
	$(CBUILD) src/blatd.c lib/libmlatnet.a -lm -lpthread -o blatd

example: demo/mlat-demo demo/mlat-demo-cpp
	if [[ -f lib/libmlat.so ]]; then cp lib/libmlat.so demo; fi

demo/mlat-demo: lib/libmlat.a
	$(CC) $(CFLAGS) -Iinclude -Llib demo/mlat-demo.c -lmlat -lm -o demo/mlat-demo

demo/mlat-demo-cpp: lib/libmlat.a
	$(CXX) $(CFLAGS) -Iinclude -Llib demo/mlat-demo.cpp -lmlat -lm -o demo/mlat-demo-cpp

check: build/bin/mlat demo/mlat-demo demo/mlat-demo-cpp
	./test.sh

clean:
	rm -f *.{o,a,gcda,gcno,gcov} lib/{.,core,aux,net}/*.{o,a,gcda,gcno,gcov,so,dylib}
	rm -rf build
	rm -f demo/mlat-demo demo/mlat-demo-cpp

