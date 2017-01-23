ROOT=.

include $(ROOT)/common.mk

targets=mlat mlat-demo 2bit blatc blatd

all: $(targets)
	mkdir -p build/{bin,include,lib}
	cp lib/*.{a,so,dylib} build/lib || true
	cp include/*.h build/include
	cp -r include/mlat build/include
	mv $(targets) build/bin

lib/libmlat.a:
	cd lib && make libmlat.a

lib/libmlatnet.a:
	cd lib && make libmlatnet.a

mlat: lib/libmlat.a
	$(CBUILD) src/mlat.c lib/libmlat.a -o mlat

mlat-demo:
	$(CC) $(CFLAGS) -Iinclude -Llib -lmlat src/mlat-demo.c -o mlat-demo

mlat-shared:
	$(CBUILD) src/mlat.c -Llib -lmlat -o mlat-shared

2bit:
	$(CBUILD) src/2bit.c lib/core/*.o lib/aux/*.o -o 2bit

blatc: lib/libmlatnet.a
	$(CBUILD) src/blatc.c lib/libmlatnet.a -o blatc

blatd: lib/libmlatnet.a
	$(CBUILD) src/blatd.c lib/libmlatnet.a -o blatd

clean:
	rm -f *.{o,a,gcda,gcno,gcov} lib/{.,core,aux,net}/*.{o,a,gcda,gcno,gcov}
	rm -rf build

