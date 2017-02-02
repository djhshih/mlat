ROOT=..

include $(ROOT)/common.mk

mlat_objs=core/*.o aux/*.o
mlatnet_objs=$(mlat_objs) net/*.o

libmlat.a:
	cd core && make
	cd aux && make
	ar rcs libmlat.a $(mlat_objs)

libmlatnet.a:
	cd core && make
	cd aux && make
	cd net && make
	ar rcs libmlatnet.a $(mlatnet_objs)

libmlat.so:
	$(CBUILD) -shared -o libmlat.so $(mlat_objs)

libmlat.dylib:
	$(CBUILD) -dynamiclib -o libmlat.dylib $(mlat_objs)

clean:
	cd core && make clean
	cd aux && make clean
	cd net && make clean
	rm -f *.{a,so,dylib}

