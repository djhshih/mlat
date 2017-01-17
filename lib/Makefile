ROOT=..

include $(ROOT)/common.mk

mlat_objs=core/*.o aux/*.o mlat.o

libmlat.a: mlat.o
	cd core && make
	cd aux && make
	ar rcus libmlat.a $(mlat_objs)

libmlat.so: mlat.o
	$(CBUILD) -shared -o libmlat.so $(mlat_objs)

libmlat.dylib: mlat.o
	$(CBUILD) -dynamiclib -o libmlat.dylib $(mlat_objs)

