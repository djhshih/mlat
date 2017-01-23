ROOT=..

include $(ROOT)/common.mk

mlat_objs=core/*.o aux/*.o
mlatnet_objs=$(mlat_objs) net/*.o

libmlat.a:
	cd core && make
	cd aux && make
	ar rcus libmlat.a $(mlat_objs)

libmlatnet.a:
	cd core && make
	cd aux && make
	cd net && make
	ar rcus libmlatnet.a $(mlatnet_objs)

libmlat.so: mlat.o
	$(CBUILD) -shared -o libmlat.so $(mlat_objs)

libmlat.dylib: mlat.o
	$(CBUILD) -dynamiclib -o libmlat.dylib $(mlat_objs)

