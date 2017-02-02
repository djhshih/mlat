ROOT = .

include $(ROOT)/common.mk

targets = mlat 2bit blatc blatd
headers = aliType.h gfDb.h gfResult.h mlatParams.h
DESTDIR ?= build

all: base
	

base:
	cd lib && make libmlat.a
	cd lib && make libmlatnet.a
	cd src && make

shared:
	cd lib && make libmlat.so

example:
	cd demo && make

check: base example
	./test.sh

coverage: base example check
	./gcov.sh

install: base
	mkdir -p $(DESTDIR)/{bin,include,include/mlat,lib}
	cp lib/libmlat.* $(DESTDIR)/lib || true
	cp include/mlat.{h,hpp} $(DESTDIR)/include
	cp -L $(addprefix include/mlat/, $(headers)) $(DESTDIR)/include/mlat
	cp -L $(addprefix src/, $(targets)) $(DESTDIR)/bin

uninstall:
	rm -f $(addprefix $(DESTDIR)/bin/, $(targets))
	rm -rf $(DESTDIR)/include/mlat
	rm -f $(DESTDIR)/include/mlat.{h,hpp}
	rm -f $(DESTDIR)/lib/libmlat.*

clean:
	cd lib && make clean
	cd src && make clean
	cd demo && make clean

