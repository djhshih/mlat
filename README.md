# MLAT: Minimal BLAT

MLAT is a tool for sequence alignment.

This is a fork of the BLAT package (v35) written by Jim Kent. The package has
been significantly restructured and stripped down to minimize code footprint.
Header and source files have been modified, and minimal C and C++ API's are now
provided.

This package itself is licensed under GPLv3. Use of BLAT (and hence MLAT) is
subject to the terms of the BLAT non-commercial license (see
`LICENSE_blat.txt`).


# Installation

There are no dependencies other than GNU build tools (e.g. `gcc` and `make`),
and you can simply

    make

Thereupon, binaries, libraries, and headers will be installed in `build/`. To
install to some other location (e.g. `/usr/local`),

    make DESTDIR=/usr/local

