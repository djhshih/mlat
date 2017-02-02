# MLAT: Minimal BLAT

[![travis-ci](https://travis-ci.org/djhshih/mlat.svg?branch=master)](https://travis-ci.org/djhshih/mlat)
[![codecov](https://codecov.io/gh/djhshih/mlat/branch/master/graph/badge.svg)](https://codecov.io/gh/djhshih/mlat)


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


# Known bugs

Code compilation was tested with 

```
gcc-4.4.7
gcc-5.4.0
gcc-6.3.0
clang-3.5
```

The code is not compatible with `-O2` and `-O3` optimization flags when compiled
with `gcc-4.4.7`. Use `-O2` or `-O3` only when compiling with `gcc-5.4.0` or
higher.

Unit tests were built using programs compiled with `clang-3.5`. Moderate
differences exist when compiled using a different compiler. Sometimes, programs
also produce moderate differences for different compilations with the same
compiler: the cause of these differences is inherited upstream from BLAT and
remains unresolved.

