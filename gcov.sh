#!/bin/sh

cd src && gcov *.c > gcov.log && cd -
cd demo && gcov *.{c,cpp} > gcov.log && cd -
cd lib/core && gcov *.c > gcov.log && cd -
cd lib/aux && gcov *.c > gcov.log && cd -
