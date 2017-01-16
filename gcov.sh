#!/bin/sh

gcov mlat.c > gcov.log
cd lib/core && gcov *.c > gcov.log && cd -
cd lib/aux && gcov *.c > gcov.log && cd -
