#!/bin/sh

clang-format -i -sort-includes=0 *.c
clang-format -i -sort-includes=0 lib/*.c
clang-format -i -sort-includes=0 lib/*/*.c
clang-format -i -sort-includes=0 include/*.h
clang-format -i -sort-includes=0 include/*/*.h
