#! /bin/bash

g++ -g -fsanitize=address -fno-omit-frame-pointer -c *.cpp \
  -I /usr/include/hdf5/serial -fopenmp
g++ -g -fsanitize=address -fno-omit-frame-pointer -o test *.o \
  -L /usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -fopenmp

