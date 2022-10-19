#!/bin/bash


mkdir ${PREFIX}/bin
cmake .
make 
mv fastnc ${PREFIX}/bin

# mkdir ${PREFIX}/bin
# g++ -std=c++11 -O3 -fopenmp -march=native -mavx -mfma -o ${PREFIX}/bin/fastnc src/fastnc.cpp src/fastnc_opts.cpp src/common.cpp -DMKL_ILP64 -m64\
 # ${PREFIX}/lib/libmkl_scalapack_ilp64.a \
 # -Wl,--start-group ${PREFIX}/lib/libmkl_cdft_core.a \
 # ${PREFIX}/lib/libmkl_intel_ilp64.a \
 # ${PREFIX}/lib/libmkl_gnu_thread.a \
 # ${PREFIX}/lib/libmkl_core.a\
 # ${PREFIX}/lib/libmkl_blacs_openmpi_ilp64.a \
 # -Wl,--end-group -lgomp -lpthread -lm -ldl
