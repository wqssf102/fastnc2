# FastNC++：Fast calculation of natural connectivity
## Installing
FastNC++ can be installed using conda or from source.
## Conda
To install through conda, use:
```
conda create -n fastnc++ -c wqssf102 fastnc++
```
## Compiling from source
```
Compiling from source requires these libraries and software:
```
```
C++11 (gcc-4.9.0+, clang-4.9.0+, etc)
OpenMP 4.0+
Boost
mkl
g++
```
```
mklroot=path_mkl
boostroot=path_boos
g++ -std=c++11 -O3 -fopenmp -march=native -mavx -mfma -o fastnc++ fastnc++.cpp common.cpp -DMKL_ILP64 -m64 \
 -Wl,--start-group \
 ${mklroot}/lib/libmkl_scalapack_ilp64.a \
 ${mklroot}/lib/libmkl_cdft_core.a \
 ${mklroot}/lib/libmkl_intel_ilp64.a \
 ${mklroot}/lib/libmkl_gnu_thread.a \
 ${mklroot}/lib/libmkl_core.a \
 ${mklroot}/lib/libmkl_blacs_openmpi_ilp64.a \
  -Wl,--end-group -lgomp -lpthread -lm -ldl -I${mklroot}/include/ -L${mklroot}/lib \
  -static -static-libgcc -static-libstdc++ -lboost_program_options -I${boostroot} \
  -L${boostroot}/stage/lib
```
**help**：
```
fastnc++ --h
Program: FastNC++ (use c++ to calculate the natural connectivity)
Contact: Qiusheng WU (565715597@qq.com)

  -i [ --input ] arg          input.txt,The result from get.adjacency function of igraph package
  -o [ --output ] arg         output.txt
  -s [ --step ] arg (=1)      number of threads (default 1)
  -t [ --thread ] arg (=0.80) the threshold for deletion of node (default: 0.8)
  -n [ --number ] arg (=1000)  number of iterations (default: 1000)
  -j [ --jobs ] arg (=4)       number of jobs (default: 4)
  ```

## USAGE and Tutorial

* Chinese version

  <u>https://wqssf102.github.io/fastnc++/</u>

* English version

  In the plan
  
  
