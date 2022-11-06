# FastNC：Fast calculation of natural connectivity
## Installing
FastNC can be installed using conda or from source.
## Conda
To install through conda, use:
```
conda create -n fastnc -c wqssf102 -c conda-forge -c intel fastnc -y
```
## Compiling from source
```
Compiling from source requires these libraries and software:
```
```
C++11 (gcc-4.9.0+, clang-4.9.0+, etc)
OpenMP 4.0+
Eigen
mkl
g++
```
```
MKLROOT=pathA
Eigendir=pathB
g++ -std=c++11 -O3 -fopenmp -march=native -mavx -mfma -o fastnc src/fastnc.cpp src/fastnc_opts.cpp src/common.cpp -DMKL_ILP64 -m64\
 -I ${Eigendir}/include \
 -Wl,--start-group ${MKLROOT}/lib/libmkl_cdft_core.a \
 ${MKLROOT}/lib/libmkl_scalapack_ilp64.a \
 ${MKLROOT}/lib/libmkl_intel_ilp64.a \
 ${MKLROOT}/lib/libmkl_gnu_thread.a \
 ${MKLROOT}/lib/libmkl_core.a\
 ${MKLROOT}/lib/libmkl_blacs_openmpi_ilp64.a \
 -Wl,--end-group -lgomp -lpthread -lm -ldl
```
**help**：
```
fastnc --h
Program: FastNC (use c++ to calculate the natural connectivity)
Contact: Qiusheng WU (565715597@qq.com)

  fastnc [options] --adj_table <file> --outfile <file>
  -c <file>, --adj_table <file>
                The result from get.adjacency function of igraph package
  -o <file>, -outfile <file>
                Result output table
Options:
  -t <float>, -threshold <float>
                The threshold for deletion of node (default: 0.8)
  -s <int>, -step <int>
                The step for deletion of node (default: 1)
  -g <float>, -edge <float>
                The threshold for deletion of edge (default: 0)
  -n <int>, -number <int>
                Number of iterations (default: 1000)
  -j <int>, -job <int>
                Number of jobs (default: 4)
Other:
  -h        --help
                Display this help and exit
  ```

## USAGE and Tutorial

* Chinese version

  <u>https://wqssf102.github.io/fastnc/</u>

* English version

  In the plan
  
  
