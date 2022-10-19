cmake_minimum_required(VERSION 3.13)
project(fastnc)
SET(CMAKE_C_COMPILER "gcc")
SET(CMAKE_CXX_COMPILER "g++")
include_directories(include)
link_directories(lib)
link_libraries("libmkl_cdft_core.a" 
"libmkl_scalapack_ilp64.a" 
"libmkl_intel_ilp64.a" 
"libmkl_gnu_thread.a" 
"libmkl_core.a" 
"libmkl_blacs_openmpi_ilp64.a" )

target_link_libraries(fastnc 
-Wl,--start-group
libmkl_cdft_core.a 
libmkl_scalapack_ilp64.a 
libmkl_intel_ilp64.a 
libmkl_gnu_thread.a 
libmkl_core.a 
libmkl_blacs_openmpi_ilp64.a
-Wl,--end-group
)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -fopenmp -march=native -mavx -mfma -DMKL_ILP64 -m64 -lgomp -lpthread -lm -ldl")

add_executable(fastnc
  ${CMAKE_SOURCE_DIR}/src/common.cpp
  ${CMAKE_SOURCE_DIR}/src/fastnc.cpp
  ${CMAKE_SOURCE_DIR}/src/fastnc_opts.cpp
)
# install(TARGETS fastnc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
