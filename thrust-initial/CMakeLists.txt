cmake_minimum_required(VERSION 3.16)
project(ThrustInitial CUDA)

find_package(CUDA REQUIRED)

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -g -gencode arch=compute_61,code=sm_61 -DCUDA_CUDART_LIBRARY=/usr/lib/x86_64-linux-gnu/libcudart.so.10.1.243)

add_definitions(-DEIGEN_MATRIX_PLUGIN="Matrix_initializer_list.h" -DEIGEN_ARRAY_PLUGIN="Array_initializer_list.h")
SET(ThrustInitialSrc main.cu)

INCLUDE_DIRECTORIES ( "/usr/include/eigen3"
                      "../repos/stl_reader"
                      "../repos/eigen-initializer_list/src"
                      "../reference" )
link_directories("../reference")

add_executable(initial main.cu)
target_link_libraries(initial -lPreprocessBezierLib)
