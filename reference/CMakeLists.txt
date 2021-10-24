cmake_minimum_required(VERSION 3.16)
project(BezierTrangleReference)

add_compile_options(-ggdb -D_GLIBCXX_DEBUG -std=c++17)

add_definitions(-DEIGEN_MATRIX_PLUGIN="Matrix_initializer_list.h" -DEIGEN_ARRAY_PLUGIN="Array_initializer_list.h")

INCLUDE_DIRECTORIES ( "/usr/include/eigen3"
                      "../repos/stl_reader"
                      "../repos/eigen-initializer_list/src" )

add_executable(test test.cpp)
add_executable(rayTracer main.cpp rayTracer.cpp )
