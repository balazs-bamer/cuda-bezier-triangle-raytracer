cmake_minimum_required(VERSION 3.16)
project(BezierTrangleReference)

add_compile_options(-ggdb -D_GLIBCXX_DEBUG -std=c++17)
INCLUDE_DIRECTORIES ( "/usr/include/eigen3"
                      "../repos/stl_reader" )

add_executable(test test.cpp)
add_executable(rayTracer main.cpp rayTracer.cpp )
