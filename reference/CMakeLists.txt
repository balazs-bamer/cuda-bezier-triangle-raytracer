cmake_minimum_required(VERSION 3.16)
project(BezierTrangleReference)

enable_testing()
find_package(GTest REQUIRED)

#add_compile_options(-ggdb -D_GLIBCXX_DEBUG -std=c++17)
add_compile_options(-O2 -std=c++17)

add_definitions(-DEIGEN_MATRIX_PLUGIN="Matrix_initializer_list.h" -DEIGEN_ARRAY_PLUGIN="Array_initializer_list.h")
SET(PreprocessBezierSrc mesh.cpp bezierTriangle.cpp bezierMesh.cpp bezierLens.cpp)

INCLUDE_DIRECTORIES ( "/usr/include/eigen3"
                      "../repos/stl_reader"
                      "../repos/eigen-initializer_list/src" )
ADD_LIBRARY (PreprocessBezierLib SHARED ${PreprocessBezierSrc})

add_executable(googleTest googleTest.cpp)
target_link_libraries(googleTest GTest::GTest GTest::Main)
add_test(google-test googleTest)

add_executable(manual-test test.cpp)
target_link_libraries(manual-test PreprocessBezierLib)
