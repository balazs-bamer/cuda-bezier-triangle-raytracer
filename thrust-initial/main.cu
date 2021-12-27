#define UTIL_CUDA_PREFIX_DEVICE __device__
#define UTIL_CUDA_PREFIX_HOST   __host__

#include "mesh.h"

#include<iomanip>
#include<iostream>


std::string const cgBaseDir("output/");

void testSharedLib() {
  std::string name{"testSharedLib.stl"};

  Mesh sphere;
  sphere.makeUnitSphere(3, 1);
  sphere.writeMesh(cgBaseDir + name);
}

int main(int argc, char **argv) {
  testSharedLib();
  return 0;
}
