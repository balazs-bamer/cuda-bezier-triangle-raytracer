#include "meshUtils.h"

#include<deque>
#include<vector>

using Real = float;

template<template<typename...> typename tContainer>
auto visualizeNormals(Mesh<Real, tContainer> const &aMesh) {
  Mesh<Real, std::deque> result;
  Mesh<Real, std::deque> sphere;
  sphere.makeUnitSphere(3, 1);
  for(auto &face : aMesh) {
    auto normal = Mesh<Real, std::deque>::getNormal(face).normalized();
    auto base = (face[0] + face[1] + face[2]) / 3.0f;
    auto average = ((face[1] - face[0]).norm() + (face[2] - face[1]).norm() + (face[0] - face[2]).norm()) / 3.0f;
    auto copy = sphere;
    auto shrink = Mesh<Real, std::deque>::Transform::Identity() * (average / 15.0f);
    copy.transform(shrink, base + normal * (average / 7.0f));
    std::copy(copy.cbegin(), copy.cend(), std::back_inserter(result));
  }
  return result;
}

void testDequeDivisor(char const * const aName, int32_t const aSectors, int32_t const aBelts, float const aRadius, int32_t const aDivisor) {
  std::string name{"test_"};
  name += aName;
  name += ".stl";

  Mesh<Real, std::deque> sphere;
  sphere.makeUnitSphere(aSectors, aBelts);
  Mesh<Real, std::deque>::Transform inflate = Mesh<Real, std::deque>::Transform::Identity() * aRadius;
  sphere *= inflate;
  sphere.writeMesh(name);

  Mesh<Real, std::deque> back;
  back.readMesh(name);
  back *= inflate;

  back.splitTriangles(aDivisor);
  Mesh<Real, std::deque>::Vector disp = Mesh<Real, std::deque>::Vector::Zero();
  disp[0] = aRadius;
  back += disp;
  auto nameBack = "back_" + name;
  back.standardizeVertices();
  back.standardizeNormals();
  back.writeMesh(nameBack);

  auto nameNorm = "norm_" + name;
  visualizeNormals<std::deque>(back).writeMesh(nameNorm);
}

void testVectorMax(char const * const aName, int32_t const aSectors, int32_t const aBelts, float const aRadius, float const aMaxSide) {
  std::string name{"test_"};
  name += aName;
  name += ".stl";

  Mesh<Real, std::vector> sphere;
  sphere.makeUnitSphere(aSectors, aBelts);
  Mesh<Real, std::vector>::Transform inflate = Mesh<Real, std::vector>::Transform::Identity() * aRadius;
  sphere *= inflate;
  sphere.writeMesh(name);

  Mesh<Real, std::vector> back;
  back.readMesh(name);
  back *= inflate;

  back.splitTriangles(aMaxSide);
  Mesh<Real, std::vector>::Vector disp = Mesh<Real, std::deque>::Vector::Zero();
  disp[0] = aRadius;
  back += disp;
  auto nameBack = "back_" + name;
  back.standardizeVertices();
  back.standardizeNormals();
  back.writeMesh(nameBack);

  auto nameNorm = "norm_" + name;
  visualizeNormals<std::vector>(back).writeMesh(nameNorm);
}

void testCustomStl(char * const aName) {
  Mesh<Real, std::vector> mesh;
  mesh.readMesh(aName);

  mesh.standardizeVertices();
  mesh.standardizeNormals();

  std::string nameBack{"back_"};
  nameBack += aName;
  mesh.writeMesh(nameBack);

  std::string nameNorm{"norm_"};
  nameNorm += aName;
  visualizeNormals<std::vector>(mesh).writeMesh(nameNorm);
}

int main(int argc, char **argv) {
  testDequeDivisor("dequeDivisor", 7, 7, 3.0f, 3);
  testVectorMax("vectorMax", 3, 1, 13.0f, 11.0f);

  if(argc > 1) {
    testCustomStl(argv[1]);
  }
  else { // Nothing to do
  }
  return 0;
}
