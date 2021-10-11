#include "meshUtils.h"

#include<deque>
#include<vector>

template<template<typename> typename tContainer>
auto visualizeNormals(meshUtils::Mesh<tContainer> const &aMesh) {
  meshUtils::Mesh<std::deque> result;
  
  auto sphere = meshUtils::makeUnitSphere<tContainer>(3, 1);
  for(auto const &face : aMesh) {
    auto normal = (face[1] - face[0]).cross(face[2] - face[0]).normalized();
    auto base = (face[0] + face[1] + face[2]) / 3.0f;
    auto average = ((face[1] - face[0]).norm() + (face[2] - face[1]).norm() + (face[0] - face[2]).norm()) / 3.0f;
    auto copy = sphere;
    auto shrink = meshUtils::Transform::Identity() * (average / 15.0f);
    meshUtils::transform(copy, shrink, base + normal * (average / 7.0f));
    std::copy(copy.cbegin(), copy.cend(), std::back_inserter(result));
  }
  return result;
}

void testDequeDivisor(char const * const aName, int32_t const aSectors, int32_t const aBelts, float const aRadius, int32_t const aDivisor) {
  std::string name{"test_"};
  name += aName;
  name += ".stl";

  auto sphere = meshUtils::makeUnitSphere<std::deque>(aSectors, aBelts);
  meshUtils::Transform inflate = meshUtils::Transform::Identity() * aRadius;
  meshUtils::transform(sphere, inflate);
  meshUtils::writeMesh(sphere, name);

  auto back = meshUtils::readMesh<std::deque>(name);
  meshUtils::transform(back, inflate);

  back = meshUtils::divideLargeTriangles(back, aDivisor);
  meshUtils::Vertex disp = meshUtils::Vertex::Zero();
  disp[0] = aRadius;
  meshUtils::transform(back, disp);
  auto nameBack = "back_" + name;
  meshUtils::standardizeVertices(back);
  meshUtils::standardizeNormals(back);
  meshUtils::writeMesh(back, nameBack);

  auto nameNorm = "norm_" + name;
  meshUtils::writeMesh(visualizeNormals(back), nameNorm);
}

void testVectorMax(char const * const aName, int32_t const aSectors, int32_t const aBelts, float const aRadius, float const aMaxSide) {
  std::string name{"test_"};
  name += aName;
  name += ".stl";

  auto sphere = meshUtils::makeUnitSphere<std::vector>(aSectors, aBelts);
  meshUtils::Transform inflate = meshUtils::Transform::Identity() * aRadius;
  meshUtils::transform(sphere, inflate);
  meshUtils::writeMesh(sphere, name);

  auto back = meshUtils::readMesh<std::vector>(name);
  meshUtils::transform(back, inflate);

  back = meshUtils::divideLargeTriangles(back, aMaxSide);
  meshUtils::Vertex disp = meshUtils::Vertex::Zero();
  disp[0] = aRadius;
  meshUtils::transform(back, disp);
  auto nameBack = "back_" + name;
  meshUtils::standardizeVertices(back);
  meshUtils::standardizeNormals(back);
  meshUtils::writeMesh(back, nameBack);

  auto nameNorm = "norm_" + name;
  meshUtils::writeMesh(visualizeNormals(back), nameNorm);
}

void testCustomStl(char * const aName) {
  auto mesh = meshUtils::readMesh<std::deque>(aName);

  meshUtils::standardizeVertices(mesh);
  meshUtils::standardizeNormals(mesh);

  std::string nameBack{"back_"};
  nameBack += aName;
  meshUtils::writeMesh(mesh, nameBack);

  std::string nameNorm{"norm_"};
  nameNorm += aName;
  meshUtils::writeMesh(visualizeNormals(mesh), nameNorm);
}

int main(int argc, char **argv) {
  testDequeDivisor("dequeDivisor", 7, 7, 3.0f, 3);
  testVectorMax("vectorMax", 3, 1, 13.0f, 11.0f);

  if(argc > 1) {
    testCustomStl(argv[1]);
  }
  else { // Nothing to do
  }
}
