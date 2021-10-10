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

  meshUtils::Vertex disp = meshUtils::Vertex::Zero();
  auto back = meshUtils::readMesh<std::deque>(name, inflate, disp);

  back = meshUtils::divideLargeTriangles(back, aDivisor);
  disp[0] = aRadius;
  meshUtils::transform(back, disp);
  auto nameBack = "back_" + name;
  meshUtils::standardizeVertices(back);
  auto face2neighbour = meshUtils::standardizeNormals(back);
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

  meshUtils::Vertex disp = meshUtils::Vertex::Zero();
  auto back = meshUtils::readMesh<std::vector>(name, inflate, disp);

  back = meshUtils::divideLargeTriangles(back, aMaxSide);
  disp[0] = aRadius;
  meshUtils::transform(back, disp);
  auto nameBack = "back_" + name;
  meshUtils::standardizeVertices(back);
  auto face2neighbour = meshUtils::standardizeNormals(back);
  meshUtils::writeMesh(back, nameBack);

  auto nameNorm = "norm_" + name;
  meshUtils::writeMesh(visualizeNormals(back), nameNorm);
}

int main() {
  testDequeDivisor("dequeDivisor", 3, 1, 3.0f, 3);
  testVectorMax("vectorMax", 7, 7, 23.0f, 1.0f);
}
