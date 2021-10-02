#include "meshUtils.h"

#include<deque>
#include<vector>

template<template<typename> typename tContainer>
void test(char const * const aName, int32_t const aSectors, int32_t const aBelts, float aRadius) {
  std::string name{"test_"};
  name += aName;
  name += ".stl";

  auto sphere = meshUtils::makeUnitSphere<tContainer>(aSectors, aBelts);
  meshUtils::Transform inflate;
  inflate << aRadius, 0.0f, 0.0f, 0.0f, aRadius, 0.0f, 0.0f, 0.0f, aRadius;
  meshUtils::Vertex    leave;
  leave << 0.0f, 0.0f, 0.0f;
  meshUtils::transform(sphere, inflate, leave);
  meshUtils::writeMesh(sphere, name);

  auto back = meshUtils::readMesh<tContainer>(name, inflate, leave);
  back = meshUtils::divideLargeTriangles(back, 11.1f);
  inflate /= aRadius;
  leave[0] = aRadius;
  meshUtils::transform(back, inflate, leave);
  name = "back_" + name;
  meshUtils::standardizeVertices(back);
  meshUtils::writeMesh(back, name);
}

int main() {
  test<std::deque>("deque", 3, 1, 1.0f);
  test<std::vector>("vector", 11, 5, 11.1f);
}
