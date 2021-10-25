#include "mesh.h"
#include "bezierMesh.h"

#include<deque>
#include<vector>

using Real = float;

auto visualizeNormals(Mesh<Real> const &aMesh) {
  Mesh<Real> result;
  Mesh<Real> sphere;
  sphere.makeUnitSphere(3, 1);
  for(auto &face : aMesh) {
    auto normal = getNormal(face).normalized();
    auto base = (face[0] + face[1] + face[2]) / 3.0f;
    auto average = ((face[1] - face[0]).norm() + (face[2] - face[1]).norm() + (face[0] - face[2]).norm()) / 3.0f;
    auto copy = sphere;
    auto shrink = Transform<Real>::Identity() * (average / 15.0f);
    copy.transform(shrink, base + normal * (average / 7.0f));
    std::copy(copy.cbegin(), copy.cend(), std::back_inserter(result));
  }
  return result;
}

auto visualizeVertexNormals(Mesh<Real> const &aMesh) { // Only for spheres
  Mesh<Real> result;
  Mesh<Real> sphere;
  sphere.makeUnitSphere(3, 1);
  auto &first = aMesh[0u];
  auto distance = ((first[1u] - first[0u]).norm() + (first[2u] - first[1u]).norm() + (first[0u] - first[2u]).norm()) / 3.0f;
  auto shrink = Transform<Real>::Identity() * (distance / 5.0f);
  sphere *= shrink;
  for(auto i : aMesh.getVertex2averageNormals()) {
    auto displacement = i.first + i.second * distance;
    auto copy = sphere;
    copy += displacement;
    std::copy(copy.cbegin(), copy.cend(), std::back_inserter(result));
  }
  return result;
}

void testDequeDivisor(char const * const aName, int32_t const aSectors, int32_t const aBelts, float const aRadius, int32_t const aDivisor) {
  std::string name{"test_"};
  name += aName;
  name += ".stl";

  Mesh<Real> sphere;
  sphere.makeUnitSphere(aSectors, aBelts);
  Transform<Real> inflate = Transform<Real>::Identity() * aRadius;
  sphere *= inflate;
  sphere.writeMesh(name);

  Mesh<Real> back;
  back.readMesh(name);
  back *= inflate;

  back.splitTriangles(aDivisor);
  Vector<Real> disp = Vector<Real>::Zero();
  disp[0] = aRadius;
  back += disp;
  auto nameBack = "back_" + name;
  back.standardizeVertices();
  back.standardizeNormals();
  back.writeMesh(nameBack);

  auto nameNorm = "norm_" + name;
  visualizeNormals(back).writeMesh(nameNorm);

  auto nameVertexNorm = "vertexNorm_" + name;
  visualizeVertexNormals(back).writeMesh(nameVertexNorm);
}

void testVectorMax(char const * const aName, int32_t const aSectors, int32_t const aBelts, float const aRadius, float const aMaxSide) {
  std::string name{"test_"};
  name += aName;
  name += ".stl";

  Mesh<Real> sphere;
  sphere.makeUnitSphere(aSectors, aBelts);
  Transform<Real> inflate = Transform<Real>::Identity() * aRadius;
  sphere *= inflate;
  sphere.writeMesh(name);

  Mesh<Real> back;
  back.readMesh(name);
  back *= inflate;

  back.splitTriangles(aMaxSide);
  auto nameBack = "back_" + name;
  back.standardizeVertices();
  back.standardizeNormals();
  back.writeMesh(nameBack);

  auto nameNorm = "norm_" + name;
  visualizeNormals(back).writeMesh(nameNorm);

  auto nameVertexNorm = "vertexNorm_" + name;
  visualizeVertexNormals(back).writeMesh(nameVertexNorm);
}

void testBarycentric2plane(char const * const aName, int32_t const aSectors, int32_t const aBelts, float const aRadius, int32_t const aDivisor) {
  std::string name{"baryorig_"};
  name += aName;
  name += ".stl";

  Mesh<Real> sphere;
  sphere.makeUnitSphere(aSectors, aBelts);
  Transform<Real> inflate = Transform<Real>::Identity() * aRadius;
  sphere *= inflate;
  sphere.standardizeVertices();
  sphere.standardizeNormals();
  sphere.writeMesh(name);

  BezierMesh<Real> bezier(sphere);
  auto planified = bezier.interpolate(aDivisor);

  name = "bary2plane_";
  name += aName;
  name += ".stl";
  planified.writeMesh(name);

  Mesh<Real> result;
  auto controlPoints = bezier.dumpControlPoints();
  sphere.makeUnitSphere(3, 1);
  inflate = Transform<Real>::Identity() * cgPi<Real> * aRadius / (aBelts + 1) / 20;
  sphere *= inflate;
  uint32_t i = 0u;
  for(auto const &controlPoint : controlPoints) {
    if(i < 5u) {
      auto copy = sphere;
      copy += controlPoint;
      std::copy(copy.cbegin(), copy.cend(), std::back_inserter(result));
    }
    else { // Nothing to do
    }
    i = (i + 1u) % 10u;
  }
  name = "baryControl_";
  name += aName;
  name += ".stl";
  result.writeMesh(name);
}

void testCustomStl(char * const aName) {
  Mesh<Real> mesh;
  mesh.readMesh(aName);

  mesh.standardizeVertices();
  mesh.standardizeNormals();

  std::string nameBack{"back_"};
  nameBack += aName;
  mesh.writeMesh(nameBack);

  std::string nameNorm{"norm_"};
  nameNorm += aName;
  visualizeNormals(mesh).writeMesh(nameNorm);
}

int main(int argc, char **argv) {
  testDequeDivisor("dequeDivisor", 7, 7, 3.0f, 3);
  testVectorMax("vectorMax", 3, 1, 13.0f, 11.0f);
  testVectorMax("vectorIdentity", 5, 2, 1.0f, 11.0f);
  testBarycentric2plane("7x3", 7, 3, 11.1f, 2);

  if(argc > 1) {
    testCustomStl(argv[1]);
  }
  else { // Nothing to do
  }
  return 0;
}
