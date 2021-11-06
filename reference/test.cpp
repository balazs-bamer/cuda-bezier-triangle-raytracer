#include "mesh.h"
#include "bezierMesh.h"

#include<deque>
#include<vector>
#include<iomanip>
#include<iostream>

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

void testBezier2plane(char const * const aName, int32_t const aSectors, int32_t const aBelts, float const aRadius, int32_t const aDivisor) {
  Mesh<Real> sphere;
  sphere.makeUnitSphere(aSectors, aBelts);
  Transform<Real> inflate = Transform<Real>::Identity() * aRadius;
  sphere *= inflate;
  sphere.standardizeVertices();
  sphere.standardizeNormals();

  std::string name{"baryOrig_"};
  name += aName;
  name += ".stl";
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
    if(i < 12u || i >= 3 || i >= 8u) {
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

void testBezierSplitTall(char const * const aName, int32_t const aSectors, int32_t const aBelts, Vector<Real> const &aSize, int32_t const aDivisor) {
  Mesh<Real> ellipsoid;
  ellipsoid.makeEllipsoid(aSectors, aBelts, aSize);
  ellipsoid.standardizeVertices();
  ellipsoid.standardizeNormals();

  std::string name{"barySplitOrig_"};
  name += aName;
  name += ".stl";
  ellipsoid.writeMesh(name);

  name = "barySplitVertexNorm_";
  name += aName;
  name += ".stl";
  visualizeVertexNormals(ellipsoid).writeMesh(name);

  BezierMesh<Real> bezier0(ellipsoid);
  Mesh<Real> split1 = bezier0.splitThickBezierTriangles();
  split1.standardizeVertices();
  split1.standardizeNormals();

  name = "barySplit1_";
  name += aName;
  name += ".stl";
  split1.writeMesh(name);

  BezierMesh<Real> bezier1(split1);
  Mesh<Real> split2 = bezier1.splitThickBezierTriangles();

  name = "barySplit2_";
  name += aName;
  name += ".stl";
  split2.writeMesh(name);
}

void measureApproximation(uint32_t const aSplitSteps, int32_t const aSectors, int32_t const aBelts, Vector<Real> const &aSize, int32_t const aDivisor) {
  Mesh<Real> ellipsoid;
  ellipsoid.makeEllipsoid(aSectors, aBelts, aSize);
  ellipsoid.standardizeVertices();
  ellipsoid.standardizeNormals();

  for(int32_t i = 0u; i < aSplitSteps; ++i) {
    BezierMesh<Real> bezier(ellipsoid);
    ellipsoid = bezier.splitThickBezierTriangles();
    ellipsoid.standardizeVertices();
    ellipsoid.standardizeNormals();
  }

  BezierMesh<Real> bezier(ellipsoid);
  auto planified = bezier.interpolate(aDivisor);
  planified.standardizeVertices();
  auto vertices = planified.getVertices();  // Due to Clough-Tocher division of small meshes, this will introduce big errors
                                            // Results are only meaningful for sufficiently fine initial meshes, which are easy
                                            // to follow with cubic Bezier triangles.
  Real sum = 0.0f;
  for(auto const &vertex : vertices) {
    Spherical spherical(vertex(0) / aSize(0), vertex(1) / aSize(1), vertex(2) / aSize(2));
    Vertex<Real> ethalon { aSize(0) * sin(spherical.mLatitude) * cos(spherical.mLongitude),
                           aSize(1) * sin(spherical.mLatitude) * sin(spherical.mLongitude),
                           aSize(2) * cos(spherical.mLatitude) };
    sum += (vertex - ethalon).squaredNorm() / ethalon.squaredNorm();
  }
  auto error = sum / vertices.size();
  std::cout << "SplitSteps: " << aSplitSteps << " Sectors: " << aSectors << " Belts: " << aBelts <<
               " Size: " << aSize(0) << ' ' << aSize(1) << ' ' << aSize(2) << " Divisor: " << aDivisor <<
               " error: " << sum << '\n';
}

void testCustomStl(char * const aName, int32_t const aDivisor) {  // TODO this does not work perfectly for complex and extreme surfaces like robot.stl
  Mesh<Real> mesh;
  mesh.readMesh(aName);

  mesh.standardizeVertices();
  mesh.standardizeNormals();

  std::string name{"back_"};
  name += aName;
  mesh.writeMesh(name);

  name = "norm_";
  name += aName;
  visualizeNormals(mesh).writeMesh(name);

  BezierMesh<Real> bezier(mesh);
  auto planified = bezier.interpolate(aDivisor);

  name = "bary2plane_";
  name += aName;
  planified.writeMesh(name);
}

int main(int argc, char **argv) {
  /*testDequeDivisor("dequeDivisor", 7, 7, 3.0f, 3);

  testVectorMax("vectorMax", 4, 2, 13.0f, 11.0f);
  testVectorMax("vectorIdentity", 4, 1, 1.0f, 11.0f);

  testBezier2plane("4x2", 4, 2, 3.0f, 4);
  testBezier2plane("7x5", 7, 5, 3.0f, 4);*/

  Vector<Real> ellipsoidAxes(1.0f, 4.0f, 2.0f);

/*  testBezierSplitTall("7x3", 7, 3, ellipsoidAxes, 1);
  testBezierSplitTall("15x5", 15, 5, ellipsoidAxes, 1);*/

  measureApproximation(0, 4, 1, ellipsoidAxes, 1);
  measureApproximation(0, 7, 3, ellipsoidAxes, 3);
  measureApproximation(0, 15, 5, ellipsoidAxes, 3);
  measureApproximation(1, 7, 3, ellipsoidAxes, 3);
  measureApproximation(1, 15, 5, ellipsoidAxes, 3);

  if(argc > 1) {
    testCustomStl(argv[1], 4);
  }
  else { // Nothing to do
  }
  return 0;
}
