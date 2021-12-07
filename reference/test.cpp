#include "mesh.h"
#include "bezierMesh.h"

#include<deque>
#include<vector>
#include<iomanip>
#include<iostream>

using Real = float;

std::string const cgBaseDir("output/");

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

auto visualizeRay(Ray<Real> const aRay, Real const aLength, Real const aRadius) {
  Mesh<Real> result;
  Vector<Real> const perpendicular0 = getAperpendicular(aRay.mDirection);
  Vector<Real> const perpendicular1 = aRay.mDirection.cross(perpendicular0);
  Vector<Real> const rib0 = aRadius * perpendicular0;  // p0 * cos0 + p1 * sin0
  Vector<Real> const rib1 = aRadius * (perpendicular0 * ::cos(2.0f / 3.0f * cgPi<Real>) + perpendicular1 * ::sin(2.0f / 3.0f * cgPi<Real>));
  Vector<Real> const rib2 = aRadius * (perpendicular0 * ::cos(4.0f / 3.0f * cgPi<Real>) + perpendicular1 * ::sin(4.0f / 3.0f * cgPi<Real>));
  Vector<Real> const end = aRay.mStart + aRay.mDirection * aLength;
  result.push_back({aRay.mStart + rib0, aRay.mStart + rib1, aRay.mStart + rib2});
  result.push_back({end + rib0, end + rib1, end + rib2});
  int32_t n = static_cast<int32_t>(::ceil(::abs(aLength) * 0.1f / aRadius));
  Real delta = aLength / n;
  for(int i = 0; i < n; ++i) {
    Vector<Real> section0 = aRay.mStart + aRay.mDirection * delta * i;
    Vector<Real> section1 = aRay.mStart + aRay.mDirection * delta * (i + 1);
    result.push_back({section0 + rib0, section1 + rib0, section0 + rib1});
    result.push_back({section0 + rib1, section1 + rib0, section1 + rib1});
    result.push_back({section0 + rib1, section1 + rib1, section0 + rib2});
    result.push_back({section0 + rib2, section1 + rib1, section1 + rib2});
    result.push_back({section0 + rib2, section1 + rib2, section0 + rib0});
    result.push_back({section0 + rib0, section1 + rib2, section1 + rib0});
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
  sphere.writeMesh(cgBaseDir + name);

  Mesh<Real> back;
  back.readMesh(cgBaseDir + name);
  back *= inflate;

  back.splitTriangles(aDivisor);
  Vector<Real> disp = Vector<Real>::Zero();
  disp[0] = aRadius;
  back += disp;
  auto nameBack = "back_" + name;
  back.standardizeVertices();
  back.standardizeNormals();
  back.writeMesh(cgBaseDir + nameBack);

  auto nameNorm = "norm_" + name;
  visualizeNormals(back).writeMesh(cgBaseDir + nameNorm);

  auto nameVertexNorm = "vertexNorm_" + name;
  visualizeVertexNormals(back).writeMesh(cgBaseDir + nameVertexNorm);
}

void testVectorMax(char const * const aName, int32_t const aSectors, int32_t const aBelts, float const aRadius, float const aMaxSide) {
  std::string name{"test_"};
  name += aName;
  name += ".stl";

  Mesh<Real> sphere;
  sphere.makeUnitSphere(aSectors, aBelts);
  Transform<Real> inflate = Transform<Real>::Identity() * aRadius;
  sphere *= inflate;
  sphere.writeMesh(cgBaseDir + name);

  Mesh<Real> back;
  back.readMesh(cgBaseDir + name);
  back *= inflate;

  back.splitTriangles(aMaxSide);
  auto nameBack = "back_" + name;
  back.standardizeVertices();
  back.standardizeNormals();
  back.writeMesh(cgBaseDir + nameBack);

  auto nameNorm = "norm_" + name;
  visualizeNormals(back).writeMesh(cgBaseDir + nameNorm);

  auto nameVertexNorm = "vertexNorm_" + name;
  visualizeVertexNormals(back).writeMesh(cgBaseDir + nameVertexNorm);
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
  sphere.writeMesh(cgBaseDir + name);

  BezierMesh<Real> bezier(sphere);
  auto planified = bezier.interpolate(aDivisor);

  name = "bary2plane_";
  name += aName;
  name += ".stl";
  planified.writeMesh(cgBaseDir + name);

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
  result.writeMesh(cgBaseDir + name);
}

void testBezierSplitTall(char const * const aName, int32_t const aSectors, int32_t const aBelts, Vector<Real> const &aSize, int32_t const aDivisor) {
  Mesh<Real> ellipsoid;
  ellipsoid.makeEllipsoid(aSectors, aBelts, aSize);
  ellipsoid.standardizeVertices();
  ellipsoid.standardizeNormals();

  std::string name{"barySplitOrig_"};
  name += aName;
  name += ".stl";
  ellipsoid.writeMesh(cgBaseDir + name);

  name = "barySplitVertexNorm_";
  name += aName;
  name += ".stl";
  visualizeVertexNormals(ellipsoid).writeMesh(cgBaseDir + name);

  BezierMesh<Real> bezier0(ellipsoid);
  Mesh<Real> split1 = bezier0.splitThickBezierTriangles();
  split1.standardizeVertices();
  split1.standardizeNormals();

  name = "barySplit1_";
  name += aName;
  name += ".stl";
  split1.writeMesh(cgBaseDir + name);

  BezierMesh<Real> bezier1(split1);
  Mesh<Real> split2 = bezier1.splitThickBezierTriangles();

  name = "barySplit2_";
  name += aName;
  name += ".stl";
  split2.writeMesh(cgBaseDir + name);
}

void dump(BezierIntersection<Real> const aIntersection) {
  std::cout << "result: "
            << std::setw(11) << std::setprecision(4) << aIntersection.mIntersection.mPoint(0)
            << std::setw(11) << std::setprecision(4) << aIntersection.mIntersection.mPoint(1)
            << std::setw(11) << std::setprecision(4) << aIntersection.mIntersection.mPoint(2) << " dist: "
            << std::setw(11) << std::setprecision(4) << aIntersection.mIntersection.mDistance << " norm: "
            << std::setw(11) << std::setprecision(4) << aIntersection.mNormal(0)
            << std::setw(11) << std::setprecision(4) << aIntersection.mNormal(1)
            << std::setw(11) << std::setprecision(4) << aIntersection.mNormal(2)              << " what: "
            << std::setw(11) << std::setprecision(4) << static_cast<unsigned>(aIntersection.mWhat) << "\n\n\n";
}

void testBezierIntersection(char const * const aName, int32_t const aSectors, int32_t const aBelts, Vector<Real> const &aSize, Vector<Real> const &aDirection) {
  std::cout << aName << '\n';
  Mesh<Real> ellipsoid;
  ellipsoid.makeEllipsoid(aSectors, aBelts, aSize);

  Mesh<Real> bullet;
  bullet.makeUnitSphere(5, 3);
  bullet *= 0.05f;
  Mesh<Real> intersections;
  Mesh<Real> objects;
  Mesh<Real> raws;

  Ray<Real> ray(Vertex<Real>{0.0f, 0.0f, 0.0f}, aDirection);
  Vector<Real> displacement{5.0f, 0.0f, 0.0f};

  auto original = ray;
  std::vector<Vertex<Real>> points;
  for(;;) {
    ellipsoid += displacement;                    // We start from outside.
    ellipsoid.standardizeVertices();
    ellipsoid.standardizeNormals();
    std::copy(ellipsoid.cbegin(), ellipsoid.cend(), std::back_inserter(raws));
    BezierMesh<Real> bezier(ellipsoid);
    auto object = bezier.interpolate(5);
    std::copy(object.cbegin(), object.cend(), std::back_inserter(objects));

    auto intersection = bezier.intersect(ray);
    dump(intersection);
    if(intersection.mWhat != BezierIntersection<Real>::What::cIntersect) {
      break;
    }
    else { // Nothing to do
    }
std::cout << "now should be inside, ray length " << intersection.mIntersection.mDistance << "\n";
std::cout << "-------------------------------------------------------------------\n\n";
    auto copy = bullet;
    copy += intersection.mIntersection.mPoint;
    std::copy(copy.cbegin(), copy.cend(), std::back_inserter(intersections));
    copy = bullet;
    copy += intersection.mIntersection.mPoint + intersection.mNormal * 1.0f;
    std::copy(copy.cbegin(), copy.cend(), std::back_inserter(intersections));
    points.push_back(intersection.mIntersection.mPoint);
    ray.mStart = intersection.mIntersection.mPoint;

    intersection = bezier.intersect(ray);         // TODO why doesn't the ray come out of the last shape? It wants a neighbouring triangle, and in turn a neighbouring one.
    dump(intersection);                           // Double doesn't help.
    if(intersection.mWhat != BezierIntersection<Real>::What::cIntersect) {
      break;
    }
    else { // Nothing to do
    }
std::cout << "now should be outside, ray length " << intersection.mIntersection.mDistance << "\n";
std::cout << "-------------------------------------------------------------------\n\n";
    copy = bullet;
    copy += intersection.mIntersection.mPoint;
    std::copy(copy.cbegin(), copy.cend(), std::back_inserter(intersections));
    copy = bullet;
    copy += intersection.mIntersection.mPoint + intersection.mNormal * 1.0f;
    std::copy(copy.cbegin(), copy.cend(), std::back_inserter(intersections));
    points.push_back(intersection.mIntersection.mPoint);
    ray.mStart = intersection.mIntersection.mPoint;
  }
  ellipsoid += displacement;                    // We start from outside.
  ellipsoid.standardizeVertices();
  ellipsoid.standardizeNormals();
  std::copy(ellipsoid.cbegin(), ellipsoid.cend(), std::back_inserter(raws));
  BezierMesh<Real> bezier(ellipsoid);
  auto object = bezier.interpolate(5);
  std::copy(object.cbegin(), object.cend(), std::back_inserter(objects));
  points.push_back(points.back() + ray.mDirection * 11.0f);

  std::cout << "error: " << original.getAverageErrorSquared(points) << '\n';

  std::string name{"intersectionObject_"};
  name += aName;
  name += ".stl";
  objects.writeMesh(cgBaseDir + name);

  name ="intersectionRaw_";
  name += aName;
  name += ".stl";
  raws.writeMesh(cgBaseDir + name);

  name = "intersectionLocation_";
  name += aName;
  name += ".stl";
  intersections.writeMesh(cgBaseDir + name);

  auto beam = visualizeRay(original, (points.back() - original.mStart).norm(), 0.02f);

  name = "intersectionRay_";
  name += aName;
  name += ".stl";
  beam.writeMesh(cgBaseDir + name);

  std::cout << '\n';
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
    Vertex<Real> ethalon { aSize(0) * sin(spherical.mInclination) * cos(spherical.mAzimuth),
                           aSize(1) * sin(spherical.mInclination) * sin(spherical.mAzimuth),
                           aSize(2) * cos(spherical.mInclination) };
    sum += (vertex - ethalon).squaredNorm() / ethalon.squaredNorm();
  }
  auto error = sum / vertices.size();
  std::cout << "SplitSteps: " << aSplitSteps << " Sectors: " << std::setw(2) << aSectors << " Belts: " << std::setw(2) << aBelts <<
               " Size: " << aSize(0) << ' ' << aSize(1) << ' ' << aSize(2) << " Divisor: " << aDivisor <<
               " error: " << std::setw(14) << std::setprecision(8) << error << '\n';
}

std::deque<BezierTriangle<float>> gFollowers;

void visualizeFollowers(char const * const aName) {
  for(uint32_t i = 0u; i < 1; ++i) {
    gFollowers.pop_front();            // Remove uninteresting items
  }
  for(auto const &bezier : gFollowers) {
    std::cout << "what out\n";
  }
}

void testCustomStl(char * const aName, int32_t const aDivisor) {  // TODO this does not work perfectly for complex and extreme surfaces like robot.stl
  Mesh<Real> mesh;
  mesh.readMesh(cgBaseDir + aName);

  mesh.standardizeVertices();
  mesh.standardizeNormals();

  std::string name{"back_"};
  name += aName;
  mesh.writeMesh(cgBaseDir + name);

  name = "norm_";
  name += aName;
  visualizeNormals(mesh).writeMesh(cgBaseDir + name);

  BezierMesh<Real> bezier(mesh);
  auto planified = bezier.interpolate(aDivisor);

  name = "bary2plane_";
  name += aName;
  planified.writeMesh(cgBaseDir + name);
}

int main(int argc, char **argv) {
  Vector<Real> ellipsoidAxes(1.0f, 4.0f, 2.0f);

/*  testDequeDivisor("dequeDivisor", 7, 7, 3.0f, 3);

  testVectorMax("vectorMax", 4, 2, 13.0f, 11.0f);
  testVectorMax("vectorIdentity", 4, 1, 1.0f, 11.0f);

  testBezier2plane("4x2", 4, 2, 3.0f, 4);
  testBezier2plane("7x5", 7, 5, 3.0f, 4);

  testBezierSplitTall("7x3", 7, 3, ellipsoidAxes, 1);
  testBezierSplitTall("15x5", 15, 5, ellipsoidAxes, 1);*/

/*  testBezierIntersection("7x3a", 7, 3, ellipsoidAxes, {1.0f, 0.05f, 0.02f});
  testBezierIntersection("7x3b", 7, 3, ellipsoidAxes, {1.0f, 0.05f, -0.022f});
  testBezierIntersection("9x4a", 9, 4, ellipsoidAxes, {1.0f, -0.03f, 0.035f});*/
  testBezierIntersection("9x4b", 9, 4, ellipsoidAxes, {1.0f, -0.03f, -0.045f});
 /* visualizeFollowers("follow");*/

/*  measureApproximation(0, 4, 1, ellipsoidAxes, 1);     // SplitSteps: 0 Sectors:  4 Belts:  1 Size: 1 4 2 Divisor: 1 error:      1.2555894
  measureApproximation(0, 7, 3, ellipsoidAxes, 3);       // SplitSteps: 0 Sectors:  7 Belts:  3 Size: 1 4 2 Divisor: 3 error:   0.0022721614
  measureApproximation(0, 15, 5, ellipsoidAxes, 3);      // SplitSteps: 0 Sectors: 15 Belts:  5 Size: 1 4 2 Divisor: 3 error:  1.9426199e-05
  measureApproximation(1, 7, 3, ellipsoidAxes, 3);       // SplitSteps: 1 Sectors:  7 Belts:  3 Size: 1 4 2 Divisor: 3 error:  0.00070956006
  measureApproximation(1, 15, 5, ellipsoidAxes, 3);      // SplitSteps: 1 Sectors: 15 Belts:  5 Size: 1 4 2 Divisor: 3 error:  0.00040229771
  measureApproximation(2, 7, 3, ellipsoidAxes, 3);       // SplitSteps: 2 Sectors:  7 Belts:  3 Size: 1 4 2 Divisor: 3 error:   0.0011259826
  measureApproximation(2, 15, 5, ellipsoidAxes, 3);      // SplitSteps: 2 Sectors: 15 Belts:  5 Size: 1 4 2 Divisor: 3 error:  6.7134395e-05
*/

  if(argc > 1) {
    testCustomStl(argv[1], 4);
  }
  else { // Nothing to do
  }
  return 0;
}
