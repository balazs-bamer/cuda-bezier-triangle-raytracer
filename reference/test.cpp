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

void testBarycentric2plane(char const * const aName, int32_t const aSectors, int32_t const aBelts, float const aRadius, int32_t const aDivisor) {
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

/*  auto parameters = bezier.getStuff();
  name = "baryPara_";
  name += aName;
  name += ".stl";
  parameters.writeMesh(name);*/
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
  testDequeDivisor("dequeDivisor", 7, 7, 3.0f, 3);
  testVectorMax("vectorMax", 4, 2, 13.0f, 11.0f);
  testVectorMax("vectorIdentity", 4, 1, 1.0f, 11.0f);
  testBarycentric2plane("4x2", 4, 2, 3.0f, 4);
  testBarycentric2plane("7x5", 7, 5, 3.0f, 4);


 /* Mesh<Real> sphere;
  sphere.makeUnitSphere(4, 2);
  sphere.standardizeVertices();
  sphere.standardizeNormals();
  auto neighbours = sphere.getFace2neighbours();
  for(uint32_t i = 0u; i < neighbours.size(); ++i) {
    auto const &neigh = neighbours[i];
    auto const &triangle = sphere[i];
    std::cout << i << '\n';
    for(uint32_t j = 0u; j < 3u; ++j) {
      std::cout << "  " << std::setw(3) << neigh.mFellowTriangles[j] << ":   " << (int)(neigh.mFellowCommonSideStarts[j]);
      std::cout << ' ' << std::setw(10) << std::setprecision(3) << (::abs(triangle[j](0)) < 0.001 ? 0.0f : triangle[j](0));
      std::cout << ' ' << std::setw(10) << std::setprecision(3) << (::abs(triangle[j](1)) < 0.001 ? 0.0f : triangle[j](1));
      std::cout << ' ' << std::setw(10) << std::setprecision(3) << (::abs(triangle[j](2)) < 0.001 ? 0.0f : triangle[j](2));
      std::cout << " -=-";

      auto fellow = sphere[neigh.mFellowTriangles[j]];
      auto fellowCommonStart = (int)(neigh.mFellowCommonSideStarts[j]);
      Vector<Real> ownSide = triangle[j] - triangle[(j+1)%3];
      Vector<Real> fellowSide = fellow[fellowCommonStart] - fellow[(fellowCommonStart+1)%3];
      for(int k = 0; k < 3; ++k) {
        if(::abs(ownSide(k)) < 0.001) ownSide(k) = 0.0;
        if(::abs(fellowSide(k)) < 0.001) fellowSide(k) = 0.0;
      }
      std::cout << ' ' << std::setw(10) << std::setprecision(3) << ownSide(0);
      std::cout << ' ' << std::setw(10) << std::setprecision(3) << ownSide(1);
      std::cout << ' ' << std::setw(10) << std::setprecision(3) << ownSide(2);
      std::cout << " -";
      std::cout << ' ' << std::setw(10) << std::setprecision(3) << fellowSide(0);
      std::cout << ' ' << std::setw(10) << std::setprecision(3) << fellowSide(1);
      std::cout << ' ' << std::setw(10) << std::setprecision(3) << fellowSide(2);
      std::cout << '\n';
      if((ownSide - fellowSide).norm() > 0.01f && (ownSide + fellowSide).norm() > 0.01f) throw "own fellow sides mismatch";
    }
  }


  for(uint32_t indexFace = 0u; indexFace < sphere.size(); ++indexFace) {
    auto const &neigh = neighbours[indexFace];
    auto const &originalTriangle = sphere[indexFace];
    auto normal = getNormal(originalTriangle);
    std::cout << std::setw(14) << normal(0) << std::setw(14) << normal(1) << std::setw(14) << normal(2) << "   " << indexFace << "\n";
    for(uint32_t indexVertex = 0u; indexVertex < 3u; ++indexVertex) {
      auto other = getNormal(sphere[neigh.mFellowTriangles[indexVertex]]);
      std::cout << std::setw(14) << other(0) << std::setw(14) << other(1) << std::setw(14) << other(2) << "   " << indexFace << ' ' << indexVertex << '\n';
    }
  }*/

/*  BezierMesh<Real> bezier(sphere);
  for(uint32_t i = 0u; i < bezier.size(); ++i) {
    auto const neighbours = bezier[i].getNeighbours();
    std::cout << i << '\n';
    for(uint32_t j = 0u; j < 3u; ++j) {
      std::cout << "  " << std::setw(3) << neighbours[j] << '\n';
    }
  }*/

  if(argc > 1) {
    testCustomStl(argv[1], 4);
  }
  else { // Nothing to do
  }
  return 0;
}
