#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_MESHUTILS
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_MESHUTILS

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t

#include "stl_reader.h"
#include <array>
#include <deque>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <Eigen/Dense>

#include <iostream>

namespace meshUtils {

using Vertex = Eigen::Vector3f;
using Transform = Eigen::Matrix3f;
using Triangle = std::array<Vertex, 3u>;

template<template<typename> typename tContainer>
using Mesh = tContainer<Triangle>;

using MeshVector = Mesh<std::vector>;

constexpr float cgPi = 3.1415926539f;
constexpr float cgStandardizeVerticesEpsilonFactor = 0.2f;

// TODO perhaps a function to split triangles having vertex on an edge. Probably not needed.

template<template<typename> typename tContainer>
void standardizeVertices(Mesh<tContainer> &aMesh);

template<template<typename> typename tContainer>
void standardizeNormals(Mesh<tContainer> &aMesh);

template<template<typename> typename tContainer>
void transform(Mesh<tContainer> &aMesh, Transform const &aTransform, Vertex const aDisplacement);

template<template<typename> typename tContainer>
auto divideLargeTriangles(Mesh<tContainer> &aMesh, float aMaxTriangleSide);

template<template<typename> typename tContainer>
auto readMesh(std::string const &aFilename, Transform const &aTransform, Vertex const aDisplacement);

template<template<typename> typename tContainer>
void writeMesh(Mesh<tContainer> const &aMesh, std::string const& aFilename);

template<template<typename> typename tContainer>
auto makeUnitSphere(int32_t const aSectors, int32_t const aBelts);

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<template<typename> typename tContainer>
void standardizeVertices(Mesh<tContainer> &aMesh) {        // Slow but was quick to implement. Moves all distinct but too close vertices to one place.
  float smallestSide = std::numeric_limits<float>::max();
  for(auto &triangle : aMesh) {
    for(uint32_t i = 0u; i < 3u; ++i) {
      smallestSide = std::min(smallestSide, (triangle[i] - triangle[(i + 1u) % 3u]).norm());
    }
  }
  auto epsilon = smallestSide * cgStandardizeVerticesEpsilonFactor;
  epsilon *= epsilon; // Spare sqrt
  for(uint32_t i = 0u; i < aMesh.size(); ++i) {
    for(uint32_t j = 0u; j < 3u; ++j) {
      for(uint32_t k = 0u; k < aMesh.size(); ++k) {
        for(uint32_t l = 0u; l < 3u; ++l) {
          if(i != k || j != l) {
            auto &v1 = aMesh[i][j];
            auto &v2 = aMesh[k][l];
            if((v1-v2).squaredNorm() < epsilon && (v1[0] < v2[0] || (v1[0] == v2[0] && v1[1] < v2[1]) || (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] < v2[2]))) {
              v1[0] = v2[0];
            }
            else { // nothing to do
            }
          }
          else { // nothing to do
          }
        }
      }
    }
  }
}

template<template<typename> typename tContainer>          // Should run after standardizeVertices
void standardizeNormals(Mesh<tContainer> &aMesh) {        // Makes all normalvectors (V[1]-V[0])x(V[2]-V[0]) point outwards.
  std::unordered_multimap<std::pair<uint32_t, uint32_t>, uint32_t> edge2face;
  std::unordered_multimap<uint32_t, uint32_t> vertex2face;
  std::unordered_map<Vertex, uint32_t> was;  // For presence testing.
  std::deque<Vertex> vertices;              // For indexing.
  float smallestX = std::numeric_limits<float>::max();
  uint32_t smallestIndexInAll;
  uint32_t indexFace = 0u;
  for(auto &face : aMesh) {
    uint32_t indexInFace = 0u;
    std::array<uint32_t, 3u> faceVertexIndices;
    for(auto &vertex : face) {
      smallestX = std::min(smallestX, vertex[0]);
      auto found = was.find(vertex);
      uint32_t indexInAll;
      if(found == was.end()) {
        indexInAll = vertices.size();
        faceVertexIndices[indexInFace] = indexInAll;
        was.insert(std::pair(vertex, indexInAll));
        vertices.push_back(vertex);
      }
      else {
        indexInAll = found->second;
        faceVertexIndices[indexInFace] = indexInAll;
      }
      if(vertex[0] < smallestX) {
        smallestX = vertex[0];
        smallestIndexInAll = indexInAll;
      }
      else { // Nothing to do
      }
      ++indexInFace;
    }
    for(uint32_t i = 0u; i < 3u; ++i) {
      vertex2face.insert(std::pair(faceVertexIndices[i], indexFace));
      auto low = faceVertexIndices[i];
      auto high = faceVertexIndices[(i + 1u) % 3u];
      if(low > high) {
        std::swap(low, high);
      }
      else { // Nothing to do
      }
      edge2face.insert(std::pair(std::pair(low, high), indexFace));
    }
    ++indexFace;
  }
  
}

template<template<typename> typename tContainer>
void transform(Mesh<tContainer> &aMesh, Transform const &aTransform, Vertex const aDisplacement) {
  for(auto &triangle : aMesh) {
    for(auto &vertex : triangle) {
      vertex = aTransform * vertex + aDisplacement;      
    }
  }
}

template<template<typename> typename tContainer>
auto divideLargeTriangles(Mesh<tContainer> &aMesh, float aMaxTriangleSide) {
  Mesh<tContainer> result;
  for(auto const &triangle : aMesh) {
    float maxSide = (triangle[0] - triangle[1]).norm();
    maxSide = std::max(maxSide, (triangle[0] - triangle[2]).norm());
    maxSide = std::max(maxSide, (triangle[1] - triangle[2]).norm());
    int32_t divisor = static_cast<int32_t>(std::ceil(maxSide / aMaxTriangleSide));
    auto vector01 = (triangle[1] - triangle[0]) / divisor;
    auto vector02 = (triangle[2] - triangle[0]) / divisor;
    auto lineBase = triangle[0];
    auto base0 = lineBase;
    auto base1 = (divisor > 1) ? (base0 + vector01) : triangle[1];
    auto base2 = (divisor > 1) ? (base0 + vector02) : triangle[2];
    for(int32_t i = 0; i < divisor - 1; ++i) {
      for(int32_t j = 0; j < divisor - i - 1; j++) {
        result.push_back({base0, base1, base2});
        auto base1next = base1 + vector02;
        result.push_back({base1, base1next, base2});
        base1 = base1next;
        base0 = base2;
        base2 += vector02;
      }
      result.push_back({base0, base1, base2});
      lineBase += vector01;
      base0 = lineBase;
      base1 = base0 + vector01;
      base2 = base0 + vector02;
    }
    result.push_back({base0, triangle[1], base2});
  }
  return result;
}

template<template<typename> typename tContainer>
auto readMesh(std::string const &aFilename, Transform const &aTransform, Vertex const aDisplacement) {
  Mesh<tContainer> work;
  stl_reader::StlMesh<float, int32_t> mesh(aFilename);
  for(int32_t indexTriangle = 0; indexTriangle < mesh.num_tris(); ++indexTriangle) {
    Triangle triangle;
    for(int32_t indexCorner = 0; indexCorner < 3; ++indexCorner) {
      float const * const coords = mesh.tri_corner_coords(indexTriangle, indexCorner);
      Eigen::Vector3f in;
      for(int32_t i = 0; i < 3; ++i) {
        in(i) = coords[i];
      }
      triangle[indexCorner] = aTransform * in + aDisplacement;
    }
    work.push_back(triangle);
  }
  return work;
}

template<template<typename> typename tContainer>
void writeMesh(Mesh<tContainer> const &aMesh, std::string const& aFilename) {
  std::ofstream out(aFilename);
  out << "solid Exported from Blender-2.82 (sub 7)\n";
  for(auto const & triangle : aMesh) {
    out << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";
    for(auto const & vertex : triangle) {
      out << "vertex " << vertex(0) << ' ' << vertex(1) << ' ' << vertex(2) << '\n';
    }
    out << "endloop\nendfacet\n";
  }
  out << "endsolid Exported from Blender-2.82 (sub 7)\n";
}

template<template<typename> typename tContainer>
auto makeUnitSphere(int32_t const aSectors, int32_t const aBelts) {
  Mesh<tContainer> result;
  float sectorAngleHalf = cgPi / aSectors;
  float sectorAngleFull = sectorAngleHalf * 2.0f;
  float beltAngle       = cgPi / (aBelts + 1.0f);
  float bias = 0.0f;
  float beltAngleUp = 0.0f;
  float beltAngleMiddle = beltAngle;
  float beltAngleDown = 2.0f * beltAngle;
  float beltRadiusUp = 0.0f;
  float beltRadiusMiddle = std::sin(beltAngleMiddle);
  float beltRadiusDown = std::sin(beltAngleDown);
  float beltZup = 1.0f;
  float beltZmiddle = std::cos(beltAngleMiddle);
  float beltZdown = std::cos(beltAngleDown);
  for(int32_t belt = 0; belt < aBelts; ++belt) {
    float sectorAngleUpDown = bias + sectorAngleHalf;
    float sectorAngleMiddle1 = bias + 0.0f;
    float sectorAngleMiddle2 = bias + sectorAngleFull;
    for(int32_t sector = 0; sector < aSectors; ++sector) {
      Vertex corner1(beltRadiusUp * std::sin(sectorAngleUpDown), beltRadiusUp * std::cos(sectorAngleUpDown), beltZup);
      Vertex corner2(beltRadiusMiddle * std::sin(sectorAngleMiddle1), beltRadiusMiddle * std::cos(sectorAngleMiddle1), beltZmiddle);
      Vertex corner3(beltRadiusMiddle * std::sin(sectorAngleMiddle2), beltRadiusMiddle * std::cos(sectorAngleMiddle2), beltZmiddle);
      result.push_back({corner1, corner2, corner3});
      corner1 = {beltRadiusDown * std::sin(sectorAngleUpDown), beltRadiusDown * std::cos(sectorAngleUpDown), beltZdown};
      result.push_back({corner2, corner3, corner1});
      sectorAngleUpDown += sectorAngleFull;
      sectorAngleMiddle1 = sectorAngleMiddle2;
      sectorAngleMiddle2 += sectorAngleFull;
    }
    beltAngleUp = beltAngleMiddle;
    beltAngleMiddle = beltAngleDown;
    beltAngleDown += beltAngle;
    beltRadiusUp = beltRadiusMiddle;
    beltRadiusMiddle = beltRadiusDown;
    beltRadiusDown = std::sin(beltAngleDown);
    beltZup = beltZmiddle;
    beltZmiddle = beltZdown;
    beltZdown = std::cos(beltAngleDown);
    bias += sectorAngleHalf;
  }
  return result;
}

}

#endif
