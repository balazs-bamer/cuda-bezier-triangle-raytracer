#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_MESHUTILS
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_MESHUTILS

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t

#include "stl_reader.h"
#include <map>
#include <list>
#include <array>
#include <deque>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
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
tContainer<std::array<uint32_t, 3u>> standardizeNormals(Mesh<tContainer> &aMesh);

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
void standardizeVertices(Mesh<tContainer> &aMesh) {        
  float smallestSide = std::numeric_limits<float>::max();
  for(auto &triangle : aMesh) {
    for(uint32_t i = 0u; i < 3u; ++i) {
      smallestSide = std::min(smallestSide, (triangle[i] - triangle[(i + 1u) % 3u]).norm());
    }
  }
  auto epsilon = smallestSide * cgStandardizeVerticesEpsilonFactor;       // Vertices closer to each other than this will be forced to one point.

  std::array<uint32_t, 3> maximums;
  using Multimap = std::multimap<float, std::pair<uint32_t, uint32_t>>;
  std::array<Multimap, 3u> projectedSorted;                               // Holds the projected indices sorted by a coordinate, separate for x, y and z.
  using Iterator = Multimap::const_iterator;
  std::array<std::deque<std::pair<Iterator, Iterator>>, 3u> intervals;    // Holds iterator pairs to the previous stuff in between the projections are closer to each other than epsilon.

  for(uint32_t dimension = 0u; dimension < 3u; ++dimension) {
    auto& currentProjectedSorted = projectedSorted[dimension];
    auto& currentIntervals = intervals[dimension];

    for(uint32_t indexFace = 0u; indexFace < aMesh.size(); ++indexFace) { // Project for each coordinate.
      auto &face = aMesh[indexFace];
      for(uint32_t indexVertex = 0u; indexVertex < 3u; ++indexVertex) {
        auto &vertex = face[indexVertex];
        currentProjectedSorted.emplace(std::make_pair(vertex[dimension], std::make_pair(indexFace, indexVertex)));
      }
    }
    bool was = false;
    float startValue;
    uint32_t max = 0u, counter;
    Iterator startIterator;
    for(Iterator i = currentProjectedSorted.cbegin(); i != currentProjectedSorted.cend(); ++i) {
      if(was) {
        auto now = i->first;
        if(now - startValue >= epsilon) {                                 // Determine intervals.
          currentIntervals.emplace_back(std::make_pair(startIterator, i));
          startValue = now;
          startIterator = i;
          max = std::max(max, counter);
          counter = 1u;
        }
        else {
          ++counter;
        }
      }
      else {
        startIterator = i;
        startValue = i->first;
        was = true;
        counter = 1u;
      }
    }
    Iterator last = currentProjectedSorted.cend();
    currentIntervals.emplace_back(std::make_pair(startIterator, last));
    max = std::max(max, counter);
    maximums[dimension] = max;
  }
  auto minIt = std::min_element(maximums.cbegin(), maximums.cend());  // Best coordinate is where the maximum interval length is minimal.
  auto whichDim = minIt - maximums.cbegin();
  auto& bestIntervals = intervals[whichDim];
  epsilon *= epsilon;                                                 // Avoid sqrt.

  for(auto& interval : bestIntervals) {                               // We need to check pairwise closeness only inside each interval.
    for(auto i = interval.first; i != interval.second; ++i) {
      auto &v1 = aMesh[i->second.first][i->second.second];
      for(auto j = interval.first; j != interval.second; ++j) {
        auto &v2 = aMesh[j->second.first][j->second.second];
        if((v1-v2).squaredNorm() < epsilon && (v1[0] < v2[0] || (v1[0] == v2[0] && v1[1] < v2[1]) || (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] < v2[2]))) {
          v1[0] = v2[0];                                              // Use lexicographic sorting to choose an unique one among the points too close to each other.
        }
        else { // nothing to do
        }
      }
    }
  }
}

struct PairHash {
  template<typename t1, typename t2>
  std::size_t operator()(std::pair<t1, t2> const &aPair) const {
    return std::hash<t1>{}(aPair.first) ^ (std::hash<t2>{}(aPair.second) << 1u);
  }
};

struct VertexHash {
  std::size_t operator()(Vertex const &aVertex) const {
    return std::hash<float>{}(aVertex[0]) ^ (std::hash<float>{}(aVertex[1]) << 1u) ^ (std::hash<float>{}(aVertex[2]) << 2u);
  }
};

Vertex normalize(Triangle &aFace, Vertex const &aDesiredVector) {
  auto normal = (aFace[1] - aFace[0]).cross(aFace[2] - aFace[0]).normalized();
  if(aDesiredVector.dot(normal) < 0.0f) {
    std::swap(aFace[0u], aFace[1u]);
    normal = -normal;
  }
  else { // Nothing to do
  }
  return normal;
}

template<template<typename> typename tContainer>                      // Should run after standardizeVertices.
tContainer<std::array<uint32_t, 3u>> standardizeNormals(Mesh<tContainer> &aMesh) {      // Makes all normalvectors (V[1]-V[0])x(V[2]-V[0]) point outwards.
  std::unordered_multimap<std::pair<uint32_t, uint32_t>, uint32_t, PairHash> edge2face; // Edge is identified by its two ordered vertex indices.
  std::unordered_map<Vertex, uint32_t, VertexHash> was;               // For presence testing, maps each vertex to its index in the deque below.
  std::deque<Vertex> vertices;                                        // Enables accessing each vertex by index. TODO consider if needed
  std::vector<std::array<uint32_t, 3u>> face2vertex;
  face2vertex.reserve(aMesh.size());
  float smallestX = std::numeric_limits<float>::max();
  uint32_t smallestIndex;
  for(uint32_t indexFace = 0u; indexFace < aMesh.size(); ++indexFace) {
    auto const &face = aMesh[indexFace];
    std::array<uint32_t, 3u> faceVertexIndices;                       // We collect here the vertex indices for the actual face.
    for(uint32_t indexInFace = 0u; indexInFace < 3u; ++indexInFace) {
      auto const vertex = face[indexInFace];
      uint32_t indexInVertices;
      auto found = was.find(vertex);
      if(found == was.end()) {
        indexInVertices = vertices.size();
        faceVertexIndices[indexInFace] = indexInVertices;
        was.emplace(std::make_pair(vertex, indexInVertices));
        vertices.push_back(vertex);                                   // Vertex indexing.
      }
      else {
        indexInVertices = found->second;
        faceVertexIndices[indexInFace] = indexInVertices;
      }
      if(vertex[0] < smallestX) {
        smallestX = vertex[0];
        smallestIndex = indexInVertices;
      }
      else { // Nothing to do
      }
    }
    face2vertex.push_back(faceVertexIndices);
    for(uint32_t i = 0u; i < 3u; ++i) {                               // Filling maps.
      auto low = faceVertexIndices[i];
      auto high = faceVertexIndices[(i + 1u) % 3u];
      if(low > high) {
        std::swap(low, high);
      }
      else { // Nothing to do
      }
std::cout << low << ' ' << high << ' ' << indexFace << std::endl;
      edge2face.emplace(std::make_pair(std::pair(low, high), indexFace));
    }
  }
  tContainer<std::array<uint32_t, 3u>> face2neighbour;                // Obtaining neighbours...
  std::unordered_set<uint32_t> facesAtSmallestX;                      // ...and the faces at the smallest x vertex
  for(uint32_t indexFace = 0u; indexFace < aMesh.size(); ++indexFace) {
    auto const &face = face2vertex[indexFace];
    std::array<uint32_t, 3u> neighbours;
    for(uint32_t indexInFace = 0u; indexInFace < 3u; ++indexInFace) {
      if(face[indexInFace] == smallestIndex) {
        facesAtSmallestX.insert(indexFace);
      }
      else { // Nothing to do
      }
      auto edge = std::make_pair(face[indexInFace], face[(indexInFace + 1u) % 3u]);
      if(edge.first > edge.second) {
        std::swap(edge.first, edge.second);
      }
      else { // Nothing to do
      }
      auto[begin, end] = edge2face.equal_range(edge);
      if(begin->second == indexFace) {                                // We need the other, neighbouring face.
        ++begin;
if(begin == end) throw 1;
      }
      else { // Nothing to do
      }
      neighbours[indexInFace] = begin->second;
    }
    face2neighbour.push_back(neighbours);
  }
  Vertex desiredVector;                                               // A known unit vector pointing outwards for a definite face.
  desiredVector << 1.0f, 0.0f, 0.0f;
  float maxAbsoluteDotProduct = std::numeric_limits<float>::max();
  uint32_t initialFaceIndex;                                          // First determine the initial face for which we want (f[1]-f[0])x(f[2]-f[0]) point outwards.
  for(auto const indexFace : facesAtSmallestX) {
    auto const &face = aMesh[indexFace];
    auto normal = (face[1] - face[0]).cross(face[2] - face[0]).normalized();
    auto absoluteDotProduct = std::abs(desiredVector.dot(normal));
    if(absoluteDotProduct > maxAbsoluteDotProduct) {
      maxAbsoluteDotProduct = absoluteDotProduct;
      initialFaceIndex = indexFace;
    }
    else { // Nothing to do
    }
  }
  struct KnownUnknownFacePair {
    uint32_t mKnownIndex;
    Vertex   mKnownNormal;
    uint32_t mUnknownIndex;
  };
  std::list<KnownUnknownFacePair> queue;                              // Faces with normals yet to be normalized.
  desiredVector = normalize(aMesh[initialFaceIndex], desiredVector);
  std::unordered_map<uint32_t, std::array<uint32_t, 3u>> remaining;
  for(uint32_t i = 0u; i < face2neighbour.size(); ++i) {
    remaining.emplace(std::make_pair(i, face2neighbour[i]));
  }
  for(auto const indexFace : remaining[initialFaceIndex]) {
    queue.emplace_back(KnownUnknownFacePair{initialFaceIndex, desiredVector, indexFace});
  }
  remaining.erase(initialFaceIndex);
  while(!queue.empty()) {
    auto actualPair = queue.front();
    queue.pop_front();
    desiredVector = normalize(aMesh[actualPair.mUnknownIndex], actualPair.mKnownNormal);
    for(auto const indexFace : remaining[actualPair.mUnknownIndex]) {
      if(remaining.find(indexFace) != remaining.end()) {
        queue.emplace_back(KnownUnknownFacePair{actualPair.mUnknownIndex, desiredVector, indexFace});
      }
      else { // Nothing to do
      }
    }
    remaining.erase(actualPair.mKnownIndex);
  }
  return face2neighbour;
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
