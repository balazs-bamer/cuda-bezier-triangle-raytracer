#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH

#include "mesh.h"
#include "bezierTriangle.h"
#include <functional>
#include <numeric>

template<typename tReal>
class BezierMesh final {
private:
  static constexpr std::array<tReal, 3u> csSampleRatiosOriginalSide = { 0.25f, 0.5f, 0.75f };
  static constexpr tReal csBezierHeightPerPerimeterLimit = 0.03f;
  static constexpr tReal csSplitBezierInterpolateFactor  = 0.7f;   // For triangle splitting, new vertex is computed as
                                                                   // barycentricSplit * (1.0 - csSBIF) + interpolate(barycentricSplit) * csSBIF

  std::vector<BezierTriangle<tReal>>    mMesh;
  typename Mesh<tReal>::Face2neighbours mOriginalNeighbours;
  
  using Vector              = ::Vector<tReal>;
  using Vertex              = ::Vertex<tReal>;
  using Matrix              = ::Matrix<tReal>;
  using Triangle            = ::Triangle<tReal>;
  using Plane               = ::Plane<tReal>;
  using BezierTriangle      = ::BezierTriangle<tReal>;
  using Neighbours          = typename Mesh<tReal>::Neighbours;
  using MissingFieldsMethod = std::function<void(BezierTriangle &, Vertex const &aOriginalCentroid, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious)>;

public:
  BezierMesh(Mesh<tReal> const &aMesh);

  auto size() const { return mMesh.size(); }
  auto cbegin() const { return mMesh.cbegin(); }
  auto cend() const { return mMesh.cend(); }
  auto const &operator[](uint32_t const aI) const { return mMesh[aI]; }

  Mesh<tReal> interpolate(int32_t const aDivisor) const;
  std::vector<Vertex> dumpControlPoints() const;
  Mesh<tReal> splitThickBezierTriangles() const;        // TODO consider automatic mechanism to iterate until either
                                                        // - no more split or under a percentage
                                                        // - absolute height left is below a value
                                                        // - minimal triangle perimeter is below a value or percentage of the original

/*Mesh<tReal> stuff;
Mesh<tReal> const& getStuff() const { return stuff; }*/

private:
  void setMissingFields(Mesh<tReal> const &aMesh, MissingFieldsMethod aMissingFieldsMethod);
  void append2split(Mesh<tReal> &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle, uint8_t const aSplit) const;
  void append3split(Mesh<tReal> &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle, uint8_t const aSplit) const;
  void append4split(Mesh<tReal> &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle) const;
  Vertex interpolate(uint32_t const aIndex, tReal const aBary0, tReal const aBary1, tReal const aBary2) const;
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tReal>
BezierMesh<tReal>::BezierMesh(Mesh<tReal> const &aMesh) : mOriginalNeighbours(aMesh.getFace2neighbours()) {
  mMesh.clear();
  mMesh.reserve(aMesh.size() * 3u);
  auto const &mesh                  = aMesh.getMesh();
  auto const &vertex2averageNormals = aMesh.getVertex2averageNormals();
  for(uint32_t indexFace = 0u; indexFace < aMesh.size(); ++indexFace) {
    auto const &neigh = mOriginalNeighbours[indexFace];
    auto const &originalTriangle = aMesh[indexFace];
    auto const originalCentroid = (originalTriangle[0u] + originalTriangle[1u] + originalTriangle[2u]) / 3.0f;
    auto normal = getNormal(originalTriangle).normalized();
    for(uint32_t indexVertex = 0u; indexVertex < 3u; ++indexVertex) { // Appending a regular triangle means doing Clough-Tocher split on it and appending each small one.
      auto const &originalCommonVertex0 = originalTriangle[indexVertex];
      auto const &originalCommonVertex1 = originalTriangle[(indexVertex + 1u) % 3u];
      auto const &averageNormal0 = vertex2averageNormals.at(originalCommonVertex0);
      auto const &averageNormal1 = vertex2averageNormals.at(originalCommonVertex1);
                                                                                 // Neighbour of edge (index, index + 1)
      Plane planeBetweenOriginalNeighbours = Plane::createFrom1vector2points(normal + getNormal(aMesh[neigh.mFellowTriangles[indexVertex]]).normalized(),
                                                                             originalCommonVertex0, originalCommonVertex1);
      auto currentBase = indexFace * 3u;
      std::array<uint32_t, 3u> newNeighbourIndices({ 3u * neigh.mFellowTriangles[indexVertex] + neigh.mFellowCommonSideStarts[indexVertex],
                                                     currentBase + (indexVertex + 1u) % 3u,
                                                     currentBase + (indexVertex + 2u) % 3u });
      mMesh.emplace_back(BezierTriangle(originalCommonVertex0, originalCommonVertex1, originalCentroid,
                                               averageNormal0, averageNormal1, planeBetweenOriginalNeighbours,
                                               newNeighbourIndices));
    }
  }
  setMissingFields(aMesh, &BezierTriangle::setMissingFields1);
  setMissingFields(aMesh, &BezierTriangle::setMissingFields2);
  setMissingFields(aMesh, &BezierTriangle::setMissingFields3);
}

template<typename tReal>
void BezierMesh<tReal>::setMissingFields(Mesh<tReal> const &aMesh, MissingFieldsMethod aMissingFieldsMethod) {
  Vertex originalCentroid;
  for(uint32_t i = 0u; i < mMesh.size(); ++i) {
    auto subTriangleIndex = i % 3u;
    auto indexBase = i - subTriangleIndex;
    auto const &triangleNext = mMesh[indexBase + (subTriangleIndex + 1u) % 3u];
    auto const &trianglePrevious = mMesh[indexBase + (subTriangleIndex + 2u) % 3u];
    if(subTriangleIndex % 3u == 0u) {
      auto const &originalTriangle = aMesh[indexBase / 3u];
      originalCentroid = (originalTriangle[0u] + originalTriangle[1u] + originalTriangle[2u]) / 3.0f;
    }
    else { // Nothing to do
    }
    aMissingFieldsMethod(mMesh[i], originalCentroid, triangleNext, trianglePrevious);
  }
}

#include<iostream> // TODO remove

template<typename tReal>
Mesh<tReal> BezierMesh<tReal>::interpolate(int32_t const aDivisor) const {
  Mesh<tReal> result;
  Triangle barycentric{Vertex{{1.0f, 0.0f, 0.0f}}, Vertex{{0.0f, 1.0f, 0.0f}}, Vertex{{0.0f, 0.0f, 1.0f}}};
  divide(barycentric, aDivisor, [&result, this](Triangle && aNewBary){
    for(auto const &bezier : mMesh) {
      result.push_back({ bezier.interpolate(aNewBary[0u]),
                         bezier.interpolate(aNewBary[1u]),
                         bezier.interpolate(aNewBary[2u]) });
    }
  } );
  return result;
}

template<typename tReal>
std::vector<Vertex<tReal>> BezierMesh<tReal>::dumpControlPoints() const {
  std::vector<Vertex> result;
  result.reserve(mMesh.size() * BezierTriangle::csControlPointsSize);
  for(auto const &bezier : mMesh) {
    for(uint32_t i = 0u; i < BezierTriangle::csControlPointsSize; ++i) {
      result.push_back(bezier.getControlPoint(i));
    }
  }
  return result;
}

template<typename tReal>
Mesh<tReal> BezierMesh<tReal>::splitThickBezierTriangles() const {
  static constexpr uint8_t csSplitMask[3u]  = { 1u, 2u, 4u };
  static constexpr uint8_t csSplitAll       = 7u;
  static constexpr uint8_t csSplitCount[8u] = { 1u, 2u, 2u, 3u, 2u, 3u, 3u, 4u };
  std::vector<uint8_t> splitSides(mOriginalNeighbours.size(), 0u);

  for(uint32_t indexOriginal = 0u; indexOriginal < mOriginalNeighbours.size(); ++indexOriginal) {
    auto indexSplit = indexOriginal * 3u;
    Triangle original { mMesh[indexSplit].getControlPoint(0u),
                        mMesh[indexSplit + 1u].getControlPoint(0u),
                        mMesh[indexSplit + 2u].getControlPoint(0u) };
    auto plane = Plane::createFromTriangle(original);
    tReal max = ::abs(plane.distance(mMesh[indexOriginal * 3u].interpolateAboveOriginalCentroid()));
    for(uint32_t i = 0; i < 3u; ++i) {
      for(auto const ratio : csSampleRatiosOriginalSide) {
        max = std::max(max, ::abs(plane.distance(mMesh[indexSplit + i].interpolate(ratio, 1.0f - ratio, 0.0f))));
      }
    }
    if(max / getPerimeter(original) > csBezierHeightPerPerimeterLimit) {
      splitSides[indexOriginal] = csSplitAll;
      auto const &neigh = mOriginalNeighbours[indexOriginal];
      for(uint32_t side = 0u; side < 3u; ++side) {
        splitSides[neigh.mFellowTriangles[side]] |= csSplitMask[neigh.mFellowCommonSideStarts[side]];
      }
    }
    else { // Nothing to do
    }
  }
  // TODO go over original edges to see if that edge is not split but other edges of the 2 containing
  // original triangles are. In this case this edge should be split as well to prevent it become a depressed edge.
  // Or perhaps inspect each remaining original edge if the convexity changes there.
  uint32_t newCount = 0u;
  newCount = std::accumulate(splitSides.cbegin(), splitSides.cend(), newCount, [](uint32_t aSofar, uint8_t aNew){ return aSofar + csSplitCount[aNew]; });
  Mesh<tReal> result;
  result.reserve(newCount);
  for(uint32_t indexOriginal = 0u; indexOriginal < mOriginalNeighbours.size(); ++indexOriginal) {
    Triangle original { mMesh[indexOriginal * 3u].getControlPoint(0u),
                        mMesh[indexOriginal * 3u + 1u].getControlPoint(0u),
                        mMesh[indexOriginal * 3u + 2u].getControlPoint(0u) };
    uint8_t split = splitSides[indexOriginal];
    newCount = csSplitCount[split];
    if(newCount == 1u) {
      result.push_back(original);
    }
    else if(newCount == 2u) {
      append2split(result, indexOriginal * 3u, original, split);
    }
    else if(newCount == 3u) {
      append3split(result, indexOriginal * 3u, original, split);
    }
    else { // newCount == 4u
      append4split(result, indexOriginal * 3u, original);
    }
  }
  return result;
}

//       1
//      /|\
//     / | \
//    /  |  \
//   /   |   \
//  /____|____\
// 0           2   index2 == 2
//
template<typename tReal>
void BezierMesh<tReal>::append2split(Mesh<tReal> &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle, uint8_t const aSplit) const {
  static constexpr uint8_t csIndexFor2onSide[8u] = { 3u, 0u, 1u, 3u, 2u, 3u, 3u, 3u };
  uint32_t index2 = csIndexFor2onSide[aSplit];
  auto splitVertex = interpolate(aIndexBase + index2, 0.5f, 0.5f, 0.0f);
  uint32_t indexAfter2 = (index2 + 1u) % 3u;
  uint32_t indexBefore2 = (index2 + 2u) % 3u;
  aResult.push_back(Triangle{ aOriginalTriangle[indexAfter2], aOriginalTriangle[indexBefore2], splitVertex });
  aResult.push_back(Triangle{ aOriginalTriangle[indexBefore2], aOriginalTriangle[index2], splitVertex });
}

//       1
//      /|\
//     / | \SA1
//    /  | /\
//   /   |/  \     or not 1-SB1 but 0-SA1, whichever is shorter 
//  /____/____\
// 0    SB1    2   index1 == 0
//
template<typename tReal>
void BezierMesh<tReal>::append3split(Mesh<tReal> &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle, uint8_t const aSplit) const {
  static constexpr uint8_t csIndexFor1onSide[8u] = { 3u, 3u, 3u, 2u, 3u, 1u, 0u, 3u };
  uint32_t index1 = csIndexFor1onSide[aSplit];
  uint32_t indexAfter1 = (index1 + 1u) % 3u;
  uint32_t indexBefore1 = (index1 + 2u) % 3u;
  auto splitVertexBefore1 = interpolate(aIndexBase + indexBefore1, 0.5f, 0.5f, 0.0f);
  auto splitVertexAfter1  = interpolate(aIndexBase + indexAfter1, 0.5f, 0.5f, 0.0f);
  aResult.push_back(Triangle{ aOriginalTriangle[indexBefore1], splitVertexBefore1, splitVertexAfter1 });
  if((aOriginalTriangle[indexAfter1] - splitVertexBefore1).norm() < (aOriginalTriangle[index1] - splitVertexAfter1).norm()) {
    aResult.push_back(Triangle{ aOriginalTriangle[indexAfter1], splitVertexAfter1, splitVertexBefore1 });
    aResult.push_back(Triangle{ aOriginalTriangle[index1], aOriginalTriangle[indexAfter1], splitVertexBefore1 });
  }
  else {
    aResult.push_back(Triangle{ aOriginalTriangle[indexAfter1], splitVertexAfter1, aOriginalTriangle[index1] });
    aResult.push_back(Triangle{ aOriginalTriangle[index1], splitVertexAfter1, splitVertexBefore1 });
  }
}

//        1
//       /\
//      /  \
//   S0/____\S1
//    /\    /\
//   /  \  /  \
//  /____\/____\
// 0     S2     2
//
template<typename tReal>
void BezierMesh<tReal>::append4split(Mesh<tReal> &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle) const {
  Triangle middle;
  for(uint32_t i = 0u; i < 3u; ++i) {
    middle[i] = interpolate(aIndexBase + i, 0.5f, 0.5f, 0.0f);
  }
  aResult.push_back(middle);
  for(uint32_t i = 0u; i < 3u; ++i) {
    aResult.push_back(Triangle{ aOriginalTriangle[i], middle[i], middle[(i + 2u) % 3u] });
  }
}

template<typename tReal>
Vertex<tReal> BezierMesh<tReal>::interpolate(uint32_t const aIndex, tReal const aBary0, tReal const aBary1, tReal const aBary2) const {
  auto const &triangle = mMesh[aIndex];
  return triangle.interpolate(aBary0, aBary1, aBary2) * csSplitBezierInterpolateFactor +
         triangle.interpolateLinear(aBary0, aBary1, aBary2) * (1.0f - csSplitBezierInterpolateFactor);
}

#endif
