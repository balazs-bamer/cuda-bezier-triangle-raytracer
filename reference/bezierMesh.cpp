#include "bezierMesh.h"


BezierMesh::BezierMesh(Mesh const &aMesh) : mOriginalNeighbours(aMesh.getFace2neighbours()) {
  mMesh.clear();
  mMesh.reserve(aMesh.size() * 3u);
  auto const &mesh                  = aMesh.getMesh();
  auto const &vertex2averageNormals = aMesh.getVertex2averageNormals();
  for(uint32_t indexFace = 0u; indexFace < aMesh.size(); ++indexFace) {
    auto const &neigh = mOriginalNeighbours[indexFace];
    auto const &originalTriangle = aMesh[indexFace];
    auto const originalCentroid = (originalTriangle[0u] + originalTriangle[1u] + originalTriangle[2u]) / 3.0f;
    auto normal = util::getNormal(originalTriangle).normalized();
    for(uint32_t indexVertex = 0u; indexVertex < 3u; ++indexVertex) { // Appending a regular triangle means doing Clough-Tocher split on it and appending each small one.
      auto const &originalCommonVertex0 = originalTriangle[indexVertex];
      auto const &originalCommonVertex1 = originalTriangle[(indexVertex + 1u) % 3u];
      auto const &averageNormal0 = vertex2averageNormals.at(originalCommonVertex0);
      auto const &averageNormal1 = vertex2averageNormals.at(originalCommonVertex1);
                                                                                 // Neighbour of edge (index, index + 1)
      Plane planeBetweenOriginalNeighbours = Plane::createFrom1vector2points(normal + util::getNormal(aMesh[neigh.mFellowTriangles[indexVertex]]).normalized(),
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

void BezierMesh::setMissingFields(Mesh const &aMesh, MissingFieldsMethod aMissingFieldsMethod) {
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

Mesh BezierMesh::interpolate(int32_t const aDivisor) const {
  Mesh result;
  Triangle barycentric{Vertex{{1.0f, 0.0f, 0.0f}}, Vertex{{0.0f, 1.0f, 0.0f}}, Vertex{{0.0f, 0.0f, 1.0f}}};
  util::divide(barycentric, aDivisor, [&result, this](Triangle && aNewBary){
    for(auto const &bezier : mMesh) {
      result.push_back({ bezier.interpolate(aNewBary[0u]),
                         bezier.interpolate(aNewBary[1u]),
                         bezier.interpolate(aNewBary[2u]) });
    }
  } );
  return result;
}

std::vector<Vertex> BezierMesh::dumpControlPoints() const {
  std::vector<Vertex> result;
  result.reserve(mMesh.size() * BezierTriangle::csControlPointsSize);
  for(auto const &bezier : mMesh) {
    for(uint32_t i = 0u; i < BezierTriangle::csControlPointsSize; ++i) {
      result.push_back(bezier.getControlPoint(i));
    }
  }
  return result;
}

Mesh BezierMesh::splitThickBezierTriangles() const {
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
    float max = ::abs(plane.distance(mMesh[indexOriginal * 3u].interpolateAboveOriginalCentroid()));
    for(uint32_t i = 0; i < 3u; ++i) {
      for(auto const ratio : csSampleRatiosOriginalSide) {
        max = std::max(max, ::abs(plane.distance(mMesh[indexSplit + i].interpolate(ratio, 1.0f - ratio, 0.0f))));
      }
    }
    if(max / util::getPerimeter(original) > csBezierHeightPerPerimeterLimit) {
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
  Mesh result;
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
void BezierMesh::append2split(Mesh &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle, uint8_t const aSplit) const {
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
void BezierMesh::append3split(Mesh &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle, uint8_t const aSplit) const {
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
void BezierMesh::append4split(Mesh &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle) const {
  Triangle middle;
  for(uint32_t i = 0u; i < 3u; ++i) {
    middle[i] = interpolate(aIndexBase + i, 0.5f, 0.5f, 0.0f);
  }
  aResult.push_back(middle);
  for(uint32_t i = 0u; i < 3u; ++i) {
    aResult.push_back(Triangle{ aOriginalTriangle[i], middle[i], middle[(i + 2u) % 3u] });
  }
}

Vertex BezierMesh::interpolate(uint32_t const aIndex, float const aBary0, float const aBary1, float const aBary2) const {
  auto const &triangle = mMesh[aIndex];
  return triangle.interpolate(aBary0, aBary1, aBary2) * csSplitBezierInterpolateFactor +
         triangle.interpolateLinear(aBary0, aBary1, aBary2) * (1.0f - csSplitBezierInterpolateFactor);
}

BezierIntersection BezierMesh::intersect(Ray const &aRay) const {
  // Brute-force for now
  BezierIntersection result;
  result.mIntersection.mDistance = std::numeric_limits<float>::max();
  result.mWhat = BezierIntersection::What::cNone;
  for(auto const &bezier : mMesh) {
    auto candidate = bezier.intersect(aRay, BezierTriangle::LimitPlaneIntersection::cThis);
    if(candidate.mWhat == BezierIntersection::What::cFollowSide0 ||
       candidate.mWhat == BezierIntersection::What::cFollowSide1 ||
       candidate.mWhat == BezierIntersection::What::cFollowSide2) {
      candidate = mMesh[bezier.getNeighbours()[static_cast<uint32_t>(candidate.mWhat)]].intersect(aRay, BezierTriangle::LimitPlaneIntersection::cNone);
    }
    else { // Nothing to do
    }
    if(candidate.mWhat == BezierIntersection::What::cIntersect && candidate.mIntersection.mDistance < result.mIntersection.mDistance) {
      result = candidate;
    }
    else { // Nothing to do
    }
  }
  return result;
}
