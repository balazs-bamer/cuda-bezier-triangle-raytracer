#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH

#include "mesh.h"
#include "bezierTriangle.h"
#include <functional>
#include <numeric>

template<typename tReal>
class BezierMesh final {
private:
  static constexpr tReal csBezierHeightPerPerimeterLimit = 0.1f;      // TODO consider

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
  Mesh<tReal> splitThickBezierTriangles() const;

/*Mesh<tReal> stuff;
Mesh<tReal> const& getStuff() const { return stuff; }*/

private:
  void setMissingFields(Mesh<tReal> const &aMesh, MissingFieldsMethod aMissingFieldsMethod);
  void append2split(Mesh<tReal> &aResult, Triangle const &aOriginalTriangle, uint8_t const aSplit) const;
  void append3split(Mesh<tReal> &aResult, Triangle const &aOriginalTriangle, uint8_t const aSplit) const;
  void append4split(Mesh<tReal> &aResult, Triangle const &aOriginalTriangle, uint8_t const aSplit) const;
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

/*#include<fstream>
#include<iostream>
#include<iomanip>*/

template<typename tReal>
BezierMesh<tReal>::BezierMesh(Mesh<tReal> const &aMesh) : mOriginalNeighbours(aMesh.getFace2neighbours()) {
  mMesh.clear();
  mMesh.reserve(aMesh.size() * 3u);
  auto const &mesh                  = aMesh.getMesh();
  auto const &vertex2averageNormals = aMesh.getVertex2averageNormals();
//std::map<Plane, uint32_t> neighPlanes;
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
/*std::cout << std::setw(3) << indexFace << std::setw(3) << indexVertex << " =" << std::setw(3) << mMesh.size() <<
std::setprecision(4) << std::setw(9) << (::abs(originalCommonVertex0(0)) < 0.001f ? 0.0f : originalCommonVertex0(0)) <<
std::setprecision(4) << std::setw(9) << (::abs(originalCommonVertex0(1)) < 0.001f ? 0.0f : originalCommonVertex0(1)) <<
std::setprecision(4) << std::setw(9) << (::abs(originalCommonVertex0(2)) < 0.001f ? 0.0f : originalCommonVertex0(2)) << "  -" <<
std::setprecision(4) << std::setw(9) << (::abs(originalCommonVertex1(0)) < 0.001f ? 0.0f : originalCommonVertex1(0)) <<
std::setprecision(4) << std::setw(9) << (::abs(originalCommonVertex1(1)) < 0.001f ? 0.0f : originalCommonVertex1(1)) <<
std::setprecision(4) << std::setw(9) << (::abs(originalCommonVertex1(2)) < 0.001f ? 0.0f : originalCommonVertex1(2)) << "  |" <<
std::setprecision(4) << std::setw(9) << (::abs(planeBetweenOriginalNeighbours.mNormal(0)) < 0.001f ? 0.0f : planeBetweenOriginalNeighbours.mNormal(0)) <<
std::setprecision(4) << std::setw(9) << (::abs(planeBetweenOriginalNeighbours.mNormal(1)) < 0.001f ? 0.0f : planeBetweenOriginalNeighbours.mNormal(1)) <<
std::setprecision(4) << std::setw(9) << (::abs(planeBetweenOriginalNeighbours.mNormal(2)) < 0.001f ? 0.0f : planeBetweenOriginalNeighbours.mNormal(2)) << " ," <<
std::setprecision(4) << std::setw(9) << (::abs(planeBetweenOriginalNeighbours.mConstant) < 0.001f ? 0.0f : planeBetweenOriginalNeighbours.mConstant) << "  |" <<
std::setprecision(4) << std::setw(9) << (::abs(averageNormal0(0)) < 0.001f ? 0.0f : averageNormal0(0)) <<
std::setprecision(4) << std::setw(9) << (::abs(averageNormal0(1)) < 0.001f ? 0.0f : averageNormal0(1)) <<
std::setprecision(4) << std::setw(9) << (::abs(averageNormal0(2)) < 0.001f ? 0.0f : averageNormal0(2)) << "  -" <<
std::setprecision(4) << std::setw(9) << (::abs(averageNormal1(0)) < 0.001f ? 0.0f : averageNormal1(0)) <<
std::setprecision(4) << std::setw(9) << (::abs(averageNormal1(1)) < 0.001f ? 0.0f : averageNormal1(1)) <<
std::setprecision(4) << std::setw(9) << (::abs(averageNormal1(2)) < 0.001f ? 0.0f : averageNormal1(2)) << '\n';
Plane v = planeBetweenOriginalNeighbours;
v.mNormal(0) = ::round(v.mNormal(0) * 128.0f) / 128.0f;
v.mNormal(1) = ::round(v.mNormal(1) * 128.0f) / 128.0f;
v.mNormal(2) = ::round(v.mNormal(2) * 128.0f) / 128.0f;
v.mConstant = ::round(v.mConstant * 128.0f) / 128.0f;
Plane vm;
vm.mNormal = -v.mNormal;
vm.mConstant = -v.mConstant;
auto found = neighPlanes.find(v);
auto foundm = neighPlanes.find(vm);
if(found != neighPlanes.end()) {
  ++found->second;
}
else if(foundm != neighPlanes.end()) {
  ++foundm->second;
}
else {
  neighPlanes.insert(std::make_pair(v, 1u));
}
Vertex v0 = originalCommonVertex0;
Vertex v1 = originalCommonVertex0 + averageNormal0;
Vertex v2 = planeBetweenOriginalNeighbours * v1;
Vertex v3 = originalCommonVertex1;
if((v2-v1).norm()<0.01f) {
  v2 = v1 + (v1-v0).cross(v3-v0)/(v3-v0).norm()/10.0f;
}
stuff.push_back({v0, v1, v2});
stuff.push_back({v0, v1, v3});
stuff.push_back({v0, v2, v3});
stuff.push_back({v1, v2, v3});*/
      mMesh.emplace_back(BezierTriangle(originalCommonVertex0, originalCommonVertex1, originalCentroid,
                                               averageNormal0, averageNormal1, planeBetweenOriginalNeighbours,
                                               newNeighbourIndices));
    }
  }
/*for(auto const &[key, value] : neighPlanes) {
  std::cout << value << ' ';
}
std::cout << '\n';*/
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
    Vertex bezierAboveOriginalCentroid = mMesh[indexOriginal * 3u].interpolateAboveOriginalCentroid();
    Triangle original { mMesh[indexOriginal * 3u].getControlPoint(0u),
                        mMesh[indexOriginal * 3u + 1u].getControlPoint(0u),
                        mMesh[indexOriginal * 3u + 2u].getControlPoint(0u) };
    if((original / bezierAboveOriginalCentroid / getPerimeter(original)) > csBezierHeightPerPerimeterLimit) {
      splitSides[indexOriginal] = csSplitAll;
      auto const &neigh = mOriginalNeighbours[indexOriginal];
      for(uint32_t side = 0u; side < 3u; ++side) {
        splitSides[neigh.mFellowTriangles[side]] |= csSplitMask[neigh.mFellowCommonSideStarts[side]];
      }
    }
    else { // Nothing to do
    }
  }

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
      append2split(result, original, split);
    }
    else if(newCount == 3u) {
      append3split(result, original, split);
    }
    else { // newCount == 4u
      append4split(result, original, split);
    }
  }
  return result;
}

template<typename tReal>
void BezierMesh<tReal>::append2split(Mesh<tReal> &aResult, Triangle const &aOriginalTriangle, uint8_t const aSplit) const {
  static constexpr uint8_t csIndexFor2onSide[8u] = { 3u, 0u, 1u, 3u, 2u, 3u, 3u, 3u };
  uint32_t index2 = csIndexFor2onSide[aSplit];

}

template<typename tReal>
void BezierMesh<tReal>::append3split(Mesh<tReal> &aResult, Triangle const &aOriginalTriangle, uint8_t const aSplit) const {

}

template<typename tReal>
void BezierMesh<tReal>::append4split(Mesh<tReal> &aResult, Triangle const &aOriginalTriangle, uint8_t const aSplit) const {

}

#endif
