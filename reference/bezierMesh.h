#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH

#include "mesh.h"
#include "bezierTriangle.h"
#include <functional>
#include <numeric>


class BezierMesh final {
private:
  static constexpr std::array<float, 3u> csSampleRatiosOriginalSide = { 0.25f, 0.5f, 0.75f };
  static constexpr float                 csBezierHeightPerPerimeterLimit = 0.03f;
  static constexpr float                 csSplitBezierInterpolateFactor  = 0.7f;   // For triangle splitting, new vertex is computed as
                                                                                   // barycentricSplit * (1.0 - csSBIF) + interpolate(barycentricSplit) * csSBIF

  std::vector<BezierTriangle>    mMesh;
  typename Mesh::Face2neighbours mOriginalNeighbours;
  
  using Neighbours          = typename Mesh::Neighbours;
  using MissingFieldsMethod = std::function<void(BezierTriangle &, Vertex const &aOriginalCentroid, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious)>;

public:
  BezierMesh(Mesh const &aMesh);

  auto size() const { return mMesh.size(); }
  auto cbegin() const { return mMesh.cbegin(); }
  auto cend() const { return mMesh.cend(); }
  auto const &operator[](uint32_t const aI) const { return mMesh[aI]; }

  Mesh interpolate(int32_t const aDivisor) const;
  std::vector<Vertex> dumpControlPoints() const;
  Mesh splitThickBezierTriangles() const;        // TODO consider automatic mechanism to iterate until either
                                                        // - no more split or under a percentage
                                                        // - absolute height left is below a value
                                                        // - minimal triangle perimeter is below a value or percentage of the original
  BezierIntersection intersect(Ray const &aRay) const;

private:
  void setMissingFields(Mesh const &aMesh, MissingFieldsMethod aMissingFieldsMethod);
  void append2split(Mesh &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle, uint8_t const aSplit) const;
  void append3split(Mesh &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle, uint8_t const aSplit) const;
  void append4split(Mesh &aResult, uint32_t const aIndexBase, Triangle const &aOriginalTriangle) const;
  Vertex interpolate(uint32_t const aIndex, float const aBary0, float const aBary1, float const aBary2) const;
};

#endif
