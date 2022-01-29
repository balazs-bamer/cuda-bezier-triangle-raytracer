#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER

#include "3dGeomUtil.h"


struct BezierIntersection final {
  enum class What : uint32_t {
    cFollowSide0 = 0u,
    cFollowSide1 = 1u,
    cFollowSide2 = 2u,
    cNone        = 3u,
    cIntersect   = 4u
  };

  Intersection mIntersection;
  Vertex       mBarycentric;
  Vector       mNormal;
  What         mWhat;
};

class BezierTriangle final {    // Cubic Bezier triangle
public:
  enum class LimitPlaneIntersection : uint32_t {
    cThis = 0u,
    cNone = 1u
  };

  static constexpr uint32_t csControlPointsSize                       = 10u;
private:
  static constexpr uint32_t csControlIndexOriginalVertex0             = 0u;
  static constexpr uint32_t csControlIndexOriginalVertex1             = 1u;
  static constexpr uint32_t csControlIndexAboveOriginalCentroid       = 2u;
  static constexpr uint32_t csControlIndexOnOriginalSide0             = 3u;
  static constexpr uint32_t csControlIndexOnOriginalSide1             = 4u;
  static constexpr uint32_t csControlIndexOnSideToOriginalCentroid0   = 5u;
  static constexpr uint32_t csControlIndexOnSideToOriginalCentroid1   = 6u;
  static constexpr uint32_t csControlIndexOnSideFromOriginalCentroid0 = 7u;
  static constexpr uint32_t csControlIndexOnSideFromOriginalCentroid1 = 8u;
  static constexpr uint32_t csControlIndexMiddle                      = 9u;

  static constexpr uint32_t csControlIndex300 = csControlIndexOriginalVertex0;
  static constexpr uint32_t csControlIndex030 = csControlIndexOriginalVertex1;
  static constexpr uint32_t csControlIndex003 = csControlIndexAboveOriginalCentroid;
  static constexpr uint32_t csControlIndex210 = csControlIndexOnOriginalSide0;
  static constexpr uint32_t csControlIndex120 = csControlIndexOnOriginalSide1;
  static constexpr uint32_t csControlIndex021 = csControlIndexOnSideToOriginalCentroid0;
  static constexpr uint32_t csControlIndex012 = csControlIndexOnSideToOriginalCentroid1;
  static constexpr uint32_t csControlIndex102 = csControlIndexOnSideFromOriginalCentroid0;
  static constexpr uint32_t csControlIndex201 = csControlIndexOnSideFromOriginalCentroid1;
  static constexpr uint32_t csControlIndex111 = csControlIndexMiddle;

  static constexpr float    csProportionControlOnOriginalSide           = 0.291f;
  static constexpr float    csProportionControlOnOriginalVertexCentroid = 0.304f;
  static constexpr float    csProportionControlOnOriginalMedian         = 0.2f;
  static constexpr float    csHeightSafetyFactor                        = 1.33333333f;
  static constexpr float    csOneThird                                  = 1.0 / 3.0;
  static constexpr uint32_t csRootSearchIterations                      = 4u;
  static constexpr int32_t  csHeightSampleDivisor                       = 5;
  static constexpr float    csMaxIntersectionDistanceFromRay            = 0.01f;
  static constexpr float    csMinimalRayDistance                        = 1.0f;
  static constexpr float    csIntersectionEstimationEpsilon             = 0.000001f;

  Plane                    mUnderlyingPlane;         // Plane span by mControlPoints[0-2], normal is normalized and points outwards.
  std::array<Plane, 3u>    mNeighbourDividerPlanes;  // For indices 1 and 2, these are about the plane going through each edge in direction the average of the adjacent side normals.
                                                     // Index 0 will be aPlaneBetweenOriginalNeighbours
                                                     // Their normal is such that for all planes any point P belonging to this Bezier triangle has plane.distance(P) >= 0
  std::array<uint32_t, 3u> mNeighbours;              // With new indices after Clough-Tocher split.
  // 0-1 identical to original triangle vertices, named aOriginalCommonVertex0, aOriginalCommonVertex1 in constructor
  // 2 somewhere above the middle of the original triangle center, to be calculated
  // 3-8 somewhere above sides of triangle (0-2): 3i+3, 3i+4 on the same side as (i, (i+1)%3), to be calculated
  // 9 in the middle, to be calculated
  std::array<Vertex, csControlPointsSize>  mControlPoints;
  // Perhaps not needed, because the iterative method to find ray and Bezier surface intersection takes long. std::array<Vertex, 12u> mDerivativeControlPoints; // T most probably more than 2*6 for each partial derivative.
  Matrix                   mBarycentricInverse;      // T = (v1 v2 v3), b = barycentric coefficient column, v = point on plane, v=Tb, b = mBI * v
                                                     // multiplication by the inverse proven to be much faster than solving the linear equation under Eigen.
  float                    mHeightInside;            // Sampled biggest distance of the Bezier triangle measured from the underlying triangle, < 0
  float                    mHeightOutside;           // this is > 0
  Vector                   mBezierDerivativeDirectionVectorA; // We calculate directional derivatives in two perpendicular directions
  Vector                   mBezierDerivativeDirectionVectorB; // and their cross product will yield the Bezier surface normal

public:
  // Vertices in argument are in order such that the normal points to the desired direction.
  // Neighbour i is neighbour of edge (i, i + 1)
  BezierTriangle(Vertex const &aOriginalCommonVertex0, Vertex const &aOriginalCommonVertex1, Vertex const &aOriginalCentroid,
                 Vector const &aAverageNormal0, Vector const &aAverageNormal1, Plane const &aPlaneBetweenOriginalNeighbours,
                 std::array<uint32_t, 3u> const &aNeighbourIndices);

  void setMissingFields1(Vertex const &aOriginalCentroid, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious);
  void setMissingFields2(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &);
  void setMissingFields3(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious);

  Vertex getControlPoint(uint32_t const aI) const { return mControlPoints[aI]; }
  std::array<uint32_t, 3u> getNeighbours() const { return mNeighbours; }           // TODO remove

  Vertex interpolateLinear(float const aBary0, float const aBary1, float const aBary2) const;
  Vertex interpolateLinear(float const aBary0, float const aBary1) const { return interpolateLinear(aBary0, aBary1, 1.0f - aBary0 - aBary1); }
  Vertex interpolateLinear(Vertex const &aBary) const { return interpolateLinear(aBary(0u), aBary(1u), aBary(2u)); }

  Vertex interpolate(float const aBary0, float const aBary1, float const aBary2) const;
  Vertex interpolate(float const aBary0, float const aBary1) const { return interpolate(aBary0, aBary1, 1.0f - aBary0 - aBary1); }
  Vertex interpolate(Vertex const &aBary) const { return interpolate(aBary(0u), aBary(1u), aBary(2u)); }
  Vertex interpolateAboveOriginalCentroid() const { return mControlPoints[csControlIndexAboveOriginalCentroid]; }

  BezierIntersection intersect(Ray const &aRay, LimitPlaneIntersection const aShouldLimitPlaneIntersection) const;

private:
  Vector getNormal(Vector const &aBarycentric) const;
};

#endif
