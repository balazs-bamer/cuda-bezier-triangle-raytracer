#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER

#include "util.h"

#include<iostream> // TODO remove

template<typename tReal>
struct BezierIntersection final {
  enum class What : uint8_t {
    cFollowSide0 = 0u,
    cFollowSide1 = 1u,
    cFollowSide2 = 2u,
    cNone        = 3u,
    cIntersect   = 4u
  };

  Intersection<tReal> mIntersection;
  Vertex<tReal>       mBarycentric;
  Vector<tReal>       mNormal;
  What                mWhat;
};

template<typename tReal>
class BezierTriangle final {    // Cubic Bezier triangle
public:
  static constexpr std::array<tReal, 3u>   csSampleRatiosOriginalSide = { 0.25f, 0.5f, 0.75f };

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

  static constexpr tReal    csProportionControlOnOriginalSide           = 0.291f;
  static constexpr tReal    csProportionControlOnOriginalVertexCentroid = 0.304f;
  static constexpr tReal    csProportionControlOnOriginalMedian         = 0.2f;
  static constexpr tReal    csHeightSafetyFactor                        = 1.33333333f;
  static constexpr tReal    csOneThird                                  = 1.0 / 3.0;
  static constexpr tReal    csRootSearchImpossibleFactor                = 0.03f;
  static constexpr tReal    csRootSearchEpsilonFactor                   = 0.001f;       // TODO consider

  using Vector             = ::Vector<tReal>;
  using Vertex             = ::Vertex<tReal>;
  using Matrix             = ::Matrix<tReal>;
  using Triangle           = ::Triangle<tReal>;
  using Ray                = ::Ray<tReal>;
  using Intersection       = ::Intersection<tReal>;
  using Plane              = ::Plane<tReal>;
  using BezierIntersection = ::BezierIntersection<tReal>;

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
  tReal                    mHeightInside;            // Sampled biggest distance of the Bezier triangle measured from the underlying triangle, < 0
  tReal                    mHeightOutside;           // this is > 0
  tReal                    mRootSearchEpsilon;
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

  Vertex interpolateLinear(tReal const aBary0, tReal const aBary1, tReal const aBary2) const;
  Vertex interpolateLinear(tReal const aBary0, tReal const aBary1) const { return interpolateLinear(aBary0, aBary1, 1.0f - aBary0 - aBary1); }
  Vertex interpolateLinear(Vertex const &aBary) const { return interpolateLinear(aBary(0u), aBary(1u), aBary(2u)); }

  Vertex interpolate(tReal const aBary0, tReal const aBary1, tReal const aBary2) const;
  Vertex interpolate(tReal const aBary0, tReal const aBary1) const { return interpolate(aBary0, aBary1, 1.0f - aBary0 - aBary1); }
  Vertex interpolate(Vertex const &aBary) const { return interpolate(aBary(0u), aBary(1u), aBary(2u)); }
  Vertex interpolateAboveOriginalCentroid() const { return mControlPoints[csControlIndexAboveOriginalCentroid]; }

  BezierIntersection intersect(Ray const &aRay, bool const aCheckBarycentricLimits) const;

private:
  BezierIntersection intersect(Ray const &aRay, tReal const aParameterCloser, tReal const aParameterFurther) const;
  std::tuple<tReal, Vertex, Vertex> getLineBezierDifferenceSignumBarySurface(Ray const &aRay, tReal const aDistance) const;
  Vector getNormal(Vector const &aBarycentric) const;
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tReal>
BezierTriangle<tReal>::BezierTriangle(Vertex const &aOriginalCommonVertex0, Vertex const &aOriginalCommonVertex1, Vertex const &aOriginalCentroid,
                                      Vector const &aAverageNormal0, Vector const &aAverageNormal1, Plane const &aPlaneBetweenOriginalNeighbours,
                                      std::array<uint32_t, 3u> const &aNeighbourIndices)
  : mNeighbours(aNeighbourIndices) {
  mControlPoints[csControlIndexOriginalVertex0] = aOriginalCommonVertex0;
  mControlPoints[csControlIndexOriginalVertex1] = aOriginalCommonVertex1;

  Plane const commonPlaneVertex0{aAverageNormal0, aOriginalCommonVertex0.dot(aAverageNormal0) };
  Plane const commonPlaneVertex1{aAverageNormal1, aOriginalCommonVertex1.dot(aAverageNormal1) };
  Plane const perpendicularToOriginalSideInProportion0 = Plane::createFrom1proportion2points(csProportionControlOnOriginalSide, aOriginalCommonVertex0, aOriginalCommonVertex1);
  Plane const perpendicularToOriginalSideInProportion1 = Plane::createFrom1proportion2points(csProportionControlOnOriginalSide, aOriginalCommonVertex1, aOriginalCommonVertex0);

  // For each control point obtained by planes intersection, the first plane is the critical one ensuring Bezier surface C1 continuity among triangles.
  mControlPoints[csControlIndexOnOriginalSide0] = Plane::intersect(commonPlaneVertex0, aPlaneBetweenOriginalNeighbours, perpendicularToOriginalSideInProportion0);
  mControlPoints[csControlIndexOnOriginalSide1] = Plane::intersect(commonPlaneVertex1, aPlaneBetweenOriginalNeighbours, perpendicularToOriginalSideInProportion1);

  Vector const originalNormal = ::getNormal(aOriginalCommonVertex0, aOriginalCommonVertex1, aOriginalCentroid);
  Plane const parallelToOriginalNormalBetweenSplitTriangles0 = Plane::createFrom1vector2points(originalNormal, aOriginalCommonVertex0, aOriginalCentroid);
  Plane const parallelToOriginalNormalBetweenSplitTriangles1 = Plane::createFrom1vector2points(originalNormal, aOriginalCommonVertex1, aOriginalCentroid);
  Plane const perpendicularToSplitBetweenTrianglesInProportion0 = Plane::createFrom1proportion2points(csProportionControlOnOriginalVertexCentroid, aOriginalCommonVertex0, aOriginalCentroid);
  Plane const perpendicularToSplitBetweenTrianglesInProportion1 = Plane::createFrom1proportion2points(csProportionControlOnOriginalVertexCentroid, aOriginalCommonVertex1, aOriginalCentroid);

  mControlPoints[csControlIndexOnSideFromOriginalCentroid1] = Plane::intersect(commonPlaneVertex0,
                                                                               parallelToOriginalNormalBetweenSplitTriangles0,
                                                                               perpendicularToSplitBetweenTrianglesInProportion0);
  mControlPoints[csControlIndexOnSideToOriginalCentroid0] = Plane::intersect(commonPlaneVertex1,
                                                                             parallelToOriginalNormalBetweenSplitTriangles1,
                                                                             perpendicularToSplitBetweenTrianglesInProportion1);

  Plane const perpendicularToPlaneBetweenOriginalNeighboursViaControlPointsInOriginalSide = Plane::createFrom1vector2points(aPlaneBetweenOriginalNeighbours.mNormal, mControlPoints[csControlIndexOnOriginalSide0], mControlPoints[csControlIndexOnOriginalSide1]);
  Plane const halfPlaneOfControlPointsInOriginalSide = Plane::createFrom1proportion2points(0.5f, mControlPoints[csControlIndexOnOriginalSide0], mControlPoints[csControlIndexOnOriginalSide1]);
  Plane const perpendicularToOriginalMedianInProportion = Plane::createFrom1proportion2points(csProportionControlOnOriginalMedian, (aOriginalCommonVertex0 + aOriginalCommonVertex1) / 2.0f, aOriginalCentroid);

  mControlPoints[csControlIndexMiddle] = Plane::intersect(perpendicularToPlaneBetweenOriginalNeighboursViaControlPointsInOriginalSide,
                                                          halfPlaneOfControlPointsInOriginalSide,
                                                          perpendicularToOriginalMedianInProportion);

  mNeighbourDividerPlanes[0u] = aPlaneBetweenOriginalNeighbours;
  mNeighbourDividerPlanes[0u].makeDistancePositive(mControlPoints[csControlIndexMiddle]);
}

template<typename tReal>
void BezierTriangle<tReal>::setMissingFields1(Vertex const &aOriginalCentroid, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious) {
  Vector const originalNormal = ::getNormal(mControlPoints[csControlIndexOriginalVertex0], mControlPoints[csControlIndexOriginalVertex1], aOriginalCentroid);
  Plane const twoMiddlesAndSplitCloseToOriginalCommonVertex0 = Plane::createFrom3points(mControlPoints[csControlIndexOnSideFromOriginalCentroid1], mControlPoints[csControlIndexMiddle], aTrianglePrevious.mControlPoints[csControlIndexMiddle]);
  Plane const twoMiddlesAndSplitCloseToOriginalCommonVertex1 = Plane::createFrom3points(mControlPoints[csControlIndexOnSideToOriginalCentroid0], aTriangleNext.mControlPoints[csControlIndexMiddle], mControlPoints[csControlIndexMiddle]);
  Plane const parallelToOriginalNormalBetweenSplitTriangles0 = Plane::createFrom1vector2points(originalNormal, mControlPoints[csControlIndexOriginalVertex0], aOriginalCentroid);
  Plane const parallelToOriginalNormalBetweenSplitTriangles1 = Plane::createFrom1vector2points(originalNormal, mControlPoints[csControlIndexOriginalVertex1], aOriginalCentroid);
  Plane const perpendicularToSplitBetweenTrianglesInProportion0 = Plane::createFrom1proportion2points(csProportionControlOnOriginalVertexCentroid, aOriginalCentroid, mControlPoints[csControlIndexOriginalVertex0]);
  Plane const perpendicularToSplitBetweenTrianglesInProportion1 = Plane::createFrom1proportion2points(csProportionControlOnOriginalVertexCentroid, aOriginalCentroid, mControlPoints[csControlIndexOriginalVertex1]);

  mControlPoints[csControlIndexOnSideFromOriginalCentroid0] = Plane::intersect(twoMiddlesAndSplitCloseToOriginalCommonVertex0,
                                                                               parallelToOriginalNormalBetweenSplitTriangles0,
                                                                               perpendicularToSplitBetweenTrianglesInProportion0);
  mControlPoints[csControlIndexOnSideToOriginalCentroid1] = Plane::intersect(twoMiddlesAndSplitCloseToOriginalCommonVertex1,
                                                                             parallelToOriginalNormalBetweenSplitTriangles1,
                                                                             perpendicularToSplitBetweenTrianglesInProportion1);
}

template<typename tReal>
void BezierTriangle<tReal>::setMissingFields2(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &) {
  mControlPoints[csControlIndexAboveOriginalCentroid] = (mControlPoints[csControlIndexOnSideFromOriginalCentroid0] +
                                                         mControlPoints[csControlIndexOnSideToOriginalCentroid1] +
                                                         aTriangleNext.mControlPoints[csControlIndexOnSideToOriginalCentroid1]) / 3.0f;

  mUnderlyingPlane = Plane::createFrom3points(mControlPoints[csControlIndexOriginalVertex0], mControlPoints[csControlIndexOriginalVertex1], mControlPoints[csControlIndexAboveOriginalCentroid]);

  auto const &column0 = mControlPoints[csControlIndexOriginalVertex0];
  auto const &column1 = mControlPoints[csControlIndexOriginalVertex1];
  auto const &column2 = mControlPoints[csControlIndexAboveOriginalCentroid];
  Matrix vertices {
    { column0(0), column1(0), column2(0) },
    { column0(1), column1(1), column2(1) },
    { column0(2), column1(2), column2(2) },
  };
  mBarycentricInverse = vertices.inverse();

  mRootSearchEpsilon = csRootSearchEpsilonFactor * ((column0 - column1).norm() + (column1 - column2).norm() + (column2 - column0).norm());

  Vertex bary{csOneThird, csOneThird, csOneThird};
  tReal distance = mUnderlyingPlane.distance(interpolate(bary));
  mHeightInside = std::min(0.0f, distance);
  mHeightOutside = std::max(0.0f, distance);
  for(uint32_t i = 0u; i < 3u; ++i) {
    for(auto const ratio : BezierTriangle::csSampleRatiosOriginalSide) {
      bary(i) = ratio;
      bary((i + 1u) % 3u) = 1.0f - ratio;
      bary((i + 2u) % 3u) = 0.0f;
      distance = mUnderlyingPlane.distance(interpolate(bary));
      mHeightInside = std::min(0.0f, distance);
      mHeightOutside = std::max(0.0f, distance);
    }
  }

  mBezierDerivativeDirectionVectorA = Vertex{1.0f, 0.0f, -1.0f};  // No matter how long.
  mBezierDerivativeDirectionVectorB = mBarycentricInverse *
      (mControlPoints[csControlIndexAboveOriginalCentroid] - mControlPoints[csControlIndexOriginalVertex0]).cross(mUnderlyingPlane.mNormal);
}

template<typename tReal>
void BezierTriangle<tReal>::setMissingFields3(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious) {
  mNeighbourDividerPlanes[1u] = Plane::createFrom1vector2points(mUnderlyingPlane.mNormal + aTriangleNext.mUnderlyingPlane.mNormal,
                                                                mControlPoints[csControlIndexOriginalVertex1],
                                                                mControlPoints[csControlIndexAboveOriginalCentroid]);
  mNeighbourDividerPlanes[2u] = Plane::createFrom1vector2points(mUnderlyingPlane.mNormal + aTrianglePrevious.mUnderlyingPlane.mNormal,
                                                                mControlPoints[csControlIndexOriginalVertex0],
                                                                mControlPoints[csControlIndexAboveOriginalCentroid]);
  mNeighbourDividerPlanes[1u].makeDistancePositive(mControlPoints[csControlIndexMiddle]);
  mNeighbourDividerPlanes[2u].makeDistancePositive(mControlPoints[csControlIndexMiddle]);
}

template<typename tReal>
Vertex<tReal> BezierTriangle<tReal>::interpolateLinear(tReal const aBary0, tReal const aBary1, tReal const aBary2) const {
  return mControlPoints[csControlIndexOriginalVertex0]       * aBary0 +
         mControlPoints[csControlIndexOriginalVertex1]       * aBary1 +
         mControlPoints[csControlIndexAboveOriginalCentroid] * aBary2;
}

template<typename tReal>
Vertex<tReal> BezierTriangle<tReal>::interpolate(tReal const aBary0, tReal const aBary1, tReal const aBary2) const {
  auto const bary0_2 = aBary0 * aBary0;
  auto const bary1_2 = aBary1 * aBary1;
  auto const bary2_2 = aBary2 * aBary2;

  return mControlPoints[csControlIndexOriginalVertex0]             * aBary0 * bary0_2 +
         mControlPoints[csControlIndexOriginalVertex1]             * aBary1 * bary1_2 +
         mControlPoints[csControlIndexAboveOriginalCentroid]       * aBary2 * bary2_2 +
         3.0f *
        (mControlPoints[csControlIndexOnOriginalSide0]             * aBary1 * bary0_2 +
         mControlPoints[csControlIndexOnOriginalSide1]             * aBary0 * bary1_2 +
         mControlPoints[csControlIndexOnSideToOriginalCentroid0]   * aBary2 * bary1_2 +
         mControlPoints[csControlIndexOnSideToOriginalCentroid1]   * aBary1 * bary2_2 +
         mControlPoints[csControlIndexOnSideFromOriginalCentroid0] * aBary0 * bary2_2 +
         mControlPoints[csControlIndexOnSideFromOriginalCentroid1] * aBary2 * bary0_2) +
         mControlPoints[csControlIndexMiddle]                      * aBary0 * aBary1 * aBary2 * 6.0f;
}

template<typename tReal>
BezierIntersection<tReal> BezierTriangle<tReal>::intersect(Ray const &aRay, bool const aCheckBarycentricLimits) const {
  auto inPlane = mUnderlyingPlane.intersect(aRay);
  BezierIntersection result;
  if(inPlane.mValid) {
    Vector barycentric = mBarycentricInverse * inPlane.mPoint;
    if(!aCheckBarycentricLimits ||
       (barycentric(0) >= 0.0f && barycentric(0) <= 1.0f &&
        barycentric(1) >= 0.0f && barycentric(1) <= 1.0f &&
        barycentric(2) >= 0.0f && barycentric(2) <= 1.0f)) {
      auto distanceInside = mHeightInside / inPlane.mCosIncidence;
      auto distanceOutside = mHeightOutside / inPlane.mCosIncidence;
      tReal parameterCloser = inPlane.mDistance + (inPlane.mCosIncidence > 0.0f ? distanceInside : distanceOutside);
      tReal parameterFurther = inPlane.mDistance + (inPlane.mCosIncidence > 0.0f ? distanceOutside : distanceInside);
      auto totalInterestingRange = parameterFurther - parameterCloser;
      if(::abs(inPlane.mDistance - parameterCloser) / totalInterestingRange > csRootSearchImpossibleFactor) {    // Worth to search the closer half first
        result = intersect(aRay, parameterCloser, inPlane.mDistance);
      }
      else { // nothing to do
      }
      if(result.mWhat == BezierIntersection::What::cNone && ::abs(parameterFurther - inPlane.mDistance) / totalInterestingRange > csRootSearchImpossibleFactor) {
        result = intersect(aRay, inPlane.mDistance, parameterFurther);
      }
      else { // nothing to do
      }
      if(result.mWhat == BezierIntersection::What::cIntersect) {
        result.mNormal = getNormal(result.mBarycentric);
        result.mIntersection.mCosIncidence = aRay.mDirection.dot(result.mNormal); // TODO find out where it should point
      }
      else { // Nothing to do
      }
    }
    else {
      result.mWhat = BezierIntersection::What::cNone;
    }
  }
  else {
    result.mWhat = BezierIntersection::What::cNone;
  }
  return result;
}

template<typename tReal>
BezierIntersection<tReal> BezierTriangle<tReal>::intersect(Ray const &aRay, tReal const aParameterCloser, tReal const aParameterFurther) const {
  BezierIntersection result;
  auto closer = aParameterCloser;
  auto further = aParameterFurther;
  auto [signumCloser, barycentricCloser, surfaceCloser] = getLineBezierDifferenceSignum(aRay, closer);
  tReal signumFurther;
  std::tie(signumFurther, result.mBarycentric, result.mIntersection.mPoint) = getLineBezierDifferenceSignum(aRay, further);
  result.mIntersection.mDistance = further;
  if(signumCloser == signumFurther) {
    result.mWhat = BezierIntersection::What::cNone;
  }
  else {
    while(further - closer > mRootSearchEpsilon) {
      auto middle = (closer + further) / 2.0f;
      tReal signumMiddle;
      std::tie(signumMiddle, result.mBarycentric, result.mIntersection.mPoint) = getLineBezierDifferenceSignum(aRay, middle);
      result.mIntersection.mDistance = middle;
      if(signumCloser == signumMiddle) {
        closer = middle;
      }
      else {
        further = middle;
      }
    }

    uint32_t outside = (mNeighbourDividerPlanes[0].distance(result.mIntersection.mPoint) < 0.0f ? 1u : 0u);
    outside |= (mNeighbourDividerPlanes[1].distance(result.mIntersection.mPoint) < 0.0f ? 2u : 0u);
    outside |= (mNeighbourDividerPlanes[2].distance(result.mIntersection.mPoint) < 0.0f ? 4u : 0u);
    if(outside == 1u) {
      result.mWhat = BezierIntersection::What::cFollowSide0;
    }
    else if(outside == 2u) {
      result.mWhat = BezierIntersection::What::cFollowSide1;
    }
    else if(outside == 4u) {
      result.mWhat = BezierIntersection::What::cFollowSide0;
    }
    else { // Most probably not possible for 2 sides at the same time. If yes, that rare case is not interesting
      result.mWhat = BezierIntersection::What::cIntersect;
    }
  }
  return result;
}

template<typename tReal>
std::tuple<tReal, Vertex<tReal>, Vertex<tReal>> BezierTriangle<tReal>::getLineBezierDifferenceSignumBarySurface(Ray const &aRay, tReal const aDistance) const {
  Vertex pointOnRay = aRay.mStart + aRay.mDirection * aDistance;
  Vertex barycentric = mBarycentricInverse * (mUnderlyingPlane.project(pointOnRay));
  Vertex pointOnSurface = interpolate(barycentric);
  return std::make_tuple(::copysign(1.0f, ::abs(mUnderlyingPlane.distance(pointOnRay)) - ::abs(mUnderlyingPlane.distance(pointOnSurface))), barycentric, pointOnSurface);
}

template<typename tReal>
Vector<tReal> BezierTriangle<tReal>::getNormal(Vector const &aBarycentric) const {
  tReal bary0_2 = aBarycentric(0) * aBarycentric(0);
  tReal bary1_2 = aBarycentric(1) * aBarycentric(1);
  tReal bary2_2 = aBarycentric(2) * aBarycentric(2);

  Vector component0 = mControlPoints[csControlIndex300] * bary0_2 +
                      mControlPoints[csControlIndex102] * bary2_2 +
                      mControlPoints[csControlIndex120] * bary1_2 +
                      2.0f *
                     (mControlPoints[csControlIndex201] * aBarycentric(0) * aBarycentric(2) +
                      mControlPoints[csControlIndex210] * aBarycentric(0) * aBarycentric(1) +
                      mControlPoints[csControlIndex111] * aBarycentric(2) * aBarycentric(1));

  Vector component1 = mControlPoints[csControlIndex030] * bary1_2 +
                      mControlPoints[csControlIndex012] * bary2_2 +
                      mControlPoints[csControlIndex210] * bary0_2 +
                      2.0f *
                     (mControlPoints[csControlIndex102] * aBarycentric(0) * aBarycentric(2) +
                      mControlPoints[csControlIndex120] * aBarycentric(0) * aBarycentric(1) +
                      mControlPoints[csControlIndex111] * aBarycentric(2) * aBarycentric(0));

  Vector component2 = mControlPoints[csControlIndex003] * bary2_2 +
                      mControlPoints[csControlIndex201] * bary0_2 +
                      mControlPoints[csControlIndex021] * bary1_2 +
                      2.0f *
                     (mControlPoints[csControlIndex102] * aBarycentric(0) * aBarycentric(2) +
                      mControlPoints[csControlIndex012] * aBarycentric(1) * aBarycentric(1) +
                      mControlPoints[csControlIndex111] * aBarycentric(0) * aBarycentric(1));

  Vector componentA = mBezierDerivativeDirectionVectorA(0) * component0 +
                      mBezierDerivativeDirectionVectorA(1) * component1 +
                      mBezierDerivativeDirectionVectorA(2) * component2;
  Vector componentB = mBezierDerivativeDirectionVectorB(0) * component0 +
                      mBezierDerivativeDirectionVectorB(1) * component1 +
                      mBezierDerivativeDirectionVectorB(2) * component2;
  return componentA.cross(componentB).normalized();
}
#endif
