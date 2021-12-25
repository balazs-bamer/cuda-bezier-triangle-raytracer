#include "bezierTriangle.h"


BezierTriangle::BezierTriangle(Vertex const &aOriginalCommonVertex0, Vertex const &aOriginalCommonVertex1, Vertex const &aOriginalCentroid,
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

void BezierTriangle::setMissingFields1(Vertex const &aOriginalCentroid, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious) {
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

void BezierTriangle::setMissingFields2(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &) {
  mControlPoints[csControlIndexAboveOriginalCentroid] = (mControlPoints[csControlIndexOnSideFromOriginalCentroid0] +
                                                         mControlPoints[csControlIndexOnSideToOriginalCentroid1] +
                                                         aTriangleNext.mControlPoints[csControlIndexOnSideToOriginalCentroid1]) / 3.0f;

  mUnderlyingPlane = Plane::createFrom3points(mControlPoints[csControlIndexOriginalVertex0], mControlPoints[csControlIndexOriginalVertex1], mControlPoints[csControlIndexAboveOriginalCentroid]);

  mBarycentricInverse = getBarycentricInverse(mControlPoints[csControlIndexOriginalVertex0], mControlPoints[csControlIndexOriginalVertex1], mControlPoints[csControlIndexAboveOriginalCentroid]);

  mHeightInside = 0.0f;
  mHeightOutside = 0.0f;
  Triangle barycentric{Vertex{{1.0f, 0.0f, 0.0f}}, Vertex{{0.0f, 1.0f, 0.0f}}, Vertex{{0.0f, 0.0f, 1.0f}}};
  divide(barycentric, csHeightSampleDivisor, [this](Triangle && aNewBary) {
    for(uint32_t i = 0u; i < 3u; ++i) {
      auto distance = mUnderlyingPlane.distance(interpolate(aNewBary[i]));  // Will cont most points 2 or 3 times, but don't care now.
      mHeightInside = std::min(mHeightInside, distance);
      mHeightOutside = std::max(mHeightOutside, distance);
    }
  } );
  mHeightInside *= csHeightSafetyFactor;
  mHeightOutside *= csHeightSafetyFactor;
  mBezierDerivativeDirectionVectorA = Vertex{1.0f, 0.0f, -1.0f}; // Parallel to the side from 2 to 0, no matter how long.
  mBezierDerivativeDirectionVectorB = mBarycentricInverse *      // Perpendicular to the previous one.
      (mControlPoints[csControlIndexAboveOriginalCentroid] - mControlPoints[csControlIndexOriginalVertex0]).cross(mUnderlyingPlane.mNormal);
}

void BezierTriangle::setMissingFields3(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious) {
  mNeighbourDividerPlanes[1u] = Plane::createFrom1vector2points(mUnderlyingPlane.mNormal + aTriangleNext.mUnderlyingPlane.mNormal,
                                                                mControlPoints[csControlIndexOriginalVertex1],
                                                                mControlPoints[csControlIndexAboveOriginalCentroid]);
  mNeighbourDividerPlanes[2u] = Plane::createFrom1vector2points(mUnderlyingPlane.mNormal + aTrianglePrevious.mUnderlyingPlane.mNormal,
                                                                mControlPoints[csControlIndexOriginalVertex0],
                                                                mControlPoints[csControlIndexAboveOriginalCentroid]);
  mNeighbourDividerPlanes[1u].makeDistancePositive(mControlPoints[csControlIndexMiddle]);
  mNeighbourDividerPlanes[2u].makeDistancePositive(mControlPoints[csControlIndexMiddle]);
}

Vertex BezierTriangle::interpolateLinear(float const aBary0, float const aBary1, float const aBary2) const {
  return mControlPoints[csControlIndexOriginalVertex0]       * aBary0 +
         mControlPoints[csControlIndexOriginalVertex1]       * aBary1 +
         mControlPoints[csControlIndexAboveOriginalCentroid] * aBary2;
}

Vertex BezierTriangle::interpolate(float const aBary0, float const aBary1, float const aBary2) const {
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

BezierIntersection BezierTriangle::intersect(Ray const &aRay, LimitPlaneIntersection const aShouldLimitPlaneIntersection) const {
  auto inPlane = mUnderlyingPlane.intersect(aRay);
  BezierIntersection result;
  if(inPlane.mValid && ::abs(inPlane.mDistance) > -mHeightInside && ::abs(inPlane.mDistance) > mHeightOutside) {  // Make sure we don't intersect the same triangle again
    Vector barycentric = mBarycentricInverse * inPlane.mPoint;
    if(aShouldLimitPlaneIntersection == LimitPlaneIntersection::cNone ||
       (barycentric(0) >= 0.0f && barycentric(0) <= 1.0f &&
        barycentric(1) >= 0.0f && barycentric(1) <= 1.0f &&
        barycentric(2) >= 0.0f && barycentric(2) <= 1.0f)) {
      auto distanceInside = mHeightInside / inPlane.mCosIncidence;
      auto distanceOutside = mHeightOutside / inPlane.mCosIncidence;
      float closer = inPlane.mDistance + (inPlane.mCosIncidence > 0.0f ? distanceInside : distanceOutside);
      float further = inPlane.mDistance + (inPlane.mCosIncidence > 0.0f ? distanceOutside : distanceInside);

      Vertex pointOnRay = aRay.mStart + aRay.mDirection * closer;
      Vertex barycentric = mBarycentricInverse * mUnderlyingPlane.project(pointOnRay);
      auto diffCloser = ::abs(mUnderlyingPlane.distance(pointOnRay)) - ::abs(mUnderlyingPlane.distance(interpolate(barycentric)));

      pointOnRay = aRay.mStart + aRay.mDirection * further;
      barycentric = mBarycentricInverse * mUnderlyingPlane.project(pointOnRay);
      auto diffFurther = ::abs(mUnderlyingPlane.distance(pointOnRay)) - ::abs(mUnderlyingPlane.distance(interpolate(barycentric)));

      float middle;
      auto denom = diffCloser - diffFurther;
      if(::abs(denom) < csIntersectionEstimationEpsilon) {
        middle = (closer + further) / 2.0f;
      }
      else {
        middle = (diffCloser * further - diffFurther * closer) / denom;
      }
      Vector projectionDirection = mUnderlyingPlane.mNormal;

      for(uint32_t i = 0u; i < csRootSearchIterations; ++i) {
        result.mIntersection.mDistance = middle;
        pointOnRay = aRay.mStart + aRay.mDirection * middle;
        auto intersection = mUnderlyingPlane.intersect(pointOnRay, projectionDirection);
        result.mBarycentric = mBarycentricInverse * intersection.mPoint;
        result.mNormal = getNormal(result.mBarycentric);
        result.mIntersection.mPoint = interpolate(result.mBarycentric);
        projectionDirection = (result.mIntersection.mPoint - intersection.mPoint).normalized();
        middle = ((result.mIntersection.mPoint - aRay.mStart).dot(result.mNormal) / aRay.mDirection.dot(result.mNormal));
      }
      if(aRay.getDistance(result.mIntersection.mPoint) > csMaxIntersectionDistanceFromRay || result.mIntersection.mDistance < (further - closer) * csMinimalRayDistance) {
        result.mWhat = BezierIntersection::What::cNone;
      }
      else {
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
          result.mWhat = BezierIntersection::What::cFollowSide2;
        }
        else { // Most probably not possible for 2 sides at the same time. If yes, that rare case is not interesting
          result.mWhat = BezierIntersection::What::cIntersect;
          result.mIntersection.mCosIncidence = aRay.mDirection.dot(result.mNormal);
        }
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

Vector BezierTriangle::getNormal(Vector const &aBarycentric) const {
  float bary0_2 = aBarycentric(0) * aBarycentric(0);
  float bary1_2 = aBarycentric(1) * aBarycentric(1);
  float bary2_2 = aBarycentric(2) * aBarycentric(2);

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
                     (mControlPoints[csControlIndex111] * aBarycentric(0) * aBarycentric(2) +
                      mControlPoints[csControlIndex120] * aBarycentric(0) * aBarycentric(1) +
                      mControlPoints[csControlIndex021] * aBarycentric(1) * aBarycentric(2));

  Vector component2 = mControlPoints[csControlIndex003] * bary2_2 +
                      mControlPoints[csControlIndex201] * bary0_2 +
                      mControlPoints[csControlIndex021] * bary1_2 +
                      2.0f *
                     (mControlPoints[csControlIndex102] * aBarycentric(0) * aBarycentric(2) +
                      mControlPoints[csControlIndex012] * aBarycentric(1) * aBarycentric(2) +
                      mControlPoints[csControlIndex111] * aBarycentric(0) * aBarycentric(1));

  Vector componentA = mBezierDerivativeDirectionVectorA(0) * component0 +
                      mBezierDerivativeDirectionVectorA(1) * component1 +
                      mBezierDerivativeDirectionVectorA(2) * component2;
  Vector componentB = mBezierDerivativeDirectionVectorB(0) * component0 +
                      mBezierDerivativeDirectionVectorB(1) * component1 +
                      mBezierDerivativeDirectionVectorB(2) * component2;
  return componentA.cross(componentB).normalized();
}
