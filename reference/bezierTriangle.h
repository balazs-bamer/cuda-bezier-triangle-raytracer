#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER

#include "util.h"

template<typename tReal>
class BezierTriangle final {    // Cubic Bezier triangle
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

  static constexpr tReal    csProportionControlOnOriginalSide           = 0.333f;
  static constexpr tReal    csProportionControlOnOriginalVertexCentroid = 0.333f;
  static constexpr tReal    csProportionControlOnOriginalMedian         = 0.333f;

  using Vector     = ::Vector<tReal>;
  using Vertex     = ::Vertex<tReal>;
  using Matrix     = ::Matrix<tReal>;
  using Triangle   = ::Triangle<tReal>;
  using Plane      = ::Plane<tReal>;

  Plane                    mUnderlyingPlane;         // Plane span by mControlPoints[0-2]
  std::array<Plane, 3u>    mNeighbourDividerPlanes;  // For indices 1 and 2, these are about the plane going through each edge in direction the average of the adjacent side normals.
                                                     // Index 0 will be aPlaneBetweenOriginalNeighbours
  std::array<uint32_t, 3u> mNeighbours;              // With new indices after Clough-Tocher split.
  // 0-1 identical to original triangle vertices, named aOriginalCommonVertex0, aOriginalCommonVertex1 in constructor
  // 2 somewhere above the middle of the original triangle center, to be calculated
  // 3-8 somewhere above sides of triangle (0-2): 3i+3, 3i+4 on the same side as (i, (i+1)%3), to be calculated
  // 9 in the middle, to be calculated
  std::array<Vertex, 10u>  mControlPoints;
  // Perhaps not needed, because the iterative method to find ray and Bezier surface intersection takes long. std::array<Vertex, 12u> mDerivativeControlPoints; // T most probably more than 2*6 for each partial derivative.
  Matrix                   mBarycentricInverse;      // T = (v1 v2 v3), b = barycentric coefficient column, v = point on plane, v=Tb, b = mBI * v
                                                     // multiplication by the inverse proven to be much faster than solving the linear equation under Eigen.

public:
  // Vertices in argument are in order such that the normal points to the desired direction.
  // Neighbour i is neighbour of edge (i, i + 1)
  BezierTriangle(Vertex const &aOriginalCommonVertex0, Vertex const &aOriginalCommonVertex1, Vertex const &aOriginalCentroid,
                 Vector const &aAverageNormal0, Vector const &aAverageNormal1, Plane const &aPlaneBetweenOriginalNeighbours,
                 std::array<uint32_t, 3u> const &aNeighbourIndices);

  void setMissingFields1(Vertex const &aOriginalCentroid, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious);
  void setMissingFields2(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &);
  void setMissingFields3(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious);

  Vertex interpolate(tReal const aBary0, tReal const aBary1, tReal const aBary2) const;
  Vertex interpolate(tReal const aBary0, tReal const aBary1) const { return resolveParameters(aBary0, aBary1, 1.0f - aBary0 - aBary1); }


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
  mControlPoints[csControlIndexOnOriginalSide0] = Plane::intersect(commonPlaneVertex0. aPlaneBetweenOriginalNeighbours, perpendicularToOriginalSideInProportion0);
  mControlPoints[csControlIndexOnOriginalSide1] = Plane::intersect(commonPlaneVertex1. aPlaneBetweenOriginalNeighbours, perpendicularToOriginalSideInProportion1);

  Vector const originalNormal = getNormal(aOriginalCommonVertex0, aOriginalCommonVertex1, aOriginalCentroid);
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
  Plane const perpendicularToOriginalMedianInProportion = Plane::createFrom1proportion2points(csProportionControlOnOriginalMedian, (aOriginalCommonVertex0, aOriginalCommonVertex1) / 2.0f, aOriginalCentroid);

  mControlPoints[csControlIndexMiddle] = Plane::intersect(perpendicularToPlaneBetweenOriginalNeighboursViaControlPointsInOriginalSide,
                                                          halfPlaneOfControlPointsInOriginalSide,
                                                          perpendicularToOriginalMedianInProportion);

  mNeighbourDividerPlanes[0u] = aPlaneBetweenOriginalNeighbours;
}

template<typename tReal>
void BezierTriangle<tReal>::setMissingFields1(Vertex const &aOriginalCentroid, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious) {
  Vector const originalNormal = getNormal(mControlPoints[csControlIndexOriginalVertex0], mControlPoints[csControlIndexOriginalVertex1], aOriginalCentroid);
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

  Matrix vertices;
  vertices.column(0) =  mControlPoints[csControlIndexOriginalVertex0];
  vertices.column(1) =  mControlPoints[csControlIndexOriginalVertex1];
  vertices.column(2) =  mControlPoints[csControlIndexAboveOriginalCentroid];
  mBarycentricInverse = vertices.inverse();
}

template<typename tReal>
void BezierTriangle<tReal>::setMissingFields3(Vertex const &, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious) {
  mNeighbourDividerPlanes[1u] = Plane::createFrom1vector2points(mUnderlyingPlane.mNormal + aTriangleNext.mUnderlyingPlane.mNormal,
                                                                mControlPoints[csControlIndexOriginalVertex1],
                                                                mControlPoints[csControlIndexAboveOriginalCentroid]);
  mNeighbourDividerPlanes[2u] = Plane::createFrom1vector2points(mUnderlyingPlane.mNormal + aTrianglePrevious.mUnderlyingPlane.mNormal,
                                                                mControlPoints[csControlIndexOriginalVertex0],
                                                                mControlPoints[csControlIndexAboveOriginalCentroid]);
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
#endif
