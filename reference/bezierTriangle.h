#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER

#include "util.h"

template<typename tReal>
class BezierTriangle final {    // Cubic Bezier triangle
private:
  static constexpr uint32_t csControlIndexOriginalVertex0           = 0u;
  static constexpr uint32_t csControlIndexOriginalVertex1           = 1u;
  static constexpr uint32_t csControlIndexAboveOriginalCenter       = 2u;
  static constexpr uint32_t csControlIndexOnOriginalSide0           = 3u;
  static constexpr uint32_t csControlIndexOnOriginalSide1           = 4u;
  static constexpr uint32_t csControlIndexOnSideToOriginalCenter0   = 5u;
  static constexpr uint32_t csControlIndexOnSideToOriginalCenter1   = 6u;
  static constexpr uint32_t csControlIndexOnSideFromOriginalCenter0 = 7u;
  static constexpr uint32_t csControlIndexOnSideFromOriginalCenter1 = 8u;
  static constexpr uint32_t csControlIndexMiddle                    = 9u;

  static constexpr tReal    csProportionControlOnOriginalSide          = 0.333f;
  static constexpr tReal    csProportionControlOnOriginalVertexCentral = 0.333f;

  using Vector     = ::Vector<tReal>;
  using Vertex     = ::Vertex<tReal>;
  using Matrix     = ::Matrix<tReal>;
  using Triangle   = ::Triangle<tReal>;
  using Plane      = ::Plane<tReal>;

  Plane                    mOriginalPlane;
  std::array<Plane, 3u>    mNeighbourDividerPlanes;  // These are about the plane going through each edge and having the normal the average of the adjacent side normals.
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
  BezierTriangle(Vertex const &aOriginalCommonVertex0, Vertex const &aOriginalCommonVertex1, Vertex const &aOriginalCentral,
                 Vector const &aAverageNormal0, Vector const &aAverageNormal1, Plane const &aPlaneBetweenOriginalNeighbours,
                 std::array<uint32_t, 3u> const &aNeighbourIndices);

  void setMissingFields(Vertex const &aOriginalCentral, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious);
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tReal>
BezierTriangle<tReal>::BezierTriangle(Vertex const &aOriginalCommonVertex0, Vertex const &aOriginalCommonVertex1, Vertex const &aOriginalCentral,
                                      Vector const &aAverageNormal0, Vector const &aAverageNormal1, Plane const &aPlaneBetweenOriginalNeighbours,
                                      std::array<uint32_t, 3u> const &aNeighbourIndices) {
  mNeighbours = aNeighbourIndices;
  mControlPoints[csControlIndexOriginalVertex0] = aOriginalCommonVertex0;
  mControlPoints[csControlIndexOriginalVertex1] = aOriginalCommonVertex1;

  Plane commonPlaneVertex0{aAverageNormal0, aOriginalCommonVertex0.dot(aAverageNormal0) };
  Plane commonPlaneVertex1{aAverageNormal1, aOriginalCommonVertex0.dot(aAverageNormal1) };
  Plane perpendicularToOriginalSideInProportion0 = Plane::createFrom1proportion2points(csProportionControlOnOriginalSide, aOriginalCommonVertex0, aOriginalCommonVertex1);
  Plane perpendicularToOriginalSideInProportion1 = Plane::createFrom1proportion2points(csProportionControlOnOriginalSide, aOriginalCommonVertex1, aOriginalCommonVertex0);

  mControlPoints[csControlIndexOnOriginalSide0] = Plane::intersect(commonPlaneVertex0. aPlaneBetweenOriginalNeighbours, perpendicularToOriginalSideInProportion0);
  mControlPoints[csControlIndexOnOriginalSide1] = Plane::intersect(commonPlaneVertex1. aPlaneBetweenOriginalNeighbours, perpendicularToOriginalSideInProportion1);

  Vector originalNormal = getNormal(aOriginalCommonVertex0, aOriginalCommonVertex1, aOriginalCentral);
  Plane betweenSplitTriangles0 = Plane::createFrom1vector2points(originalNormal, aOriginalCommonVertex0, aOriginalCentral);
  Plane betweenSplitTriangles1 = Plane::createFrom1vector2points(originalNormal, aOriginalCommonVertex1, aOriginalCentral);
  Plane perpendicularToSplitBetweenTrianglesInProportion0 = Plane::createFrom1proportion2points(csProportionControlOnOriginalVertexCentral, aOriginalCommonVertex0, aOriginalCentral);
  Plane perpendicularToSplitBetweenTrianglesInProportion1 = Plane::createFrom1proportion2points(csProportionControlOnOriginalVertexCentral, aOriginalCommonVertex1, aOriginalCentral);

  mControlPoints[csControlIndexOnSideToOriginalCenter0] = Plane::intersect(commonPlaneVertex0. betweenSplitTriangles0, perpendicularToSplitBetweenTrianglesInProportion0);
  mControlPoints[csControlIndexOnSideToOriginalCenter1] = Plane::intersect(commonPlaneVertex1. betweenSplitTriangles1, perpendicularToSplitBetweenTrianglesInProportion1);

  Plane parallelToPlaneBetweenOriginalNeighbours; // TODO because of using commonPlaneVertex* this will not be parallel to aPlaneBetweenOriginalNeighbours
  Plane halfPlaneOfControlPointsInOriginalSide = Plane::createFrom1proportion2points(0.5f, mControlPoints[csControlIndexOnOriginalSide0], mControlPoints[csControlIndexOnOriginalSide1]);

  // legacy TODO remove
/*  mNormal = Mesh<tReal>::getNormal({aVertex0, aVertex1, aVertex2}).normalized();
  mNormalEquationConstant = mNormal.dot(aVertex0);
  Matrix vertices;
  vertices.column(0) = aVertex0; 
  vertices.column(1) = aVertex1; 
  vertices.column(2) = aVertex2; 
  mBarycentricInverse = vertices.inverse();*/
}

template<typename tReal>
void BezierTriangle<tReal>::setMissingFields(Vertex const &aOriginalCentral, BezierTriangle const &aTriangleNext, BezierTriangle const &aTrianglePrevious) {
}

#endif
