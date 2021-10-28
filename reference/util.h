#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_UTIL
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_UTIL

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t

#include <Eigen/Dense>

template<typename tReal>
constexpr tReal cgPi = 3.14159265358979323846;

template<typename tReal>
using Vector         = Eigen::Matrix<tReal, 3, 1>;

template<typename tReal>
using Vertex         = Eigen::Matrix<tReal, 3, 1>;

template<typename tReal>
using Matrix         = Eigen::Matrix<tReal, 3, 3>;

template<typename tReal>
using Transform      = Eigen::Matrix<tReal, 3, 3>;

template<typename tReal>
using Triangle       = std::array<Vertex<tReal>, 3u>;

template<typename tReal>
Vector<tReal> getNormal(Triangle<tReal> const &aFace) {
  return (aFace[1] - aFace[0]).cross(aFace[2] - aFace[0]);
}

template<typename tReal>
Vector<tReal> getNormal(Vertex<tReal> const &aVertex0, Vertex<tReal> const &aVertex1, Vertex<tReal> const &aVertex2) {
  return (aVertex1 - aVertex0).cross(aVertex2 - aVertex0);
}

// Assumes sum(aB*) == 1
template<typename tReal>
Vertex<tReal> barycentric2cartesian(Triangle<tReal> const &aTriangle, tReal const aB0, tReal const aB1, tReal const aB2) {
  return aTriangle[0u] * aB0 + aTriangle[1u] * aB1 + aTriangle[2u] * aB2;
}

template<typename tReal>
Vertex<tReal> barycentric2cartesian(Triangle<tReal> const &aTriangle, tReal const aB0, tReal const aB1) {
  return aTriangle[0u] * aB0 + aTriangle[1u] * aB1 + aTriangle[2u] * (1.0f - aB0 - aB1);
}

// Assumes sum(aB*) == 1
template<typename tReal>
Vertex<tReal> barycentric2cartesian(Vertex<tReal> const &aV0, Vertex<tReal> const &aV1, Vertex<tReal> const &aV2, tReal const aB0, tReal const aB1, tReal const aB2) {
  return aV0 * aB0 + aV1 * aB1 + aV2 * aB2;
}

template<typename tReal>
Vertex<tReal> barycentric2cartesian(Vertex<tReal> const &aV0, Vertex<tReal> const &aV1, Vertex<tReal> const &aV2, tReal const aB0, tReal const aB1) {
  return aV0 * aB0 + aV1 * aB1 + aV2 * (1.0f - aB0 - aB1);
}

template<typename tReal>
struct Plane final {
  Vector<tReal> mNormal;
  tReal         mConstant;

  Plane() = default;
  Plane(Vector<tReal> const &aNormal, tReal const aConstant) : mNormal(aNormal), mConstant(aConstant) {}

  static Plane         createFrom1proportion2points(tReal const aProportion, Vertex<tReal> const &aPoint0, Vertex<tReal> const &aPoint1);
  static Plane         createFrom3points(Vertex<tReal> const &aPoint0, Vertex<tReal> const &aPoint1, Vertex<tReal> const &aPoint2);
  static Plane         createFrom1vector2points(Vector<tReal> const &aDirection, Vertex<tReal> const &aPoint0, Vertex<tReal> const &aPoint1);
  static Plane         createFrom2vectors1point(Vertex<tReal> const &aDirection0, Vertex<tReal> const &aDirection1, Vertex<tReal> const &aPoint);
  static Vertex<tReal> intersect(Plane const &aPlane0, Plane const &aPlane1, Plane const &aPlane2);
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tReal, typename tLambda>
void divide(Triangle<tReal> const &aTriangle, int32_t const aDivisor, tLambda &&aCollector) {
  auto vector01 = (aTriangle[1] - aTriangle[0]) / aDivisor;
  auto vector02 = (aTriangle[2] - aTriangle[0]) / aDivisor;
  auto lineBase = aTriangle[0];
  auto base0 = lineBase;
  auto base1 = (aDivisor > 1) ? (base0 + vector01) : aTriangle[1];
  auto base2 = (aDivisor > 1) ? (base0 + vector02) : aTriangle[2];
  for(int32_t i = 0; i < aDivisor - 1; ++i) {
    for(int32_t j = 0; j < aDivisor - i - 1; j++) {
      aCollector({base0, base1, base2});
      auto base1next = base1 + vector02;
      aCollector({base1, base1next, base2});
      base1 = base1next;
      base0 = base2;
      base2 += vector02;
    }
    aCollector({base0, base1, base2});
    lineBase += vector01;
    base0 = lineBase;
    base1 = base0 + vector01;
    base2 = base0 + vector02;
  }
  aCollector({base0, aTriangle[1], base2});
}

// If aProportion < 0.5, the result will be closer to aPoint0
template<typename tReal>
Plane<tReal> Plane<tReal>::createFrom1proportion2points(tReal const aProportion, Vertex<tReal> const &aPoint0, Vertex<tReal> const &aPoint1) {
  Plane result;
  result.mNormal = (aPoint1 - aPoint0).normalized();
  result.mConstant = result.mNormal.dot(aPoint1 * aProportion + aPoint0 * (1.0f - aProportion));
  return result;
}

template<typename tReal>
Plane<tReal> Plane<tReal>::createFrom3points(Vertex<tReal> const &aPoint0, Vertex<tReal> const &aPoint1, Vertex<tReal> const &aPoint2) {
  Plane result;
  result.mNormal = (aPoint1 - aPoint0).cross(aPoint2 - aPoint0).normalized();
  result.mConstant = result.mNormal.dot(aPoint0);
  return result;
}

template<typename tReal>
Plane<tReal> Plane<tReal>::createFrom1vector2points(Vector<tReal> const &aDirection, Vertex<tReal> const &aPoint0, Vertex<tReal> const &aPoint1) {
  Plane result;
  result.mNormal = aDirection.cross(aPoint1 - aPoint0).normalized();
  result.mConstant = result.mNormal.dot(aPoint0);
  return result;
}

template<typename tReal>
Plane<tReal> Plane<tReal>::createFrom2vectors1point(Vertex<tReal> const &aDirection0, Vertex<tReal> const &aDirection1, Vertex<tReal> const &aPoint) {
  Plane result;
  result.mNormal = aDirection0.cross(aDirection1).normalized();
  result.mConstant = result.mNormal.dot(aPoint);
  return result;
}

template<typename tReal>
Vertex<tReal> Plane<tReal>::intersect(Plane const &aPlane0, Plane const &aPlane1, Plane const &aPlane2) {
  Matrix<tReal> matrix {
    { aPlane0.mNormal(0), aPlane0.mNormal(1), aPlane0.mNormal(2) },
    { aPlane1.mNormal(0), aPlane1.mNormal(1), aPlane1.mNormal(2) },
    { aPlane2.mNormal(0), aPlane2.mNormal(1), aPlane2.mNormal(2) }
  };
  Vector<tReal> vector{ aPlane0.mConstant, aPlane1.mConstant, aPlane2.mConstant };
  Vertex<tReal> result = matrix.inverse() * vector;
  // TODO remove
  auto error =::abs((matrix * result - vector).norm() / vector.norm());
  static_assert(std::is_same_v<decltype(error), tReal>);
  if(error > 0.1f) {
    throw "Plane::intersect rel error";
  }
  return result;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif
