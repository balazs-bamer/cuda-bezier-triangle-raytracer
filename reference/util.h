#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_UTIL
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_UTIL

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t

#include <Eigen/Dense>

template<typename tReal>
constexpr tReal cgPi = 3.14159265358979323846;

template<typename tReal>
constexpr tReal cgGeneralEpsilon = 1.0e-5;

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

template<typename tReal>
Vector<tReal> getAperpendicular(Vector<tReal> const &aVector);

template<typename tReal>
tReal getPerimeter(Triangle<tReal> const &aTriangle) {
  return (aTriangle[0u] - aTriangle[1u]).norm() + (aTriangle[1u] - aTriangle[2u]).norm() + (aTriangle[2u] - aTriangle[0u]).norm();
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
Matrix<tReal> getBarycentricInverse(Vertex<tReal> const &aVertex0, Vertex<tReal> const &aVertex1, Vertex<tReal> const &aVertex2) {
  Matrix<tReal> vertices {
    { aVertex0(0), aVertex1(0), aVertex2(0) },
    { aVertex0(1), aVertex1(1), aVertex2(1) },
    { aVertex0(2), aVertex1(2), aVertex2(2) },
  };
  return vertices.inverse();
}

template<typename tReal>
Vector<tReal> getAltitude(Vertex<tReal> const &aCommon1, Vertex<tReal> const &aCommon2, Vertex<tReal> const &aIndependent);

// Takes barycentric ends of a vector with aStart inside the triangle, finds out which side it will intersect.
// 0 between 300 and 030
// 1 between 030 and 003
// 2 between 003 and 300
template<typename tReal>
uint32_t toWhichSide(Vertex<tReal> const &aStart, Vertex<tReal> const &aEnd);

template<typename tReal>
struct Ray final {
  Vertex<tReal> mStart;
  Vector<tReal> mDirection; // normalized

  Ray() = default;
  Ray(Vertex<tReal> const &aStart, Vector<tReal> const &aDirection)
  : mStart(aStart)
  , mDirection(aDirection.normalized()) {}

  tReal getAverageErrorSquared(std::vector<Vertex<tReal>> const &aPoints) const;
};

template<typename tReal>
struct Intersection final {
  bool          mValid;         // Will be true even if distance < 0, because it can be important.
  Vertex<tReal> mPoint;
  tReal         mCosIncidence;  // Negative if the ray comes from outside, so the surface normal and the ray point in opposite direction.
  tReal         mDistance;      // Distance of ray source point and intersection point.
};

// Plane equation is in the form of point.dot(mNormal) == mConstant where point is any point in the plane.
template<typename tReal>
struct Plane final {
  static constexpr tReal csRayPlaneIntersectionEpsilon = 0.00001f;

  Vector<tReal> mNormal;        // normalized
  tReal         mConstant;

  Plane() = default;
  // aNormal must be normalized.
  Plane(Vector<tReal> const &aNormal, tReal const aConstant) : mNormal(aNormal), mConstant(aConstant) {}

  static Plane         createFrom1proportion2points(tReal const aProportion, Vertex<tReal> const &aPoint0, Vertex<tReal> const &aPoint1);
  static Plane         createFrom3points(Vertex<tReal> const &aPoint0, Vertex<tReal> const &aPoint1, Vertex<tReal> const &aPoint2);
  static Plane         createFromTriangle(Triangle<tReal> const &aTriangle) { return createFrom3points(aTriangle[0u], aTriangle[1u], aTriangle[2u]); }
  static Plane         createFrom1vector2points(Vector<tReal> const &aDirection, Vertex<tReal> const &aPoint0, Vertex<tReal> const &aPoint1);
  static Plane         createFrom2vectors1point(Vertex<tReal> const &aDirection0, Vertex<tReal> const &aDirection1, Vertex<tReal> const &aPoint);
  static Vertex<tReal> intersect(Plane const &aPlane0, Plane const &aPlane1, Plane const &aPlane2);
  Intersection<tReal>  intersect(Ray<tReal> const &aRay) const;

  // Point projection on this plane
  Vector<tReal>        project(Vector<tReal> const &aPoint) const { return aPoint - mNormal * (aPoint.dot(mNormal) - mConstant); }

  // Distance of point and this plane, >0 if the normal points towards the point.
  tReal                distance(Vector<tReal> const &aPoint) const { return aPoint.dot(mNormal) - mConstant; }

  void                 makeDistancePositive(Vector<tReal> const aPoint);
  void                 makeDistanceNegative(Vector<tReal> const aPoint);

  bool operator<(Plane<tReal> const &aOther) const { return mNormal(0) < aOther.mNormal(0) // TODO remove: only for debugging
                                                         || mNormal(0) == aOther.mNormal(0) && mNormal(1) < aOther.mNormal(1)
                                                         || mNormal(0) == aOther.mNormal(0) && mNormal(1) == aOther.mNormal(1) && mNormal(2) < aOther.mNormal(2)
                                                         || mNormal(0) == aOther.mNormal(0) && mNormal(1) == aOther.mNormal(1) && mNormal(2) == aOther.mNormal(2) && mConstant < aOther.mConstant; }
};

template<typename tReal>
struct Spherical final {
  tReal mR;
  tReal mAzimuth; // or sector
  tReal mInclination;  // or belt

  // Assumes the conversion can be done.
  Spherical(tReal const aX, tReal const aY, tReal const aZ)
  : mR(::sqrt(aX * aX + aY * aY + aZ * aZ))
  , mInclination(::acos(aZ / mR))
  , mAzimuth(::atan2(aY, aX)) {}
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tReal>
Vector<tReal> getAperpendicular(Vector<tReal> const &aVector) {
  static constexpr tReal csEpsilon = 1e-10;

  Vector<tReal> result;
  result(0) = 0.0f;
  if(::abs(aVector(1)) < csEpsilon && ::abs(aVector(2)) < csEpsilon) {
    result(1) = 1.0f;
    result(2) = 0.0f;
  }
  else {
    auto denominator = ::sqrt(aVector(1) * aVector(1) + aVector(2) * aVector(2));
    result(1) = -aVector(2) / denominator;
    result(2) =  aVector(1) / denominator;
  }
  return result;
}

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

template<typename tReal>
Vector<tReal> getAltitude(Vertex<tReal> const &aCommon1, Vertex<tReal> const &aCommon2, Vertex<tReal> const &aIndependent) {
  auto commonVector = aCommon2 - aCommon1;
  auto independentVector = aIndependent - aCommon1;
  auto footFactor = commonVector.dot(independentVector) / commonVector.squaredNorm();
  return independentVector - commonVector * footFactor;
}

template<typename tReal>
uint32_t toWhichSide(Vertex<tReal> const &aStart, Vertex<tReal> const &aEnd) {
  uint32_t result = 3u;
  tReal denom = aStart(0) - aEnd(0) + aStart(1) - aEnd(1);
  if(::abs(denom) > cgGeneralEpsilon<tReal>) {
    auto ratio = ((aStart(0) - 1.0f) * aEnd(1) - aStart(1) * (aEnd(0) - 1.0f)) / denom;
    auto direction = (aStart(0) + aStart(1) - 1.0f) / denom;
    result = (ratio > -cgGeneralEpsilon<tReal> && ratio < 1.0f + cgGeneralEpsilon<tReal> && direction > 0.0f) ? 0u : result;
  }
  else { // nothing to do
  }
  denom = aStart(1) - aEnd(1) + aStart(2) - aEnd(2);
  if(::abs(denom) > cgGeneralEpsilon<tReal>) {
    auto ratio = ((aStart(1) - 1.0f) * aEnd(2) - aStart(2) * (aEnd(1) - 1.0f)) / denom;
    auto direction = (aStart(1) + aStart(2) - 1.0f) / denom;
    result = (ratio > -cgGeneralEpsilon<tReal> && ratio < 1.0f + cgGeneralEpsilon<tReal> && direction > 0.0f) ? 1u : result;
  }
  else { // nothing to do
  }
  denom = aStart(2) - aEnd(2) + aStart(0) - aEnd(0);
  if(::abs(denom) > cgGeneralEpsilon<tReal>) {
    auto ratio = ((aStart(2) - 1.0f) * aEnd(0) - aStart(0) * (aEnd(2) - 1.0f)) / denom;
    auto direction = (aStart(2) + aStart(0) - 1.0f) / denom;
    result = (ratio > -cgGeneralEpsilon<tReal> && ratio < 1.0f + cgGeneralEpsilon<tReal> && direction > 0.0f) ? 2u : result;
  }
  else { // nothing to do
  }
  return result;
}

template<typename tReal>
tReal Ray<tReal>::getAverageErrorSquared(std::vector<Vertex<tReal>> const &aPoints) const {
  tReal sum = 0.0f;
  for(auto const &point : aPoints) {
    sum += (point - mStart - (point - mStart).dot(mDirection) * mDirection).squaredNorm();
  }
  return aPoints.size() == 0u ? 0u : sum / aPoints.size();
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

template<typename tReal>
Intersection<tReal> Plane<tReal>::intersect(Ray<tReal> const &aRay) const {
  Intersection<tReal> result;
  result.mCosIncidence = aRay.mDirection.dot(mNormal);
  if(::abs(result.mCosIncidence) >= csRayPlaneIntersectionEpsilon) {
    result.mDistance = (mConstant - mNormal.dot(aRay.mStart)) / result.mCosIncidence;
    result.mValid = true;
    result.mPoint = aRay.mStart + result.mDistance * aRay.mDirection;
  }
  else {
    result.mValid = false;
  }
  return result;
}

template<typename tReal>
void Plane<tReal>::makeDistancePositive(Vector<tReal> const aPoint) {
  if(distance(aPoint) < 0.0f) {
    mNormal = -mNormal;
    mConstant = -mConstant;
  }
  else { // Nothing to do
  }
}

template<typename tReal>
void Plane<tReal>::makeDistanceNegative(Vector<tReal> const aPoint) {
  if(distance(aPoint) > 0.0f) {
    mNormal = -mNormal;
    mConstant = -mConstant;
  }
  else { // Nothing to do
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
