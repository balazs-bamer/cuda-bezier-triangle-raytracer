#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_UTIL
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_UTIL

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t

#include <Eigen/Dense>


constexpr float cgPi = 3.14159265358979323846;
constexpr float cgGeneralEpsilon = 1.0e-5;

using Vector         = Eigen::Matrix<float, 3, 1>;
using Vertex         = Eigen::Matrix<float, 3, 1>;
using Matrix         = Eigen::Matrix<float, 3, 3>;
using Transform      = Eigen::Matrix<float, 3, 3>;
using Triangle       = std::array<Vertex, 3u>;

Vector getNormal(Triangle const &aFace);
Vector getNormal(Vertex const &aVertex0, Vertex const &aVertex1, Vertex const &aVertex2);

Vector getAperpendicular(Vector const &aVector);

float getPerimeter(Triangle const &aTriangle);

// Assumes sum(aB*) == 1
Vertex barycentric2cartesian(Triangle const &aTriangle, float const aB0, float const aB1, float const aB2);
Vertex barycentric2cartesian(Triangle const &aTriangle, float const aB0, float const aB1);

// Assumes sum(aB*) == 1
Vertex barycentric2cartesian(Vertex const &aV0, Vertex const &aV1, Vertex const &aV2, float const aB0, float const aB1, float const aB2);
Vertex barycentric2cartesian(Vertex const &aV0, Vertex const &aV1, Vertex const &aV2, float const aB0, float const aB1);

Matrix getBarycentricInverse(Vertex const &aVertex0, Vertex const &aVertex1, Vertex const &aVertex2);

Vector getAltitude(Vertex const &aCommon1, Vertex const &aCommon2, Vertex const &aIndependent);

template<typename tLambda>
void divide(Triangle const &aTriangle, int32_t const aDivisor, tLambda &&aCollector);

// Takes barycentric ends of a vector with aStart inside the triangle, finds out which side it will intersect.
// 0 between 300 and 030
// 1 between 030 and 003
// 2 between 003 and 300
uint32_t toWhichSide(Vertex const &aStart, Vertex const &aEnd);

struct Ray final {
  Vertex mStart;
  Vector mDirection; // normalized

  Ray() = default;
  Ray(Vertex const &aStart, Vector const &aDirection)
  : mStart(aStart)
  , mDirection(aDirection.normalized()) {}

  Vector getPerpendicularTo(Vertex const &aPoint) const {
    return aPoint - mStart - (aPoint - mStart).dot(mDirection) * mDirection;
  }

  float getDistance(Vertex const &aPoint) const {
    return getPerpendicularTo(aPoint).norm();
  }

  float getDistance2(Vertex const &aPoint) const {
    return getPerpendicularTo(aPoint).squaredNorm();
  }

  float getAverageErrorSquared(std::vector<Vertex> const &aPoints) const;
};

struct Intersection final {
  bool          mValid;         // Will be true even if distance < 0, because it can be important.
  Vertex mPoint;
  float         mCosIncidence;  // Negative if the ray comes from outside, so the surface normal and the ray point in opposite direction.
  float         mDistance;      // Distance of ray source point and intersection point.
};

// Plane equation is in the form of point.dot(mNormal) == mConstant where point is any point in the plane.
struct Plane final {
  static constexpr float csRayPlaneIntersectionEpsilon = 0.00001f;

  Vector mNormal;        // normalized
  float         mConstant;

  Plane() = default;
  // aNormal must be normalized.
  Plane(Vector const &aNormal, float const aConstant) : mNormal(aNormal), mConstant(aConstant) {}

  static Plane         createFrom1proportion2points(float const aProportion, Vertex const &aPoint0, Vertex const &aPoint1);
  static Plane         createFrom3points(Vertex const &aPoint0, Vertex const &aPoint1, Vertex const &aPoint2);
  static Plane         createFromTriangle(Triangle const &aTriangle) { return createFrom3points(aTriangle[0u], aTriangle[1u], aTriangle[2u]); }
  static Plane         createFrom1vector2points(Vector const &aDirection, Vertex const &aPoint0, Vertex const &aPoint1);
  static Plane         createFrom2vectors1point(Vertex const &aDirection0, Vertex const &aDirection1, Vertex const &aPoint);
  static Vertex intersect(Plane const &aPlane0, Plane const &aPlane1, Plane const &aPlane2);
  Intersection  intersect(Vertex const &aStart, Vector const aDirection) const;
  Intersection  intersect(Ray const &aRay) const { return intersect(aRay.mStart, aRay.mDirection); }

  // Point projection on this plane
  Vector        project(Vector const &aPoint) const { return aPoint - mNormal * (aPoint.dot(mNormal) - mConstant); }

  // Distance of point and this plane, >0 if the normal points towards the point.
  float                distance(Vector const &aPoint) const { return aPoint.dot(mNormal) - mConstant; }

  void                 makeDistancePositive(Vector const aPoint);
  void                 makeDistanceNegative(Vector const aPoint);

  bool operator<(Plane const &aOther) const { return mNormal(0) < aOther.mNormal(0) // TODO remove: only for debugging
                                                         || mNormal(0) == aOther.mNormal(0) && mNormal(1) < aOther.mNormal(1)
                                                         || mNormal(0) == aOther.mNormal(0) && mNormal(1) == aOther.mNormal(1) && mNormal(2) < aOther.mNormal(2)
                                                         || mNormal(0) == aOther.mNormal(0) && mNormal(1) == aOther.mNormal(1) && mNormal(2) == aOther.mNormal(2) && mConstant < aOther.mConstant; }
};

struct Spherical final {
  float mR;
  float mAzimuth; // or sector
  float mInclination;  // or belt

  // Assumes the conversion can be done.
  Spherical(float const aX, float const aY, float const aZ)
  : mR(::sqrt(aX * aX + aY * aY + aZ * aZ))
  , mInclination(::acos(aZ / mR))
  , mAzimuth(::atan2(aY, aX)) {}
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tLambda>
void divide(Triangle const &aTriangle, int32_t const aDivisor, tLambda &&aCollector) {
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

#endif
