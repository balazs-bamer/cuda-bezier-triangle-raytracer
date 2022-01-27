#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_UTIL
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_UTIL

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t

#include <Eigen/Dense>

#ifndef UTIL_CUDA_PREFIX_DEVICE  // When defined, should be __device__
#define UTIL_CUDA_PREFIX_DEVICE
#endif

#ifndef UTIL_CUDA_PREFIX_HOST  // When defined, should be __host__
#define UTIL_CUDA_PREFIX_HOST
#endif


constexpr float cgPi = 3.14159265358979323846;
constexpr float cgGeneralEpsilon = 1.0e-5;


using Vector         = Eigen::Matrix<float, 3, 1>;
using Vertex         = Eigen::Matrix<float, 3, 1>;
using Matrix         = Eigen::Matrix<float, 3, 3>;
using Transform      = Eigen::Matrix<float, 3, 3>;
using Triangle       = std::array<Vertex, 3u>;


// Namespace-like struct to allow this file act as a header-only library for CPP and CUDA.
struct util {
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Vector getNormal(Triangle const &aFace) {
    return (aFace[1] - aFace[0]).cross(aFace[2] - aFace[0]);
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Vector getNormal(Vertex const &aVertex0, Vertex const &aVertex1, Vertex const &aVertex2) {
    return (aVertex1 - aVertex0).cross(aVertex2 - aVertex0);
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static float getPerimeter(Triangle const &aTriangle) {
    return (aTriangle[0u] - aTriangle[1u]).norm() + (aTriangle[1u] - aTriangle[2u]).norm() + (aTriangle[2u] - aTriangle[0u]).norm();
  }

  // Assumes sum(aB*) == 1
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Vertex barycentric2cartesian(Triangle const &aTriangle, float const aB0, float const aB1, float const aB2) {
    return aTriangle[0u] * aB0 + aTriangle[1u] * aB1 + aTriangle[2u] * aB2;
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Vertex barycentric2cartesian(Triangle const &aTriangle, float const aB0, float const aB1) {
    return aTriangle[0u] * aB0 + aTriangle[1u] * aB1 + aTriangle[2u] * (1.0f - aB0 - aB1);
  }

  // Assumes sum(aB*) == 1
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Vertex barycentric2cartesian(Vertex const &aV0, Vertex const &aV1, Vertex const &aV2, float const aB0, float const aB1, float const aB2) {
    return aV0 * aB0 + aV1 * aB1 + aV2 * aB2;
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Vertex barycentric2cartesian(Vertex const &aV0, Vertex const &aV1, Vertex const &aV2, float const aB0, float const aB1) {
    return aV0 * aB0 + aV1 * aB1 + aV2 * (1.0f - aB0 - aB1);
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Matrix getBarycentricInverse(Vertex const &aVertex0, Vertex const &aVertex1, Vertex const &aVertex2) {
    Matrix vertices {
      { aVertex0(0), aVertex1(0), aVertex2(0) },
      { aVertex0(1), aVertex1(1), aVertex2(1) },
      { aVertex0(2), aVertex1(2), aVertex2(2) },
    };
    return vertices.inverse();
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Vector getAperpendicular(Vector const &aVector) {
    static constexpr float csEpsilon = 1e-10;

    Vector result;
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

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  template<typename tLambda>
  static void divide(Triangle const &aTriangle, int32_t const aDivisor, tLambda &&aCollector) {
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

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Vector getAltitude(Vertex const &aCommon1, Vertex const &aCommon2, Vertex const &aIndependent) {
    auto commonVector = aCommon2 - aCommon1;
    auto independentVector = aIndependent - aCommon1;
    auto footFactor = commonVector.dot(independentVector) / commonVector.squaredNorm();
    return independentVector - commonVector * footFactor;
  }

  // Takes barycentric ends of a vector with aStart inside the triangle, finds out which side it will intersect.
  // 0 between 300 and 030
  // 1 between 030 and 003
  // 2 between 003 and 300
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static uint32_t toWhichSide(Vertex const &aStart, Vertex const &aEnd) {
    uint32_t result = 3u;
    float denom = aStart(0) - aEnd(0) + aStart(1) - aEnd(1);
    if(::abs(denom) > cgGeneralEpsilon) {
      auto ratio = ((aStart(0) - 1.0f) * aEnd(1) - aStart(1) * (aEnd(0) - 1.0f)) / denom;
      auto direction = (aStart(0) + aStart(1) - 1.0f) / denom;
      result = (ratio > -cgGeneralEpsilon && ratio < 1.0f + cgGeneralEpsilon && direction > 0.0f) ? 0u : result;
    }
    else { // nothing to do
    }
    denom = aStart(1) - aEnd(1) + aStart(2) - aEnd(2);
    if(::abs(denom) > cgGeneralEpsilon) {
      auto ratio = ((aStart(1) - 1.0f) * aEnd(2) - aStart(2) * (aEnd(1) - 1.0f)) / denom;
      auto direction = (aStart(1) + aStart(2) - 1.0f) / denom;
      result = (ratio > -cgGeneralEpsilon && ratio < 1.0f + cgGeneralEpsilon && direction > 0.0f) ? 1u : result;
    }
    else { // nothing to do
    }
    denom = aStart(2) - aEnd(2) + aStart(0) - aEnd(0);
    if(::abs(denom) > cgGeneralEpsilon) {
      auto ratio = ((aStart(2) - 1.0f) * aEnd(0) - aStart(0) * (aEnd(2) - 1.0f)) / denom;
      auto direction = (aStart(2) + aStart(0) - 1.0f) / denom;
      result = (ratio > -cgGeneralEpsilon && ratio < 1.0f + cgGeneralEpsilon && direction > 0.0f) ? 2u : result;
    }
    else { // nothing to do
    }
    return result;
  }
};


struct Ray final {
  Vertex mStart;
  Vector mDirection; // normalized

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Ray() = default;

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Ray(Vertex const &aStart, Vector const &aDirection)
  : mStart(aStart)
  , mDirection(aDirection.normalized()) {}

  // Does not consider that a ray is a half-line.
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Vector getPerpendicularTo(Vertex const &aPoint) const {
    return aPoint - mStart - (aPoint - mStart).dot(mDirection) * mDirection;
  }

  // Does not consider that a ray is a half-line.
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  float getDistance(Vertex const &aPoint) const {
    return getPerpendicularTo(aPoint).norm();
  }

  // Does not consider that a ray is a half-line.
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  float getDistance2(Vertex const &aPoint) const {
    return getPerpendicularTo(aPoint).squaredNorm();
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  float getAverageErrorSquared(std::vector<Vertex> const &aPoints) const {
    float sum = 0.0f;
    for(auto const &point : aPoints) {
      sum += getDistance2(point);
    }
    return aPoints.size() == 0u ? 0u : sum / aPoints.size();
  }
};


struct Intersection final {
  bool   mValid;         // Will be true even if distance < 0, because it can be important.
  Vertex mPoint;
  float  mCosIncidence;  // Negative if the ray comes from outside, so the surface normal and the ray point in opposite direction.
  float  mDistance;      // Distance of ray source point and intersection point.
};


// Plane equation is in the form of point.dot(mNormal) == mConstant where point is any point in the plane.
struct Plane final {
  static constexpr float csRayPlaneIntersectionEpsilon = 0.00001f;

  Vector mNormal;        // normalized
  float  mConstant;

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Plane() = default;

  // aNormal must be normalized.
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Plane(Vector const &aNormal, float const aConstant) : mNormal(aNormal), mConstant(aConstant) {}

  // If aProportion < 0.5, the result will be closer to aPoint0
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Plane         createFrom1proportion2points(float const aProportion, Vertex const &aPoint0, Vertex const &aPoint1) {
    Plane result;
    result.mNormal = (aPoint1 - aPoint0).normalized();
    result.mConstant = result.mNormal.dot(aPoint1 * aProportion + aPoint0 * (1.0f - aProportion));
    return result;
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Plane         createFrom3points(Vertex const &aPoint0, Vertex const &aPoint1, Vertex const &aPoint2) {
    Plane result;
    result.mNormal = (aPoint1 - aPoint0).cross(aPoint2 - aPoint0).normalized();
    result.mConstant = result.mNormal.dot(aPoint0);
    return result;
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Plane         createFromTriangle(Triangle const &aTriangle) { return createFrom3points(aTriangle[0u], aTriangle[1u], aTriangle[2u]); }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Plane         createFrom1vector2points(Vector const &aDirection, Vertex const &aPoint0, Vertex const &aPoint1) {
    Plane result;
    result.mNormal = aDirection.cross(aPoint1 - aPoint0).normalized();
    result.mConstant = result.mNormal.dot(aPoint0);
    return result;
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Plane         createFrom2vectors1point(Vertex const &aDirection0, Vertex const &aDirection1, Vertex const &aPoint) {
    Plane result;
    result.mNormal = aDirection0.cross(aDirection1).normalized();
    result.mConstant = result.mNormal.dot(aPoint);
    return result;
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  static Vertex intersect(Plane const &aPlane0, Plane const &aPlane1, Plane const &aPlane2) {
    Matrix matrix {
      { aPlane0.mNormal(0), aPlane0.mNormal(1), aPlane0.mNormal(2) },
      { aPlane1.mNormal(0), aPlane1.mNormal(1), aPlane1.mNormal(2) },
      { aPlane2.mNormal(0), aPlane2.mNormal(1), aPlane2.mNormal(2) }
    };
    Vector vector{ aPlane0.mConstant, aPlane1.mConstant, aPlane2.mConstant };
    return matrix.inverse() * vector;
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Intersection  intersect(Vertex const &aStart, Vector const aDirection) const {
    Intersection result;
    result.mCosIncidence = aDirection.dot(mNormal);
    if(::abs(result.mCosIncidence) >= csRayPlaneIntersectionEpsilon) {
      result.mDistance = (mConstant - mNormal.dot(aStart)) / result.mCosIncidence;
      result.mValid = true;
      result.mPoint = aStart + result.mDistance * aDirection;
    }
    else {
      result.mValid = false;
    }
    return result;
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Intersection  intersect(Ray const &aRay) const { return intersect(aRay.mStart, aRay.mDirection); }

  // Point projection on this plane
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Vector        project(Vector const &aPoint) const { return aPoint - mNormal * (aPoint.dot(mNormal) - mConstant); }

  // Distance of point and this plane, >0 if the normal points towards the point.
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  float                distance(Vector const &aPoint) const { return aPoint.dot(mNormal) - mConstant; }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  void                 makeDistancePositive(Vector const aPoint) {
    if(distance(aPoint) < 0.0f) {
      mNormal = -mNormal;
      mConstant = -mConstant;
    }
    else { // Nothing to do
    }
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  void                 makeDistanceNegative(Vector const aPoint) {
    if(distance(aPoint) > 0.0f) {
      mNormal = -mNormal;
      mConstant = -mConstant;
    }
    else { // Nothing to do
    }
  }

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
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
  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Spherical(float const aX, float const aY, float const aZ)
  : mR(::sqrt(aX * aX + aY * aY + aZ * aZ))
  , mInclination(::acos(aZ / mR))
  , mAzimuth(::atan2(aY, aX)) {}
};


struct Sphere final {  // TODO consider if we need this for bounding sphere intersection or do it by hand. Problem is ray is a half-line and this would also need to store r2.
  Vector mCenter;
  float  mRadius;

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Sphere() = default;

  UTIL_CUDA_PREFIX_DEVICE UTIL_CUDA_PREFIX_HOST
  Sphere(Vector const &aCenter, float const aRadius) : mCenter(aCenter), mRadius(aRadius) {}

  bool doesIntersect(Ray const aRay);
};


template<uint32_t tDegree>
class PolynomApprox final {
  static_assert(tDegree >= 1u);

private:
  static constexpr uint32_t csDegree1 = tDegree + 1u;

  Eigen::Matrix<float, tDegree + 1u, 1u> mCoefficients;
  float                                  mRrmsError;    // https://stats.stackexchange.com/questions/413209/is-there-something-like-a-root-mean-square-relative-error-rmsre-or-what-is-t

public:
  PolynomApprox(std::vector<float> const& aSamplesX, std::vector<float> const aSamplesY) {
    std::array<float, tDegree + csDegree1>     fitX;
    Eigen::Matrix<float, csDegree1, csDegree1> fitA;
    Eigen::Matrix<float, csDegree1 , 1u>       fitB;

    auto sampleCount = std::min(aSamplesX.size(), aSamplesY.size());
    auto degree1 = tDegree + 1u;

    if(sampleCount >= csDegree1) {
      for(uint32_t i = 0u; i < csDegree1 + tDegree; ++i) {
        fitX[i] = 0.0f;
        for(uint32_t j = 0u; j < sampleCount; ++j) {
          fitX[i] += ::pow(aSamplesX[j], i);
        }
      }
      for(uint32_t i = 0u; i < csDegree1; ++i) {
        for(uint32_t j = 0u; j < csDegree1; ++j) {
          fitA(i, j) =fitX[i + j];
        }
      }
      for(uint32_t i = 0u; i < csDegree1; ++i) {
        float tmp = 0.0f;
        for(uint32_t j = 0u; j < sampleCount; ++j) {
          tmp += ::pow(aSamplesX[j], i) * aSamplesY[j];
        }
        fitB(i) = tmp;
      }
      mCoefficients = fitA.colPivHouseholderQr().solve(fitB);
      auto diffs = 0.0f;
      auto desireds = 0.0f;
      for(uint32_t i = 0u; i < sampleCount; ++i) {
        auto desired = aSamplesY[i];
        auto diff = (desired - eval(aSamplesX[i]));
        diffs += diff * diff;
        desireds += desired * desired;
      }
      mRrmsError = (desireds > 0.0f ? ::sqrt(diffs / desireds / sampleCount) : 0.0f);
      }
    else {
      mCoefficients = Eigen::Matrix<float, tDegree + 1u, 1u>::Zero();
      mRrmsError = std::numeric_limits<float>::max();
    }
  }

  float getRrmsError() const { return mRrmsError; }

  float eval(float const aX) const {
    float result = mCoefficients(tDegree) * aX;
    for(uint32_t i = tDegree - 1u; i > 0u; --i) {
      result = (result + mCoefficients(i)) * aX;
    }
    result += mCoefficients(0);
    return result;
  }
};

#endif
