#include "util.h"


Vector getNormal(Triangle const &aFace) {
  return (aFace[1] - aFace[0]).cross(aFace[2] - aFace[0]);
}

Vector getNormal(Vertex const &aVertex0, Vertex const &aVertex1, Vertex const &aVertex2) {
  return (aVertex1 - aVertex0).cross(aVertex2 - aVertex0);
}

Vector getAperpendicular(Vector const &aVector) {
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

float getPerimeter(Triangle const &aTriangle) {
  return (aTriangle[0u] - aTriangle[1u]).norm() + (aTriangle[1u] - aTriangle[2u]).norm() + (aTriangle[2u] - aTriangle[0u]).norm();
}

// Assumes sum(aB*) == 1
Vertex barycentric2cartesian(Triangle const &aTriangle, float const aB0, float const aB1, float const aB2) {
  return aTriangle[0u] * aB0 + aTriangle[1u] * aB1 + aTriangle[2u] * aB2;
}

Vertex barycentric2cartesian(Triangle const &aTriangle, float const aB0, float const aB1) {
  return aTriangle[0u] * aB0 + aTriangle[1u] * aB1 + aTriangle[2u] * (1.0f - aB0 - aB1);
}

// Assumes sum(aB*) == 1
Vertex barycentric2cartesian(Vertex const &aV0, Vertex const &aV1, Vertex const &aV2, float const aB0, float const aB1, float const aB2) {
  return aV0 * aB0 + aV1 * aB1 + aV2 * aB2;
}

Vertex barycentric2cartesian(Vertex const &aV0, Vertex const &aV1, Vertex const &aV2, float const aB0, float const aB1) {
  return aV0 * aB0 + aV1 * aB1 + aV2 * (1.0f - aB0 - aB1);
}

Matrix getBarycentricInverse(Vertex const &aVertex0, Vertex const &aVertex1, Vertex const &aVertex2) {
  Matrix vertices {
    { aVertex0(0), aVertex1(0), aVertex2(0) },
    { aVertex0(1), aVertex1(1), aVertex2(1) },
    { aVertex0(2), aVertex1(2), aVertex2(2) },
  };
  return vertices.inverse();
}

Vector getAltitude(Vertex const &aCommon1, Vertex const &aCommon2, Vertex const &aIndependent) {
  auto commonVector = aCommon2 - aCommon1;
  auto independentVector = aIndependent - aCommon1;
  auto footFactor = commonVector.dot(independentVector) / commonVector.squaredNorm();
  return independentVector - commonVector * footFactor;
}

uint32_t toWhichSide(Vertex const &aStart, Vertex const &aEnd) {
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

float Ray::getAverageErrorSquared(std::vector<Vertex> const &aPoints) const {
  float sum = 0.0f;
  for(auto const &point : aPoints) {
    sum += getDistance2(point);
  }
  return aPoints.size() == 0u ? 0u : sum / aPoints.size();
}

// If aProportion < 0.5, the result will be closer to aPoint0
Plane Plane::createFrom1proportion2points(float const aProportion, Vertex const &aPoint0, Vertex const &aPoint1) {
  Plane result;
  result.mNormal = (aPoint1 - aPoint0).normalized();
  result.mConstant = result.mNormal.dot(aPoint1 * aProportion + aPoint0 * (1.0f - aProportion));
  return result;
}

Plane Plane::createFrom3points(Vertex const &aPoint0, Vertex const &aPoint1, Vertex const &aPoint2) {
  Plane result;
  result.mNormal = (aPoint1 - aPoint0).cross(aPoint2 - aPoint0).normalized();
  result.mConstant = result.mNormal.dot(aPoint0);
  return result;
}

Plane Plane::createFrom1vector2points(Vector const &aDirection, Vertex const &aPoint0, Vertex const &aPoint1) {
  Plane result;
  result.mNormal = aDirection.cross(aPoint1 - aPoint0).normalized();
  result.mConstant = result.mNormal.dot(aPoint0);
  return result;
}

Plane Plane::createFrom2vectors1point(Vertex const &aDirection0, Vertex const &aDirection1, Vertex const &aPoint) {
  Plane result;
  result.mNormal = aDirection0.cross(aDirection1).normalized();
  result.mConstant = result.mNormal.dot(aPoint);
  return result;
}

Vertex Plane::intersect(Plane const &aPlane0, Plane const &aPlane1, Plane const &aPlane2) {
  Matrix matrix {
    { aPlane0.mNormal(0), aPlane0.mNormal(1), aPlane0.mNormal(2) },
    { aPlane1.mNormal(0), aPlane1.mNormal(1), aPlane1.mNormal(2) },
    { aPlane2.mNormal(0), aPlane2.mNormal(1), aPlane2.mNormal(2) }
  };
  Vector vector{ aPlane0.mConstant, aPlane1.mConstant, aPlane2.mConstant };
  Vertex result = matrix.inverse() * vector;
  // TODO remove
  auto error =::abs((matrix * result - vector).norm() / vector.norm());
  static_assert(std::is_same_v<decltype(error), float>);
  if(error > 0.1f) {
    throw "Plane::intersect rel error";
  }
  return result;
}

Intersection  Plane::intersect(Vertex const &aStart, Vector const aDirection) const {
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

void Plane::makeDistancePositive(Vector const aPoint) {
  if(distance(aPoint) < 0.0f) {
    mNormal = -mNormal;
    mConstant = -mConstant;
  }
  else { // Nothing to do
  }
}

void Plane::makeDistanceNegative(Vector const aPoint) {
  if(distance(aPoint) > 0.0f) {
    mNormal = -mNormal;
    mConstant = -mConstant;
  }
  else { // Nothing to do
  }
}
