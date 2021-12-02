#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERLENS
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERLENS

#include "bezierMesh.h"

#include <iostream> // TODO REMOVE
#include <iomanip> // TODO REMOVE

enum class RefractionResult : uint8_t {
  cNone    = 0u,
  cInside  = 1u,
  cOutside = 2u
};

// This class first will be used only to simulate discrete lenses. Compounds will come later, if any.
template<typename tReal>
class BezierLens final {
private:
  static constexpr tReal csMinCos2refraction = 0.01f; // Squared, approximately 84 degrees.
  static constexpr tReal csQuadCoeffEpsilon  = 1e-7f;

  tReal             mRefractiveIndex;
  BezierMesh<tReal> mMesh;

  using Vector              = ::Vector<tReal>;
  using Vertex              = ::Vertex<tReal>;
  using Plane               = ::Plane<tReal>;
  using Ray                 = ::Ray<tReal>;
  using BezierIntersection  = ::BezierIntersection<tReal>;

public:
  BezierLens(tReal const aRi, BezierMesh<tReal> const & aMesh) : mRefractiveIndex(aRi), mMesh(aMesh) {}
  BezierLens(tReal const aRi, BezierMesh<tReal> && aMesh)      : mRefractiveIndex(aRi), mMesh(std::move(aMesh)) {}

  // bool indicates wether the ray is inside the lens AFTER refraction
  std::pair<Ray, RefractionResult> refract(Ray const &aRay) const;
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

// The ray either is refracted, or gets lost (if misses the lens or suffers total reflection).
template<typename tReal>
std::pair<Ray<tReal>, RefractionResult> BezierLens<tReal>::refract(Ray const &aRay) const {
  Ray result;
  RefractionResult statusLater;
  auto const intersect = mMesh.intersect(aRay);
  if(intersect.mWhat == BezierIntersection::What::cIntersect) {
std::cout << intersect.mIntersection.mCosIncidence << '\n';
    statusLater = (intersect.mIntersection.mCosIncidence < 0.0f ? RefractionResult::cInside : RefractionResult::cOutside);
    result.mStart = intersect.mIntersection.mPoint;
    auto const directionFactor = (statusLater == RefractionResult::cInside ? -1.0f : 1.0f);
    auto const effectiveRefractiveFactor = (statusLater == RefractionResult::cInside ? 1.0f / mRefractiveIndex : mRefractiveIndex);
    auto const cos2refraction = 1.0f + effectiveRefractiveFactor * effectiveRefractiveFactor * (intersect.mIntersection.mCosIncidence * intersect.mIntersection.mCosIncidence - 1.0f);
    if(cos2refraction > csMinCos2refraction) {
      auto const cosIncidence = ::abs(intersect.mIntersection.mCosIncidence);
      auto const cosRefraction = ::sqrt(cos2refraction);
      auto const &incidence = aRay.mDirection;
      auto const normal = intersect.mNormal * directionFactor;
      //  a * incidence + b * normal      will be the result direction, we look for a and b
      //  cosRefraction  =  normal . (a * incidence + b * normal)
      //  cosRefraction  =  a * cosIncidence + b
      //  | a * incidence + b * normal |  =  1
      //  | a * (incidence - cosIncidence * normal) + cosRefraction * normal |  =  1
      //  so this gives a quadratic equation for a when we take the x, y, z components.
      auto quadCoeffA = 0.0f;
      auto quadCoeffB = 0.0f;
      auto quadCoeffC = 0.0f;
      for(uint32_t i = 0u; i < 3u; ++i) {
        auto const tmp = incidence(i) - cosIncidence * normal(i);
        quadCoeffA += tmp * tmp;
        quadCoeffB += tmp * normal(i);
        quadCoeffC += normal(i) * normal(i);
      }
      quadCoeffB *= 2.0f * cosRefraction;
      quadCoeffC *= cosRefraction * cosRefraction;
      quadCoeffC -= 1.0f;
      tReal a;
      tReal b;
      if(::abs(quadCoeffA) < csQuadCoeffEpsilon) {
        a = -quadCoeffC / quadCoeffB;                // both quadCoeffA and quadCoeffB can't be 0
        b = cosRefraction - a * cosIncidence;
      }
      else {
        auto const discriminantSqrt = ::sqrt(quadCoeffB * quadCoeffB - 4.0f * quadCoeffA * quadCoeffC); // Will definitely have a solution.
        auto const quadCoeffAdouble = 2.0f * quadCoeffA;
        a = (-quadCoeffB - discriminantSqrt) / quadCoeffAdouble;
        b = cosRefraction - a * cosIncidence;
        if(a < 0.0f || b < 0.0f) {
          a = (-quadCoeffB + discriminantSqrt) / quadCoeffAdouble;
          b = cosRefraction - a * cosIncidence;
        }
        else { // Nothing to do
        }
      }
if(a < 0.0f || b < 0.0f) { throw "TODO REMOVE invalid refraction direction coefficients"; }
      result.mDirection = a * incidence + b * normal;
    }
    else {
      statusLater = RefractionResult::cNone;  // Total reflection, or close to it. Ignored.
    }
  }
  else {
    statusLater = RefractionResult::cNone;
  }
  return std::make_pair(result, statusLater);
}


#endif // BEZIERLENS_H
