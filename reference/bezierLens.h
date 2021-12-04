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
  static constexpr tReal csMaxSin2refraction = 0.99f;  // Squared, approximately 82 degrees.
  static constexpr tReal csMinSin2refraction = 1e-12f; // Don't calculate refraction below this, let the ray pass as it came.

  tReal             mRefractiveIndex;
  BezierMesh<tReal> mMesh;

  using Vector              = ::Vector<tReal>;
  using Vertex              = ::Vertex<tReal>;
  using Matrix              = ::Matrix<tReal>;
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
std::cout << "ray: " << aRay.mDirection(0) << ' ' << aRay.mDirection(1) << ' ' << aRay.mDirection(2) << '\n';
std::cout << "cos: " << intersect.mIntersection.mCosIncidence << " acos: " << 180.0f/3.1415f * ::acos(intersect.mIntersection.mCosIncidence);
    statusLater = (intersect.mIntersection.mCosIncidence < 0.0f ? RefractionResult::cInside : RefractionResult::cOutside);
    result.mStart = intersect.mIntersection.mPoint;
    auto const effectiveRefractiveFactor = (statusLater == RefractionResult::cInside ? 1.0f / mRefractiveIndex : mRefractiveIndex);
    auto const sin2refraction = effectiveRefractiveFactor * effectiveRefractiveFactor * (1.0f - intersect.mIntersection.mCosIncidence * intersect.mIntersection.mCosIncidence);
    if(sin2refraction < csMaxSin2refraction) {
      if(sin2refraction > csMinSin2refraction) {
        auto const directionFactor = (statusLater == RefractionResult::cInside ? 1.0f : -1.0f);
std::cout << " df: " << directionFactor << '\n';
        auto const normal = intersect.mNormal * directionFactor;
std::cout << "normal: " << normal(0) << ' ' << normal(1) << ' ' << normal(2) << '\n';
        auto cos1 = ::abs(intersect.mIntersection.mCosIncidence);
        auto cos2 = ::sqrt(1.0f - sin2refraction);
        result.mDirection = (aRay.mDirection * effectiveRefractiveFactor + normal * (effectiveRefractiveFactor * cos1 - cos2)).normalized();
std::cout << "direct: " << result.mDirection(0) << ' ' << result.mDirection(1) << ' ' << result.mDirection(2) << '\n';
std::cout << "teta1: " << 180.0f/3.1415f * ::asin(::sqrt(1.0f - cos1 * cos1)) <<
            " teta2: " << 180.0f/3.1415f * ::acos(cos2) <<
            " teta2: " << 180.0f/3.1415f * ::acos(result.mDirection.dot(-normal)) << '\n';
      }
      else {
        result.mDirection = aRay.mDirection;        // Too little incidence, continue in the same direction.
      }
    }
    else {
      statusLater = RefractionResult::cNone;        // Total reflection, or close to it. Ignored.
    }
  }
  else {
    statusLater = RefractionResult::cNone;
  }
std::cout << " sl: " << static_cast<int>(statusLater) << '\n';
  return std::make_pair(result, statusLater);
}


#endif // BEZIERLENS_H
