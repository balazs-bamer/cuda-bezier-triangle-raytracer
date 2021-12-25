#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERLENS
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERLENS

#include "bezierMesh.h"


enum class RefractionResult : uint32_t {
  cNone    = 0u,
  cInside  = 1u,
  cOutside = 2u
};

// This class first will be used only to simulate discrete lenses. Compounds will come later, if any.
class BezierLens final {
private:
  static constexpr float csMaxSin2refraction = 0.99f;  // Squared, approximately 82 degrees.
  static constexpr float csMinSin2refraction = 1e-12f; // Don't calculate refraction below this, let the ray pass as it came.

  float      mRefractiveIndex;
  BezierMesh mMesh;

public:
  BezierLens(float const aRi, BezierMesh const & aMesh) : mRefractiveIndex(aRi), mMesh(aMesh) {}
  BezierLens(float const aRi, BezierMesh&& aMesh)      : mRefractiveIndex(aRi), mMesh(std::move(aMesh)) {}

  // The ray either is refracted, or gets lost (if misses the lens or suffers total reflection).
  std::pair<Ray, RefractionResult> refract(Ray const &aRay, RefractionResult const aExpected) const;
};

#endif // BEZIERLENS_H
