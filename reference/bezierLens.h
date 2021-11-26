#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERLENS
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERLENS

#include "bezierMesh.h"

#include <iostream> // TODO REMOVE
#include <iomanip> // TODO REMOVE

// This class first will be used only to simulate discrete lenses. Compounds will come later, if any.
template<typename tReal>
class BezierLens final {
private:
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
  std::pair<Ray, bool> refract(Ray const &aRay) const;
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tReal>
std::pair<Ray<tReal>, bool> BezierLens<tReal>::refract(Ray const &aRay) const {
  Ray result;
  bool inside;
  auto intersect = mMesh.intersect(aRay);
  if(intersect.mWhat == BezierIntersection::What::cIntersect) {
std::cout << intersect.mIntersection.mCosIncidence << '\n';
    inside = intersect.mIntersection.mCosIncidence < 0.0f;
    result.mStart = intersect.mIntersection.mPoint;  // TODO implement
    result.mDirection = aRay.mDirection;
  }
  else {
    result = aRay;
    inside = false;
  }
  return std::make_pair(result, inside);
}


#endif // BEZIERLENS_H
