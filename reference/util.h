#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_UTIL
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_UTIL

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t

#include <Eigen/Dense>

template<typename tReal>
constexpr tReal gsPi = 3.14159265358979323846;

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
struct Plane final {
  Vector<tReal> mNormal;
  tReal         mConstant;

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

template<typename tReal>
Plane<tReal> Plane<tReal>::createFrom1proportion2points(tReal const aProportion, Vertex<tReal> const &aPoint0, Vertex<tReal> const &aPoint1) {
  Plane result;
  result.mNormal = (aPoint1 - aPoint0).normalized();
  result.mConstant = result.mNormal.dot(aPoint0 * aProportion + aPoint1 * (1.0f - aProportion));
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
  Matrix<tReal> matrix{ aPlane0.mNormal, aPlane1.mNormal, aPlane2.mNormal };
  Vector<tReal> vector{ aPlane0.mConstant, aPlane1.mConstant, aPlane2.mConstant };
  return matrix.inverse() * vector;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif
