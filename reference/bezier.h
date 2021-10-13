#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER

#include "mesh.h"

template<typename tReal>
class BezierTraingle;

template<typename tReal>
using BezierMesh = std::vector<BezierTraingle<tReal>>;

template<typename tReal>
BezierMesh<tReal> transform(Mesh<tReal> const &aMesh) {

};

template<typename tReal>
class BezierTriangle final {    // Cubic Bezier triangle
private:
  using Vector     = Mesh<tReal>::Vector;
  using Vertex     = Mesh<tReal>::Vertex;
  using Matrix     = Mesh<tReal>::Matrix;
  using Triangle   = Mesh<tReal>::Triangle;
  using Neighbours = Mesh<tReal>::Neighbours;

  Vector                  mNormal;                 // normalized
  tReal                   mNormalEquationConstant; // mNormal[0] * x + mNormal[1] * y + mNormal[2] * z = mNormalEquationContant
  std::array<Vector, 3u>  mNeighbourNormals;
  std::array<tReal, 3u>   mNeighbourNormalEquationConstant;
  Neighbours              mNeighbours;             // With new indices after Clough-Tocher split.
  // 0-2 identical to vertices 0-2
  // 3-8 on sides: 3+3i, 4+3i opposite to vertex i, in same order on each side as vertices follow
  // 9 in the middle
  std::array<Vertex, 10u> mControlPoints;
  // Perhaps not needed, because the iterative method to find ray and Bezier surface intersection takes long. std::array<Vertex, 12u> mDerivativeControlPoints; // T most probably more than 2*6 for each partial derivative.
  Matrix                  mBarycentricInverse;      // T = (v1 v2 v3), b = barycentric coefficient column, v = point on plane, v=Tb, b = mBI * v
};

#endif
