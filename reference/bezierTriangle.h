#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIER

#include "mesh.h"

template<typename tReal>
class BezierTriangle final {    // Cubic Bezier triangle
private:
  using Vector     = typename Mesh<tReal>::Vector;
  using Vertex     = typename Mesh<tReal>::Vertex;
  using Matrix     = typename Mesh<tReal>::Matrix;
  using Triangle   = typename Mesh<tReal>::Triangle;

  Vector                   mNormal;                  // normalized
  tReal                    mNormalEquationConstant;  // mNormal[0] * x + mNormal[1] * y + mNormal[2] * z = mNormalEquationContant
  std::array<Vector, 3u>   mNeighbourDividerNormals; // These are about the plane going through each edge and having the normal the average of the adjacent side normals.
  std::array<tReal, 3u>    mNeighbourDividerNormalEquationConstant;
  std::array<uint32_t, 3u> mNeighbours;              // With new indices after Clough-Tocher split.
  // 0-2 identical to vertices 0-2
  // 3-8 on sides: 3+3i, 4+3i opposite to vertex i, in same order on each side as vertices follow
  // 9 in the middle
  std::array<Vertex, 10u>  mControlPoints;
  // Perhaps not needed, because the iterative method to find ray and Bezier surface intersection takes long. std::array<Vertex, 12u> mDerivativeControlPoints; // T most probably more than 2*6 for each partial derivative.
  Matrix                   mBarycentricInverse;       // T = (v1 v2 v3), b = barycentric coefficient column, v = point on plane, v=Tb, b = mBI * v

public:
  // Vertices in argument are in order such that the normal points to the desired direction.
  // Neighbour i is Neighbour of edge (i, i + 1)
  BezierTriangle(Vertex const &aVertex0, Vertex const &aVertex1, Vertex const &aVertex2,
		 uint32_t const aNeigh0, uint32_t const aNeigh1, uint32_t const aNeigh2);
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tReal>
BezierTriangle<tReal>::BezierTriangle(Vertex const &aVertex0, Vertex const &aVertex1, Vertex const &aVertex2,
		                             uint32_t const aNeigh0, uint32_t const aNeigh1, uint32_t const aNeigh2) {
  mControlPoints[0u] = aVertex0;   
  mControlPoints[1u] = aVertex1;   
  mControlPoints[2u] = aVertex2;
  mNeighbours[0u] = aNeigh0;
  mNeighbours[1u] = aNeigh0;
  mNeighbours[2u] = aNeigh0;
  mNormal = Mesh<tReal>::getNormal({aVertex0, aVertex1, aVertex2}).normalized();
  mNormalEquationConstant = mNormal.dot(aVertex0);
  Matrix<tReal> vertices;
  vertices.column(0) = aVertex0; 
  vertices.column(1) = aVertex1; 
  vertices.column(2) = aVertex2; 
  mBarycentricInverse = vertices.inverse();
}
#endif
