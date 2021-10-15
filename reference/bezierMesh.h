#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH

#include "bezierTriangle.h"
#include <functional>

template<typename tReal>
class BezierMesh final {
private:
  std::vector<BezierTriangle<tReal>> mMesh;
  
public:
  using Triangle   = typename Mesh<tReal>::Triangle;
  using Neighbours = typename Mesh<tReal>::Neighbours;

  BezierMesh(Mesh<tReal> const &aMesh) {
    mMesh.clear();
    mMesh.reserve(aMesh.size() * 3u);
    auto const &mesh = aMesh.getMesh();
    auto const &neighbours = aMesh.getNeighbours();
    for(uint32_t i = 0u; i < aMesh.size(); ++i) {
      auto const &neigh = neighbours[i];
      std::array<std::reference_wrapper<Triangle const>, 3u> fellows = { mesh[neigh.mFellowTriangles[0u]], 
	                                                                 mesh[neigh.mFellowTriangles[1u]], 
									 mesh[neigh.mFellowTriangles[2u]] };
      append(mesh[i], fellows, neigh);
    }
  }

private:
  void append(Triangle const& aTriangle, std::array<std::reference_wrapper<Triangle const>, 3u> const &aFellows, Neighbours const aNeigh);
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tReal>
void BezierMesh<tReal>::append(Triangle const& aTriangle, std::array<std::reference_wrapper<Triangle const>, 3u> const &aFellows, Neighbours const aNeigh) {
  auto center = (aTriangle[0] + aTriangle[1] + aTriangle[2]) / 3.0f;
  auto currentBase = mMesh.size();
  for(uint32_t i = 0u; i < 3u; ++i) {
    mMesh.emplace_back(BezierTriangle<tReal>(aTriangle[i], aTriangle[(i + 1u) % 3u], center,
  	               3u * aNeigh.mFellowTriangles[i] + aNeigh.mFellowCommonSideStarts[i], currentBase + (i + 1u) % 3u, currentBase + (i + 2u) % 3u)); 
  }
  std::array<std::reference_wrapper<BezierTriangle<tReal>>, 3u> triangles = { mMesh[currentBase], 
                                                                              mMesh[currentBase + 1u], 
                                                                              mMesh[currentBase + 2u] };
}

#endif
