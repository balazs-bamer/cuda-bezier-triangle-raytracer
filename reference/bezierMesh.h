#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_BEZIERMESH

#include "bezierTriangle.h"
#include <functional>

template<typename tReal>
class BezierMesh final {
private:
  std::vector<BezierTriangle<tReal>> mMesh;
  
public:
  using Vector     = typename Mesh<tReal>::Vector;
  using Triangle   = typename Mesh<tReal>::Triangle;
  using Plane      = typename Mesh<tReal>::Plane;
  using Neighbours = typename Mesh<tReal>::Neighbours;

  BezierMesh(Mesh<tReal> const &aMesh);

private:
  void append(Triangle const& aTriangle, std::array<std::reference_wrapper<Triangle const>, 3u> const &aFellows, Neighbours const aNeigh);
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tReal>
BezierMesh<tReal>::BezierMesh(Mesh<tReal> const &aMesh) {
  mMesh.clear();
  mMesh.reserve(aMesh.size() * 3u);
  auto const &mesh                  = aMesh.getMesh();
  auto const &neighbours            = aMesh.getNeighbours();
  auto const &vertex2averageNormals = aMesh.getVertex2averageNormals();
  std::vector<std::array<Vector, 3u>> normalAveragesAtOriginalVertices = getNormalAveragesAtOriginalVertices(aMesh);
  for(uint32_t indexFace = 0u; indexFace < aMesh.size(); ++indexFace) {
    auto const &neigh = neighbours[indexFace];
    auto const &originalTriangle = aMesh[indexFace];
    auto const originalCentral = (originalTriangle[0u] + originalTriangle[1u] + originalTriangle[2u]) / 3.0f;
    auto normal = Mesh<tReal>::getNormal(originalTriangle);
    for(uint32_t indexVertex = 0u; indexVertex < 3u; ++indexVertex) { // Appending a regular triangle means doing Clough-Tocher split on it and appending each one.
      auto const &originalCommonVertex0 = originalTriangle[indexVertex];
      auto const &originalCommonVertex1 = originalTriangle[(indexVertex + 1u) % 3u)];
      auto const &averageNormal0 = vertex2averageNormals[originalCommonVertex0];
      auto const &averageNormal1 = vertex2averageNormals[originalCommonVertex1];
                                                                                 // Neighbour of edge (index, index + 1)
      Plane planeBetweenOriginalNeighbours(normal + Mesh<tReal>::getNormal(aMesh[neigh.mFellowTriangles[indexVertex]]), originalCommonVertex0, originalCommonVertex1);
      auto currentBase = indexFace * 3u;
      std::array<uint32_t, 3u> newNeighbourIndices = { 3u * neigh.mFellowTriangles[i] + neigh.mFellowCommonSideStarts[i], currentBase + (indexVertex + 1u) % 3u, currentBase + (indexVertex + 2u) % 3u) };
      mMesh.emplace_back(BezierTriangle<tReal>(originalCommonVertex0, originalCommonVertex1, originalCentral,
                                               averageNormal0, averageNormal1, planeBetweenOriginalNeighbours),
                                               newNeighbourIndices);
    }
  }
  Vertex originalCentral;
  for(uint32_t i = 0u; i < mMesh.size(); ++i) {
    auto subTriangleIndex = i % 3u;
    auto indexBase = i - subTriangleIndex;
    auto const &triangleNext = mMesh[indexBase + (subTriangleIndex + 1u) % 3u];
    auto const &trianglePrevious = mMesh[indexBase + (subTriangleIndex + 2u) % 3u];
    if(subTriangleIndex % 3u == 0u) {
      auto const &originalTriangle = aMesh[indexBase / 3u];
      originalCentral = (originalTriangle[0u] + originalTriangle[1u] + originalTriangle[2u]) / 3.0f;
    }
    else { // Nothing to do
    }
    mMesh[i].setMissingFields(originalCentral, triangleNext, trianglePrevious);
  }
}

#endif
