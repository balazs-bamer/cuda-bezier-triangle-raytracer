#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_MESH
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_MESH

#include "util.h"
#include "stl_reader.h"
#include <map>
#include <list>
#include <array>
#include <deque>
#include <vector>
#include <string>
#include <fstream>
#include <functional>
#include <unordered_map>
#include <unordered_set>


class Mesh final {
private:
  static constexpr float csStandardizeVerticesEpsilonFactor        = 0.2;
  static constexpr float csStandardizeNormalsEpsilon               = 0.01;
  static constexpr float csStandardizeNormalsIndependentMoveFactor = 0.2;

public:
  using TheMesh         = std::vector<Triangle>;
  using value_type      = Triangle;
  
  struct Neighbours final {
    std::array<uint32_t, 3u> mFellowTriangles;        // Triangle indices: neighbour of edge (index, index + 1)
    std::array<uint8_t, 3u>  mFellowCommonSideStarts; // Vertice index in each fellow triangle where the common side starts,
                                                      // such that the side in the neighbouring triangle is (i, i+1)
  };

  struct VertexHash {
    std::size_t operator()(Vertex const &aVertex) const {
      return std::hash<float>{}(aVertex[0]) ^ (std::hash<float>{}(aVertex[1]) << 1u) ^ (std::hash<float>{}(aVertex[2]) << 2u);
    }
  };

  using Face2neighbours       = std::vector<Neighbours>;
  using Vertex2averageNormals = std::unordered_map<Vertex, Vector, VertexHash>;
  using Vertices              = std::unordered_set<Vertex, VertexHash>;

private:
  TheMesh               mMesh;
  Face2neighbours       mFace2neighbours;
  Vertex2averageNormals mVertex2averageNormals;

public:
  Mesh() = default;
  Mesh(Mesh &&) = default;
  Mesh(Mesh const&) = default;
  Mesh& operator=(Mesh &&) = default;
  Mesh& operator=(Mesh const&) = default;

  auto size() const { return mMesh.size(); }
  auto begin() const { return mMesh.begin(); }
  auto end() const { return mMesh.end(); }
  auto cbegin() const { return mMesh.cbegin(); }
  auto cend() const { return mMesh.cend(); }
  void reserve(uint32_t const aSize) { mMesh.reserve(aSize); }
  void push_back(Triangle const &aTriangle) { mMesh.push_back(aTriangle); }
  void clear() { mMesh.clear(); }
  auto& operator[](uint32_t const aI) { return mMesh[aI]; }
  auto const &operator[](uint32_t const aI) const { return mMesh[aI]; }

  TheMesh const&               getMesh() const                  { return mMesh; }
  Face2neighbours const&       getFace2neighbours() const       { return mFace2neighbours; }
  Vertex2averageNormals const& getVertex2averageNormals() const { return mVertex2averageNormals; } // TODO perhaps implement a way to import these values
                                                                                                   // along with the triangles when the original surface is known
                                                                                                   // and dervatives are present.
  void standardizeVertices();
  Vertices getVertices() const;

/// Makes all normalvectors (V[1]-V[0])x(V[2]-V[0]) point outwards.
/// Should run after standardizeVertices.
/// But not after auto divideLargeTriangles(float const aMaxTriangleSide), because it will probably put vertices on edges of other triangles.
/// Only this method populates the member variable mFace2neighbours and mVertex2averageNormals
  void standardizeNormals();

  void transform(Transform const &aTransform, Vertex const aDisplacement);

  Mesh& operator+=(Vector const aDisplacement)  { transform(Transform::Identity(),           aDisplacement);  return *this; }
  Mesh& operator*=(Transform const &aTransform) { transform(aTransform,                      Vertex::Zero()); return *this; }
  Mesh& operator*=(float const &aFactor)        { transform(Transform::Identity() * aFactor, Vertex::Zero()); return *this; }

/// WARNING! Likely to put new vertices on edges!
  void splitTriangles(float const aMaxTriangleSide);

/// Splits every triangle such that each will become (aDivisor + 1)(aDivisor + 2)/2 new small ones.
  void splitTriangles(int32_t const aDivisor);

  void readMesh(std::string const &aFilename);

  void writeMesh(std::string const &aFilename) const;

  /// aEnvelope is a C2 function in [-1, 1], with prefereably vertical tangents in -1 and 1 and values of 0 in -1 and 1.
  void makeSolidOfRevolution(int32_t const aSectors, int32_t const aBelts, std::function<float(float)> aEnvelope, Vector const &aSize);
  void makeEllipsoid(int32_t const aSectors, int32_t const aBelts, Vector const &aSize) { makeSolidOfRevolution(aSectors, aBelts, [](float const aX){ return ::sqrt(1 - aX * aX); }, aSize); }
  void makeUnitSphere(int32_t const aSectors, int32_t const aBelts) { makeEllipsoid(aSectors, aBelts, Vector(1.0f, 1.0f, 1.0f)); }

private:
  float getSmallestSide() const;

  using ProjectedIndices = std::multimap<float, std::pair<uint32_t, uint32_t>>;
  using Iterator = typename ProjectedIndices::const_iterator;
  using Intervals = std::deque<std::pair<Iterator, Iterator>>;

  void projectVertices(ProjectedIndices &aCurrentProjectedSorted, int32_t const aDimension) const;
  float makeProximityIntervals(ProjectedIndices const &aCurrentProjectedSorted, Intervals &aCurrentIntervals, float const aEpsilon) const;
  void standardizeInIntervals(Intervals const &aIntervals, float const aEpsilonSquared);
  static uint32_t getIndependentFrom(Triangle const &aFaceTarget, Triangle const &aFaceOther);

  struct PairHash {
    template<typename t1, typename t2>
    std::size_t operator()(std::pair<t1, t2> const &aPair) const {
      return std::hash<t1>{}(aPair.first) ^ (std::hash<t2>{}(aPair.second) << 1u);
    }
  };

  using Edge2face = std::unordered_multimap<std::pair<uint32_t, uint32_t>, uint32_t, PairHash>;
  using Face2vertex = std::vector<std::array<uint32_t, 3u>>;
  using Vertex2index = std::unordered_map<Vertex, uint32_t, VertexHash>;

  std::tuple<Edge2face, Face2vertex, Vertex2index> createEdge2faceFace2vertex() const;
  std::pair<float, std::unordered_set<uint32_t>> getSmallestXstuff(Face2vertex const &aFace2vertex, Vertex2index const &) const;
  void createFace2neighbour(Edge2face const &aEdge2face, Face2vertex const &aFace2vertex);
  uint32_t getInitialFaceIndex(std::unordered_set<uint32_t> const &aFacesAtSmallestX, Vertex const &aDesiredVector) const;
  static void normalize(Triangle &aFace, Vertex const &aDesiredVector);
  static void normalize(Triangle const &aFaceKnown, Triangle &aFaceUnknown);
  void calculateNormalAverages4vertices();
  static void divideTriangle(TheMesh &aMesh, Triangle const &aTriangle, int32_t const aDivisor);
};

#endif
