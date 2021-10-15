#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_MESHUTILS
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_MESHUTILS

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t

#include "stl_reader.h"
#include <map>
#include <list>
#include <array>
#include <deque>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <Eigen/Dense>


template<typename tReal>
class Mesh final {
private:
  static constexpr tReal csPi                                      = 3.14159265358979323846;
  static constexpr tReal csStandardizeVerticesEpsilonFactor        = 0.2;
  static constexpr tReal csStandardizeNormalsEpsilon               = 0.01;
  static constexpr tReal csStandardizeNormalsIndependentMoveFactor = 0.2;

public:
  using Vector          = Eigen::Matrix<tReal, 3, 1>;
  using Vertex          = Eigen::Matrix<tReal, 3, 1>;
  using Matrix          = Eigen::Matrix<tReal, 3, 3>;
  using Transform       = Eigen::Matrix<tReal, 3, 3>;
  using Triangle        = std::array<Vertex, 3u>;
  using TheMesh         = std::vector<Triangle>;
  using value_type      = Triangle;
  
  struct Neighbours final {
    std::array<uint32_t, 3u> mFellowTriangles;       // Triangle indices.
    std::array<uint8_t, 3u>  mFellowCommonSideStarts; // Vertice indice in each triangle where the common side starts, such that the side in the neighbouring triangle is (i, i+1)
  };
  
  using Face2neighbours = std::vector<Neighbours>;

private:
  TheMesh         mMesh;
  Face2neighbours mFace2neighbours;

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
  void push_back(Triangle const &aTriangle) { mMesh.push_back(aTriangle); }

  TheMesh const&         getMesh() const            { return mMesh; }
  Face2neighbours const& getFace2neighbours() const { return mFace2neighbours; }

// TODO perhaps a function to split triangles having vertex on an edge. Probably not needed.
// TODO later a function to split triangles using Clough-Tocher method if the Bezier Triangle is too high above the triangle.

  void standardizeVertices();

/// Makes all normalvectors (V[1]-V[0])x(V[2]-V[0]) point outwards.
/// Should run after standardizeVertices.
/// But not after auto divideLargeTriangles(tReal const aMaxTriangleSide), because it will probably put vertices on edges of other triangles.
/// Only this method populates the member variable mFace2neighbours.
  void standardizeNormals();

  void transform(Transform const &aTransform, Vertex const aDisplacement);

  Mesh& operator+=(Vector const aDisplacement) {
    transform(Transform::Identity(), aDisplacement);
    return *this;
  }

  Mesh& operator*=(Transform const &aTransform) {
    transform(aTransform, Vertex::Zero());
    return *this;
  }

/// WARNING! Likely to put new vertices on edges!
  void splitTriangles(tReal aMaxTriangleSide);

/// Splits every triangle such that each will become (aDivisor + 1)(aDivisor + 2)/2 new small ones.
  void splitTriangles(int32_t const aDivisor);

  void readMesh(std::string const &aFilename);

  void writeMesh(std::string const &aFilename) const;

  void makeUnitSphere(int32_t const aSectors, int32_t const aBelts);

  static Vector getNormal(Triangle const &aFace) {
    return (aFace[1] - aFace[0]).cross(aFace[2] - aFace[0]);
  }

private:
  tReal getSmallestSide() const;

  using ProjectedIndices = std::multimap<tReal, std::pair<uint32_t, uint32_t>>;
  using Iterator = typename ProjectedIndices::const_iterator;
  using Intervals = std::deque<std::pair<Iterator, Iterator>>;

  void projectVertices(ProjectedIndices &aCurrentProjectedSorted, int32_t const aDimension) const;

  tReal makeProximityIntervals(ProjectedIndices const &aCurrentProjectedSorted, Intervals &aCurrentIntervals, tReal const aEpsilon) const;

  void standardizeInIntervals(Intervals const &aIntervals, tReal const aEpsilonSquared);

  static uint32_t getIndependentFrom(Triangle const &aFaceTarget, Triangle const &aFaceOther);

  static Vector getAltitude(Vertex const &aCommon1, Vertex const &aCommon2, Vertex const &aIndependent);

  struct PairHash {
    template<typename t1, typename t2>
    std::size_t operator()(std::pair<t1, t2> const &aPair) const {
      return std::hash<t1>{}(aPair.first) ^ (std::hash<t2>{}(aPair.second) << 1u);
    }
  };

  struct VertexHash {
    std::size_t operator()(Vertex const &aVertex) const {
      return std::hash<tReal>{}(aVertex[0]) ^ (std::hash<tReal>{}(aVertex[1]) << 1u) ^ (std::hash<tReal>{}(aVertex[2]) << 2u);
    }
  };

  using Edge2face = std::unordered_multimap<std::pair<uint32_t, uint32_t>, uint32_t, PairHash>;
  using Face2vertex = std::vector<std::array<uint32_t, 3u>>;

  std::pair<tReal, uint32_t> createEdge2faceFace2vertex(Edge2face &aEdge2face, Face2vertex &aFace2vertex) const;

  void createFace2neighbourFacesAtSmallestX(Edge2face const &aEdge2face, Face2vertex const &aFace2vertex, uint32_t const aSmallestXverticeIndex,
                                            Face2neighbours &aFace2neighbour, std::unordered_set<uint32_t> &aFacesAtSmallestX) const;

  uint32_t getInitialFaceIndex(std::unordered_set<uint32_t> const &aFacesAtSmallestX, Vertex const &aDesiredVector) const;

  static void normalize(Triangle &aFace, Vertex const &aDesiredVector);

  static void normalize(Triangle const &aFaceKnown, Triangle &aFaceUnknown);

  static void divideTriangle(TheMesh &aMesh, Triangle const &aTriangle, int32_t const aDivisor);
};

/////////////////////////////////
//       IMPLEMENTATION        //
/////////////////////////////////

template<typename tReal>
tReal Mesh<tReal>::getSmallestSide() const {
  tReal smallestSide = std::numeric_limits<tReal>::max();
  for(auto &triangle : mMesh) {
    for(uint32_t i = 0u; i < 3u; ++i) {
      smallestSide = std::min(smallestSide, (triangle[i] - triangle[(i + 1u) % 3u]).norm());
    }
  }
  return smallestSide;
}

template<typename tReal>
void Mesh<tReal>::projectVertices(ProjectedIndices &aCurrentProjectedSorted, int32_t const aDimension) const {
  for(uint32_t indexFace = 0u; indexFace < mMesh.size(); ++indexFace) {
    auto &face = mMesh[indexFace];
    for(uint32_t indexVertex = 0u; indexVertex < 3u; ++indexVertex) {
      auto &vertex = face[indexVertex];
      aCurrentProjectedSorted.emplace(std::make_pair(vertex[aDimension], std::make_pair(indexFace, indexVertex)));
    }
  }
}

template<typename tReal>
tReal Mesh<tReal>::makeProximityIntervals(ProjectedIndices const &aCurrentProjectedSorted, Intervals &aCurrentIntervals, tReal const aEpsilon) const {
  bool was = false;
  tReal startValue;
  uint32_t max = 0u, counter;
  Iterator startIterator = aCurrentProjectedSorted.cend();;
  for(Iterator i = aCurrentProjectedSorted.cbegin(); i != aCurrentProjectedSorted.cend(); ++i) {
    if(was) {
      auto valueNow = i->first;
      if(valueNow - startValue >= aEpsilon) {
        aCurrentIntervals.emplace_back(std::make_pair(startIterator, i));
        startValue = valueNow;
        startIterator = i;
        max = std::max(max, counter);
        counter = 1u;
      }
      else {
        ++counter;
      }
    }
    else {
      startIterator = i;
      startValue = i->first;
      was = true;
      counter = 1u;
    }
  }
  Iterator last = aCurrentProjectedSorted.cend();
  aCurrentIntervals.emplace_back(std::make_pair(startIterator, last));
  max = std::max(max, counter);
  return max;
}

template<typename tReal>
void Mesh<tReal>::standardizeInIntervals(Intervals const &aIntervals, tReal const aEpsilonSquared) {
  for(auto const& interval : aIntervals) {                               // We need to check pairwise proximity only inside each interval.
    for(auto i = interval.first; i != interval.second; ++i) {
      auto &v1 = mMesh[i->second.first][i->second.second];
      for(auto j = interval.first; j != interval.second; ++j) {
        auto &v2 = mMesh[j->second.first][j->second.second];
        if((v1-v2).squaredNorm() < aEpsilonSquared && (v1[0] < v2[0] || (v1[0] == v2[0] && v1[1] < v2[1]) || (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] < v2[2]))) {
          v1 = v2;                                                        // Use lexicographic sorting to choose an unique one among the points too close to each other.
        }
        else { // nothing to do
        }
      }
    }
  }
}

template<typename tReal>
void Mesh<tReal>::standardizeVertices() {
  auto epsilon = getSmallestSide() * csStandardizeVerticesEpsilonFactor; // Vertices closer to each other than this will be forced to one point.

  std::array<uint32_t, 3> maximums;
  std::array<ProjectedIndices, 3u> projectedSorted;                           // Holds the projected indices sorted by a coordinate, separate for x, y and z.
  std::array<Intervals, 3u> proximityIntervals;                               // Holds iterator pairs to the previous stuff in between the projections are closer to each other than epsilon.

  for(uint32_t dimension = 0u; dimension < 3u; ++dimension) {
    auto& currentProjectedSorted = projectedSorted[dimension];
    auto& currentIntervals = proximityIntervals[dimension];

    projectVertices(currentProjectedSorted, dimension);
    maximums[dimension] = makeProximityIntervals(currentProjectedSorted, currentIntervals, epsilon);
  }
  auto minIt = std::min_element(maximums.cbegin(), maximums.cend());          // Best coordinate is where the maximum interval length is minimal.
  auto whichDim = minIt - maximums.cbegin();
  auto const& bestIntervals = proximityIntervals[whichDim];

  standardizeInIntervals(bestIntervals, epsilon * epsilon);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tReal>
uint32_t Mesh<tReal>::getIndependentFrom(Triangle const &aFaceTarget, Triangle const &aFaceOther) {
  for(uint32_t i = 0u; i < 3u; ++i) {
    if(std::find(aFaceOther.cbegin(), aFaceOther.cend(), aFaceTarget[i]) == aFaceOther.cend()) {
      return i;
    }
    else { // Nothing to do
    }
  }
  return 3u;                                                          // Won't get here, but avoid compiler complaint.
}

template<typename tReal>
typename Mesh<tReal>::Vector Mesh<tReal>::getAltitude(Vertex const &aCommon1, Vertex const &aCommon2, Vertex const &aIndependent) {
  auto commonVector = aCommon2 - aCommon1;
  auto independentVector = aIndependent - aCommon1;
  auto footFactor = commonVector.dot(independentVector) / commonVector.squaredNorm();
  return independentVector - commonVector * footFactor; 
}

template<typename tReal>
std::pair<tReal, uint32_t> Mesh<tReal>::createEdge2faceFace2vertex(Edge2face &aEdge2face, Face2vertex &aFace2vertex) const {
  std::unordered_map<Vertex, uint32_t, VertexHash> was;               // For presence testing, maps each vertex to its index in the deque below.
  tReal smallestX = std::numeric_limits<tReal>::max();
  uint32_t smallestXverticeIndex;
  for(uint32_t indexFace = 0u; indexFace < mMesh.size(); ++indexFace) {
    auto const &face = mMesh[indexFace];
    std::array<uint32_t, 3u> faceVertexIndices;                       // We collect here the vertex indices for the actual face.
    for(uint32_t indexInFace = 0u; indexInFace < 3u; ++indexInFace) {
      auto const vertex = face[indexInFace];
      uint32_t indexInVertices;
      auto found = was.find(vertex);
      if(found == was.end()) {
        indexInVertices = was.size();
        faceVertexIndices[indexInFace] = indexInVertices;
        was.emplace(std::make_pair(vertex, indexInVertices));
      }
      else {
        indexInVertices = found->second;
        faceVertexIndices[indexInFace] = indexInVertices;
      }
      if(vertex[0] < smallestX) {
        smallestX = vertex[0];
        smallestXverticeIndex = indexInVertices;
      }
      else { // Nothing to do
      }
    }
    aFace2vertex.push_back(faceVertexIndices);
    for(uint32_t i = 0u; i < 3u; ++i) {
      auto low = faceVertexIndices[i];
      auto high = faceVertexIndices[(i + 1u) % 3u];
      if(low > high) {
        std::swap(low, high);
      }
      else { // Nothing to do
      }
      aEdge2face.emplace(std::make_pair(std::pair(low, high), indexFace));
    }
  }
  return std::make_pair(smallestX, smallestXverticeIndex);
}

template<typename tReal>
void Mesh<tReal>::createFace2neighbourFacesAtSmallestX(Edge2face const &aEdge2face,
                                                       Face2vertex const &aFace2vertex,
                                                       uint32_t const aSmallestXverticeIndex,
                                                       Face2neighbours &aFace2neighbour,
                                                       std::unordered_set<uint32_t> &aFacesAtSmallestX) const {
  for(uint32_t indexFace = 0u; indexFace < mMesh.size(); ++indexFace) {
    auto const &face = aFace2vertex[indexFace];
    Neighbours neighbours;
    for(uint32_t indexInFace = 0u; indexInFace < 3u; ++indexInFace) {
      if(face[indexInFace] == aSmallestXverticeIndex) {
        aFacesAtSmallestX.insert(indexFace);
      }
      else { // Nothing to do
      }
      auto otherVertexIndex = face[(indexInFace + 1u) % 3u];
      auto edge = std::make_pair(face[indexInFace], otherVertexIndex);
      if(edge.first > edge.second) {
        std::swap(edge.first, edge.second);
      }
      else { // Nothing to do
      }
      auto[begin, end] = aEdge2face.equal_range(edge);
      if(begin->second == indexFace) {                                // We need the other, neighbouring face.
        ++begin;
        if(begin == end) {
          throw "Vertex on edge detected.";
        }
        else { // Nothing to do
        }
      }
      else { // Nothing to do
      }
      auto otherIndexFace = begin->second;
      neighbours.mFellowTriangles[indexInFace] = otherIndexFace;      // Neighbour of edge (index, index + 1)
      auto const &otherFace = aFace2vertex[otherIndexFace];
      auto const otherIndexInFace = std::find(otherFace.cbegin(), otherFace.cend(), otherVertexIndex) - otherFace.cbegin(); 
      neighbours.mFellowCommonSideStarts[indexInFace] = otherIndexInFace;
    }
    aFace2neighbour.push_back(neighbours);
  }
}

template<typename tReal>
uint32_t Mesh<tReal>::getInitialFaceIndex(std::unordered_set<uint32_t> const &aFacesAtSmallestX, Vertex const &aDesiredVector) const {
  tReal maxAbsoluteDotProduct = -std::numeric_limits<tReal>::max();
  uint32_t initialFaceIndex;                                          // First determine the initial face for which we want (f[1]-f[0])x(f[2]-f[0]) point outwards.
  for(auto const indexFace : aFacesAtSmallestX) {
    auto const &face = mMesh[indexFace];
    auto normal = getNormal(face).normalized();
    auto absoluteDotProduct = std::abs(aDesiredVector.dot(normal));
    if(absoluteDotProduct > maxAbsoluteDotProduct) {
      maxAbsoluteDotProduct = absoluteDotProduct;
      initialFaceIndex = indexFace;
    }
    else { // Nothing to do
    }
  }
  return initialFaceIndex;
}

template<typename tReal>
void Mesh<tReal>::normalize(Triangle &aFace, Vertex const &aDesiredVector) {
  auto normal = getNormal(aFace);
  if(aDesiredVector.dot(normal) < 0.0f) {
    std::swap(aFace[0u], aFace[1u]);
  }
  else { // Nothing to do
  }
}

template<typename tReal>
void Mesh<tReal>::normalize(Triangle const &aFaceKnown, Triangle &aFaceUnknown) {  // Calculates everything twice for each triangle, but won't cache now.
  uint32_t independentIndexKnown = getIndependentFrom(aFaceKnown, aFaceUnknown);
  uint32_t commonIndex1known = (independentIndexKnown + 1u) % 3u;
  uint32_t commonIndex2known = (independentIndexKnown + 2u) % 3u;
  uint32_t independentIndexUnknown = getIndependentFrom(aFaceUnknown, aFaceKnown);
  uint32_t commonIndex1unknown = (independentIndexUnknown + 1u) % 3u;
  uint32_t commonIndex2unknown = (independentIndexUnknown + 2u) % 3u;

  auto altitudeKnown = getAltitude(aFaceKnown[commonIndex1known], aFaceKnown[commonIndex2known], aFaceKnown[independentIndexKnown]);
  Vertex altitudeUnknown = getAltitude(aFaceUnknown[commonIndex1unknown], aFaceUnknown[commonIndex2unknown], aFaceUnknown[independentIndexUnknown]);
  auto dotFaceAltitudes = altitudeKnown.dot(altitudeUnknown);
  
  auto normalKnown = getNormal(aFaceKnown);
  auto normalUnknown = getNormal(aFaceUnknown);
  auto knownDotUnknown = normalKnown.dot(normalUnknown);
  if(std::abs(knownDotUnknown / (normalKnown.norm() * normalUnknown.norm())) < csStandardizeNormalsEpsilon) {
    auto newIndependentUnknown = aFaceUnknown[independentIndexUnknown] + csStandardizeNormalsIndependentMoveFactor * (aFaceKnown[independentIndexKnown] - (aFaceKnown[commonIndex1known] + aFaceKnown[commonIndex2known]) / 2.0f);
    auto newFaceUnknown = aFaceUnknown;
    newFaceUnknown[independentIndexUnknown] = newIndependentUnknown;

    altitudeUnknown = getAltitude(aFaceUnknown[commonIndex1unknown], aFaceUnknown[commonIndex2unknown], newIndependentUnknown);
    dotFaceAltitudes = altitudeKnown.dot(altitudeUnknown);
    normalUnknown = getNormal(newFaceUnknown);
    knownDotUnknown = normalKnown.dot(normalUnknown);
  }
  else { // Nothing to do
  }
  if(dotFaceAltitudes * knownDotUnknown > 0.0f) {
    std::swap(aFaceUnknown[commonIndex1unknown], aFaceUnknown[commonIndex2unknown]);
  }
  else { // Nothing to do
  }
}

template<typename tReal>
void Mesh<tReal>::standardizeNormals() {
  Edge2face edge2face;                                                // Edge is identified by its two ordered vertex indices.
  Face2vertex face2vertex;
  face2vertex.reserve(mMesh.size());
  mFace2neighbours.clear();
  mFace2neighbours.reserve(mMesh.size());

  auto [smallestX, smallestXverticeIndex] = createEdge2faceFace2vertex(edge2face, face2vertex);

  std::unordered_set<uint32_t> facesAtSmallestX;
  
  createFace2neighbourFacesAtSmallestX(edge2face, face2vertex, smallestXverticeIndex,
                                       mFace2neighbours, facesAtSmallestX);

  Vertex desiredVector;                                               // A known unit vector pointing outwards for a definite face.
  desiredVector << -1.0f, 0.0f, 0.0f;
  uint32_t initialFaceIndex = getInitialFaceIndex(facesAtSmallestX, desiredVector);

  struct KnownUnknownFacePair {
    uint32_t mKnownIndex;
    uint32_t mUnknownIndex;
  };
  std::list<KnownUnknownFacePair> queue;                              // Faces with normals yet to be normalized.
  normalize(mMesh[initialFaceIndex], desiredVector);
  std::vector<bool> remaining;
  remaining.reserve(mFace2neighbours.size());
  std::fill_n(std::back_inserter(remaining), mFace2neighbours.size(), true);
  for(auto const indexFace : mFace2neighbours[initialFaceIndex].mFellowTriangles) {
    queue.emplace_back(KnownUnknownFacePair{initialFaceIndex, indexFace});
  }
  remaining[initialFaceIndex] = false;
  while(!queue.empty()) {
    auto actualPair = queue.back();
    queue.pop_back();
    if(remaining[actualPair.mUnknownIndex]) {
      normalize(mMesh[actualPair.mKnownIndex], mMesh[actualPair.mUnknownIndex]);
    }
    else { // Nothing to do
    }
    remaining[actualPair.mUnknownIndex] = false;
    for(auto const indexFace : mFace2neighbours[actualPair.mUnknownIndex].mFellowTriangles) {
      if(remaining[indexFace] && actualPair.mUnknownIndex != indexFace) {
        queue.emplace_back(KnownUnknownFacePair{actualPair.mUnknownIndex, indexFace});
      }
      else { // Nothing to do
      }
    }
  }
  calculateNormalAverages4vertices();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tReal>
void Mesh<tReal>::transform(Transform const &aTransform, Vertex const aDisplacement) {
  for(auto &triangle : mMesh) {
    for(auto &vertex : triangle) {
      vertex = aTransform * vertex + aDisplacement;      
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tReal>
void Mesh<tReal>::divideTriangle(TheMesh &result, Triangle const &aTriangle, int32_t const aDivisor) {
  auto vector01 = (aTriangle[1] - aTriangle[0]) / aDivisor;
  auto vector02 = (aTriangle[2] - aTriangle[0]) / aDivisor;
  auto lineBase = aTriangle[0];
  auto base0 = lineBase;
  auto base1 = (aDivisor > 1) ? (base0 + vector01) : aTriangle[1];
  auto base2 = (aDivisor > 1) ? (base0 + vector02) : aTriangle[2];
  for(int32_t i = 0; i < aDivisor - 1; ++i) {
    for(int32_t j = 0; j < aDivisor - i - 1; j++) {
      result.push_back({base0, base1, base2});
      auto base1next = base1 + vector02;
      result.push_back({base1, base1next, base2});
      base1 = base1next;
      base0 = base2;
      base2 += vector02;
    }
    result.push_back({base0, base1, base2});
    lineBase += vector01;
    base0 = lineBase;
    base1 = base0 + vector01;
    base2 = base0 + vector02;
  }
  result.push_back({base0, aTriangle[1], base2});
}

template<typename tReal>
void Mesh<tReal>::splitTriangles(tReal const aMaxTriangleSide) {
  TheMesh result;
  for(auto const &triangle : mMesh) {
    tReal maxSide = (triangle[0] - triangle[1]).norm();
    maxSide = std::max(maxSide, (triangle[0] - triangle[2]).norm());
    maxSide = std::max(maxSide, (triangle[1] - triangle[2]).norm());
    int32_t divisor = static_cast<int32_t>(std::ceil(maxSide / aMaxTriangleSide));
    divideTriangle(result, triangle, divisor);
  }
  mMesh = result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tReal>
void Mesh<tReal>::splitTriangles(int32_t const aDivisor) {
  TheMesh result;
  for(auto const &triangle : mMesh) {
    divideTriangle(result, triangle, aDivisor);
  }
  mMesh = result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tReal>
void Mesh<tReal>::readMesh(std::string const &aFilename) {
  mMesh.clear();
  mFace2neighbours.clear();
  stl_reader::StlMesh<tReal, int32_t> mesh(aFilename);
  for(int32_t indexTriangle = 0; indexTriangle < mesh.num_tris(); ++indexTriangle) {
    Triangle triangle;
    for(int32_t indexCorner = 0; indexCorner < 3; ++indexCorner) {
      tReal const * const coords = mesh.tri_corner_coords(indexTriangle, indexCorner);
      Vertex in;
      for(int32_t i = 0; i < 3; ++i) {
        in(i) = coords[i];
      }
      triangle[indexCorner] = in;
    }
    mMesh.push_back(triangle);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tReal>
void Mesh<tReal>::writeMesh(std::string const& aFilename) const {
  std::ofstream out(aFilename);
  out << "solid Exported from Blender-2.82 (sub 7)\n";
  for(auto const & triangle : mMesh) {
    out << "facet normal 0.000000 0.000000 0.000000\nouter loop\n";
    for(auto const & vertex : triangle) {
      out << "vertex " << vertex(0) << ' ' << vertex(1) << ' ' << vertex(2) << '\n';
    }
    out << "endloop\nendfacet\n";
  }
  out << "endsolid Exported from Blender-2.82 (sub 7)\n";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename tReal>
void Mesh<tReal>::makeUnitSphere(int32_t const aSectors, int32_t const aBelts) {
  mMesh.clear();
  mMesh.reserve(static_cast<uint32_t>(2 * aSectors * aBelts));
  mFace2neighbours.clear();
  tReal sectorAngleHalf = csPi / aSectors;
  tReal sectorAngleFull = sectorAngleHalf * 2.0f;
  tReal beltAngle       = csPi / (aBelts + 1.0f);
  tReal bias = 0.0f;
  tReal beltAngleUp = 0.0f;
  tReal beltAngleMiddle = beltAngle;
  tReal beltAngleDown = 2.0f * beltAngle;
  tReal beltRadiusUp = 0.0f;
  tReal beltRadiusMiddle = std::sin(beltAngleMiddle);
  tReal beltRadiusDown = std::sin(beltAngleDown);
  tReal beltZup = 1.0f;
  tReal beltZmiddle = std::cos(beltAngleMiddle);
  tReal beltZdown = std::cos(beltAngleDown);
  for(int32_t belt = 0; belt < aBelts; ++belt) {
    tReal sectorAngleUpDown = bias + sectorAngleHalf;
    tReal sectorAngleMiddle1 = bias + 0.0f;
    tReal sectorAngleMiddle2 = bias + sectorAngleFull;
    for(int32_t sector = 0; sector < aSectors; ++sector) {
      Vertex corner1(beltRadiusUp * std::sin(sectorAngleUpDown), beltRadiusUp * std::cos(sectorAngleUpDown), beltZup);
      Vertex corner2(beltRadiusMiddle * std::sin(sectorAngleMiddle1), beltRadiusMiddle * std::cos(sectorAngleMiddle1), beltZmiddle);
      Vertex corner3(beltRadiusMiddle * std::sin(sectorAngleMiddle2), beltRadiusMiddle * std::cos(sectorAngleMiddle2), beltZmiddle);
      mMesh.push_back({corner1, corner2, corner3});
      corner1 = {beltRadiusDown * std::sin(sectorAngleUpDown), beltRadiusDown * std::cos(sectorAngleUpDown), beltZdown};
      mMesh.push_back({corner2, corner3, corner1});
      sectorAngleUpDown += sectorAngleFull;
      sectorAngleMiddle1 = sectorAngleMiddle2;
      sectorAngleMiddle2 += sectorAngleFull;
    }
    beltAngleUp = beltAngleMiddle;
    beltAngleMiddle = beltAngleDown;
    beltAngleDown += beltAngle;
    beltRadiusUp = beltRadiusMiddle;
    beltRadiusMiddle = beltRadiusDown;
    beltRadiusDown = std::sin(beltAngleDown);
    beltZup = beltZmiddle;
    beltZmiddle = beltZdown;
    beltZdown = std::cos(beltAngleDown);
    bias += sectorAngleHalf;
  }
}

#endif
