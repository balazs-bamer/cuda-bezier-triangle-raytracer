#include "util.h"
#include "mesh.h"
#include "bezierMesh.h"

#include<deque>
#include<vector>
#include<iostream>

#include "gtest/gtest.h"

using Real = float;

constexpr Real cgEpsilon = 0.0001f;

std::ostream& operator<<(std::ostream &aOut, Plane<Real> const &aPlane) {
  aOut << aPlane.mNormal(0) << ' ' << aPlane.mNormal(1) << ' ' << aPlane.mNormal(2) << ' ' << " : " << aPlane.mConstant;
  return aOut;
}


Vertex<Real> getPlaneIntersectionNormals(Vertex<Real> aPoint,
                                         Vector<Real> aDir1,
                                         Vector<Real> aDir2,
                                         Vector<Real> aDir3) {
  Vector<Real> normal1 = aDir1.normalized();
  Plane<Real> plane1(normal1, normal1.dot(aPoint));
  Vector<Real> normal2 = aDir2.normalized();
  Plane<Real> plane2(normal2, normal2.dot(aPoint));
  Vector<Real> normal3 = aDir3.normalized();
  Plane<Real> plane3(normal3, normal3.dot(aPoint));
  return Plane<Real>::intersect(plane1, plane2, plane3);
}

TEST(vector, getAperpendicular) {
  static constexpr Real csEpsilon = 1e-10;
  {
    Vector<Real> vector{1.0f, 0.0f, 0.0f};
    EXPECT_LT(::abs(getAperpendicular(vector).dot(vector)), csEpsilon);
  }
  {
    Vector<Real> vector{1.0f, 1.0f, 0.0f};
    vector.normalize();
    EXPECT_LT(::abs(getAperpendicular(vector).dot(vector)), csEpsilon);
  }
  {
    Vector<Real> vector{1.0f, 0.0f, 1.0f};
    vector.normalize();
    EXPECT_LT(::abs(getAperpendicular(vector).dot(vector)), csEpsilon);
  }
  {
    Vector<Real> vector{1.0f, -1.0f, -1.0f};
    vector.normalize();
    EXPECT_LT(::abs(getAperpendicular(vector).dot(vector)), csEpsilon);
  }
}

TEST(ray, averageErrorSquared) {
  {
    Ray<Real> ray({0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f});
    std::vector<Vertex<Real>> points{};
    EXPECT_EQ(ray.getAverageErrorSquared(points), 0.0f);
  }
  {
    Ray<Real> ray({0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f});
    std::vector<Vertex<Real>> points{{2.0f, 0.0f, 0.0f}, {-3.0f, 0.0f, 0.0f}};
    EXPECT_EQ(ray.getAverageErrorSquared(points), 0.0f);
  }
  {
    Ray<Real> ray({0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f});
    std::vector<Vertex<Real>> points{{2.0f, 1.0f, 0.0f}, {-3.0f, 0.0f, 1.0f}};
    EXPECT_GT(ray.getAverageErrorSquared(points), 0.0f);
  }
}

TEST(planeIntersection, Normals) {
  {
    Vertex<Real> common{1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionNormals(common, {1.0, 2.0f, 3.0f}, {3.0f, 1.0f, 2.0f}, {3.0f, 2.0f, 1.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{3.0f, -2.0f, 1.0f};
    auto result = getPlaneIntersectionNormals(common, {1.0, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{3.0f, -2.0f, -1.0f};
    auto result = getPlaneIntersectionNormals(common, {1.0, -2.0f, 3.0f}, {-1.0f, 2.0f, 3.0f}, {1.0f, 2.0f, -3.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
}

Vertex<Real> getPlaneIntersectionProportionPoints(Vertex<Real> aCommon,
                                                  Real aProportion1, Vertex<Real> aOne1,
                                                  Real aProportion2, Vertex<Real> aOne2,
                                                  Real aProportion3, Vertex<Real> aOne3) {
  Plane<Real> plane1 = Plane<Real>::createFrom1proportion2points(aProportion1, aOne1, aOne1 + 1.0f / aProportion1 * (aCommon - aOne1));
  Plane<Real> plane2 = Plane<Real>::createFrom1proportion2points(aProportion2, aOne2, aOne2 + 1.0f / aProportion2 * (aCommon - aOne2));
  Plane<Real> plane3 = Plane<Real>::createFrom1proportion2points(aProportion3, aOne3, aOne3 + 1.0f / aProportion3 * (aCommon - aOne3));
  return Plane<Real>::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, Proportion) {
  {
    Vertex<Real> common{0.0f, 0.0f, 0.0f};
    auto result = getPlaneIntersectionProportionPoints(common, 0.5f, {1.0f, 0.0f, 0.0f}, 0.5f, {0.0f, 1.0f, 0.0f}, 0.5f, {0.0f, 0.0f, 1.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{0.0f, 0.0f, 0.0f};
    auto result = getPlaneIntersectionProportionPoints(common, 0.5f, {1.0f, 0.0f, 0.0f}, 0.2f, {0.0f, 1.0f, 0.0f}, 0.1f, {0.0f, 0.0f, 1.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{-1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionProportionPoints(common, 0.1f, {10.0f, 10.0f, 0.0f}, 0.2f, {0.0f, 10.0f, 10.0f}, 0.3f, {10.0f, 0.0f, 10.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{-1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionProportionPoints(common, 0.1f, {-10.0f, 10.0f, 0.0f}, 0.2f, {0.0f, -10.0f, 10.0f}, 0.3f, {10.0f, 0.0f, 10.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{-1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionProportionPoints(common, 0.1f, {10.0f, 10.0f, 0.0f}, 0.2f, {0.0f, -10.0f, 10.0f}, 0.3f, {-10.0f, 0.0f, 10.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
}

Vertex<Real> getPlaneIntersectionVertices(Vertex<Real> aCommon,
                                         Vertex<Real> aOther1,
                                         Vertex<Real> aOther2,
                                         Vertex<Real> aOther3) {
  Plane<Real> plane1 = Plane<Real>::createFrom3points(aOther1, aOther2, aCommon);
  Plane<Real> plane2 = Plane<Real>::createFrom3points(aOther2, aOther3, aCommon);
  Plane<Real> plane3 = Plane<Real>::createFrom3points(aOther1, aOther3, aCommon);
  return Plane<Real>::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, Vertices) {
  {
    Vertex<Real> common{1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionVertices(common, {10.0f, 0.0f, 0.0f}, {0.0f, 10.0f, 0.0f}, {0.0f, 0.0f, 10.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionVertices(common, {-10.0f, 0.0f, 0.0f}, {0.0f, 10.0f, 0.0f}, {0.0f, 0.0f, 10.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionVertices(common, {-10.0f, 0.0f, 0.0f}, {0.0f, -10.0f, 0.0f}, {0.0f, 0.0f, 10.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionVertices(common, {-10.0f, 0.0f, 0.0f}, {0.0f, -10.0f, 0.0f}, {0.0f, 0.0f, -10.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
}

Vertex<Real> getPlaneIntersection1vector2points(Vertex<Real> aCommon,
                                                Vertex<Real> aOne1, Vector<Real> aDir1,
                                                Vertex<Real> aOne2, Vector<Real> aDir2,
                                                Vertex<Real> aOne3, Vector<Real> aDir3) {
  Plane<Real> plane1 = Plane<Real>::createFrom1vector2points(aDir1, aOne1, aCommon);
  Plane<Real> plane2 = Plane<Real>::createFrom1vector2points(aDir2, aOne2, aCommon);
  Plane<Real> plane3 = Plane<Real>::createFrom1vector2points(aDir3, aOne3, aCommon);
  return Plane<Real>::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, VectorPoints) {
  {
    Vertex<Real> common{1.0f, 2.0f, -3.0f};
    auto result = getPlaneIntersection1vector2points(common, {10.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f},
                                                             {0.0f, 10.0f, 0.0f}, {0.0f, 0.0f, 1.0f},
                                                             {0.0f, 0.0f, 10.0f}, {1.0f, 0.0f, 0.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{1.0f, 2.0f, -3.0f};
    auto result = getPlaneIntersection1vector2points(common, {10.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 1.0f},
                                                             {0.0f, 10.0f, 0.0f}, {1.0f, 0.0f, -1.0f},
                                                             {0.0f, 0.0f, 10.0f}, {1.0f, 1.0f, 0.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{1.0f, 2.0f, -3.0f};
    auto result = getPlaneIntersection1vector2points(common, {10.0f, 0.0f, 0.0f}, {-4.0f, 1.0f, 1.0f},
                                                             {0.0f, 10.0f, 0.0f}, {1.0f, -4.0f, -1.0f},
                                                             {0.0f, 0.0f, 10.0f}, {1.0f, 1.0f, -4.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
}

Vertex<Real> getPlaneIntersection2vectors1point(Vertex<Real> aCommon,
                                                Vector<Real> aOne1, Vector<Real> aOther1,
                                                Vector<Real> aOne2, Vector<Real> aOther2,
                                                Vector<Real> aOne3, Vector<Real> aOther3) {
  Plane<Real> plane1 = Plane<Real>::createFrom2vectors1point(aOne1, aOther1, aCommon);
  Plane<Real> plane2 = Plane<Real>::createFrom2vectors1point(aOne2, aOther2, aCommon);
  Plane<Real> plane3 = Plane<Real>::createFrom2vectors1point(aOne3, aOther3, aCommon);
  return Plane<Real>::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, VectorsPoint) {
  {
    Vertex<Real> common{1.0f, 2.0f, -3.0f};
    auto result = getPlaneIntersection1vector2points(common, {10.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f},
                                                             {0.0f, 10.0f, 0.0f}, {0.0f, 0.0f, 1.0f},
                                                             {0.0f, 0.0f, 10.0f}, {1.0f, 0.0f, 0.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
  {
    Vertex<Real> common{1.0f, 2.0f, -3.0f};
    auto result = getPlaneIntersection1vector2points(common, {10.0f, 1.0f, 0.0f}, {1.0f, 10.0f, 0.0f},
                                                             {0.0f, 10.0f, 1.0f}, {0.0f, 1.0f, 10.0f},
                                                             {1.0f, 0.0f, 10.0f}, {10.0f, 0.0f, 1.0f});
    EXPECT_LT((common - result).norm(), cgEpsilon);
  }
}

TEST(planeIntersection, Ray) {
  {
    Ray<Real> ray(Vertex<Real>{1.0f, 2.0f, -3.0f}, Vector<Real>{1.0f, 1.0f, 1.0f});
    Plane<Real> plane = Plane<Real>::createFrom3points(Vertex<Real>{10.0f, 1.0f, 2.0f}, Vertex<Real>{11.0f, 11.1f, 2.0f}, Vertex<Real>{12.0f, 1.1f, 4.4f});
    auto result = plane.intersect(ray);
    EXPECT_TRUE(result.mValid);
  }
  {
    Ray<Real> ray(Vertex<Real>{1.0f, 2.0f, -3.0f}, Vector<Real>{-1.0f, 2.0f, 3.0f});
    Plane<Real> plane = Plane<Real>::createFrom3points(Vertex<Real>{10.0f, 1.0f, 2.0f}, Vertex<Real>{11.0f, 11.1f, 2.0f}, Vertex<Real>{12.0f, 1.1f, 4.4f});
    auto result = plane.intersect(ray);
    EXPECT_TRUE(result.mValid);
    EXPECT_LT(result.mDistance, 0.0f);
  }
  {
    Ray<Real> ray(Vertex<Real>{1.0f, 2.0f, -3.0f}, Vector<Real>{0.0f, 2.0f, 0.0f});
    Plane<Real> plane = Plane<Real>::createFrom3points(Vertex<Real>{10.0f, 1.0f, 2.0f}, Vertex<Real>{10.0f, 11.1f, 2.0f}, Vertex<Real>{10.0f, 1.1f, 4.4f});
    auto result = plane.intersect(ray);
    EXPECT_FALSE(result.mValid);
  }
  {
    Ray<Real> ray(Vertex<Real>{1.0f, 2.0f, -3.0f}, Vector<Real>{0.0f, 2.0f, 0.0f});
    Plane<Real> plane = Plane<Real>::createFrom3points(Vertex<Real>{10.0f, 10.0f, 2.0f}, Vertex<Real>{0.0f, 10.0f, 2.0f}, Vertex<Real>{10.0f, 10.0f, 10.4f});
    auto result = plane.intersect(ray);
    EXPECT_TRUE(result.mValid);
    EXPECT_LT((result.mPoint - Vertex<Real>{1.0f, 10.0f, -3.0f}).norm(), 0.00001f);
    EXPECT_GT(::abs(result.mCosIncidence), 0.9999f);
  }
}

TEST(planeProjection, Point) {
  {
    Vertex<Real> point{0.0f, 0.0f, 0.0f};
    Plane<Real> plane = Plane<Real>::createFrom3points({2.0f, 0.0f, 0.0f}, {0.0f, 2.0f, 0.0f}, {0.0f, 0.0f, 2.0f});
    auto projected = plane.project(point);
    EXPECT_LT((projected - Vertex<Real>(0.666666f, 0.666666f, 0.666666f)).norm(), cgEpsilon);
  }
  {
    Vertex<Real> point{0.0f, 0.0f, 0.0f};
    Plane<Real> plane = Plane<Real>::createFrom3points({2.0f, 0.0f, 0.0f}, {2.0f, 1.0f, 0.0f}, {2.0f, 0.0f, 1.0f});
    auto projected = plane.project(point);
    EXPECT_LT((projected - Vertex<Real>(2.0f, 0.0f, 0.0f)).norm(), cgEpsilon);
  }
  {
    Vertex<Real> point{1.0f, 2.0f, 3.0f};
    Plane<Real> plane = Plane<Real>::createFrom3points({3.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 3.0f}, {1.0f, 2.0f, 5.0f});
    auto projected = plane.project(point);
    EXPECT_LT((projected - Vertex<Real>(1.666666f, 2.666666f, 3.666666f)).norm(), cgEpsilon);
  }
  {
    Vertex<Real> point{-1.0f, -2.0f, 3.0f};
    Plane<Real> plane = Plane<Real>::createFrom3points({1.0f, -2.0f, 3.0f}, {1.0f, -3.0f, 3.0f}, {1.0f, -2.0f, 4.0f});
    auto projected = plane.project(point);
    EXPECT_LT((projected - Vertex<Real>(1.0f, -2.0f, 3.0f)).norm(), cgEpsilon);
  }
  {
    Vertex<Real> point{1.666666f, 2.666666f, 3.666666f};
    Plane<Real> plane = Plane<Real>::createFrom3points({3.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 3.0f}, {1.0f, 2.0f, 5.0f});
    auto projected = plane.project(point);
    EXPECT_LT((projected - point).norm(), cgEpsilon);
  }
}

TEST(planeDistance, Point) {
  {
    Vertex<Real> point{0.0f, 0.0f, 0.0f};
    Plane<Real> plane = Plane<Real>::createFrom3points({2.0f, 0.0f, 0.0f}, {0.0f, 2.0f, 0.0f}, {0.0f, 0.0f, 2.0f});
    auto dist = ::abs(plane.distance(point));
    EXPECT_LT(::abs(dist - 1.15468f), cgEpsilon);
  }
  {
    Vertex<Real> point{0.0f, 0.0f, 0.0f};
    Plane<Real> plane = Plane<Real>::createFrom3points({2.0f, 0.0f, 0.0f}, {2.0f, 1.0f, 0.0f}, {2.0f, 0.0f, 1.0f});
    auto dist = ::abs(plane.distance(point));
    EXPECT_LT(::abs(dist - 2.0f), cgEpsilon);
  }
  {
    Vertex<Real> point{1.0f, 2.0f, 3.0f};
    Plane<Real> plane = Plane<Real>::createFrom3points({3.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 3.0f}, {1.0f, 2.0f, 5.0f});
    auto dist = ::abs(plane.distance(point));
    EXPECT_LT(::abs(dist - 1.15468f), cgEpsilon);
  }
  {
    Vertex<Real> point{-1.0f, -2.0f, 3.0f};
    Plane<Real> plane = Plane<Real>::createFrom3points({1.0f, -2.0f, 3.0f}, {1.0f, -3.0f, 3.0f}, {1.0f, -2.0f, 4.0f});
    auto dist = ::abs(plane.distance(point));
    EXPECT_LT(::abs(dist - 2.0f), cgEpsilon);
  }
  {
    Vertex<Real> point{1.666666f, 2.666666f, 3.666666f};
    Plane<Real> plane = Plane<Real>::createFrom3points({3.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 3.0f}, {1.0f, 2.0f, 5.0f});
    auto dist = ::abs(plane.distance(point));
    EXPECT_LT(::abs(dist), cgEpsilon);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
