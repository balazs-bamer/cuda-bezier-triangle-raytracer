#include "util.h"

#include<deque>
#include<vector>
#include<iostream>

#include "gtest/gtest.h"


constexpr float cgEpsilon = 0.0001f;

std::ostream& operator<<(std::ostream &aOut, Plane const &aPlane) {
  aOut << aPlane.mNormal(0) << ' ' << aPlane.mNormal(1) << ' ' << aPlane.mNormal(2) << ' ' << " : " << aPlane.mConstant;
  return aOut;
}

bool eq(Vector const &aV1, Vector const &aV2, float const aEpsilon) {
  return (aV1 - aV2).norm() < aEpsilon;
}

bool eq(Vector const &aV1, Vector const &aV2) {
  return eq(aV1, aV2, cgEpsilon);
}

bool eq(float const aF1, float const aF2, float const aEpsilon) {
  return ::abs(aF1 - aF2) < aEpsilon;
}

bool eq(float const aF1, float const aF2) {
  return eq(aF1, aF2, cgEpsilon);
}

Vertex getPlaneIntersectionNormals(Vertex aPoint,
                                         Vector aDir1,
                                         Vector aDir2,
                                         Vector aDir3) {
  Vector normal1 = aDir1.normalized();
  Plane plane1(normal1, normal1.dot(aPoint));
  Vector normal2 = aDir2.normalized();
  Plane plane2(normal2, normal2.dot(aPoint));
  Vector normal3 = aDir3.normalized();
  Plane plane3(normal3, normal3.dot(aPoint));
  return Plane::intersect(plane1, plane2, plane3);
}

TEST(vector, getAperpendicular) {
  static constexpr float csEpsilon = 1e-10;
  {
    Vector vector{1.0f, 0.0f, 0.0f};
    EXPECT_LT(::abs(util::getAperpendicular(vector).dot(vector)), csEpsilon);
  }
  {
    Vector vector{1.0f, 1.0f, 0.0f};
    vector.normalize();
    EXPECT_LT(::abs(util::getAperpendicular(vector).dot(vector)), csEpsilon);
  }
  {
    Vector vector{1.0f, 0.0f, 1.0f};
    vector.normalize();
    EXPECT_LT(::abs(util::getAperpendicular(vector).dot(vector)), csEpsilon);
  }
  {
    Vector vector{1.0f, -1.0f, -1.0f};
    vector.normalize();
    EXPECT_LT(::abs(util::getAperpendicular(vector).dot(vector)), csEpsilon);
  }
}

TEST(ray, averageErrorSquared) {
  {
    Ray ray({0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f});
    std::vector<Vertex> points{};
    EXPECT_EQ(ray.getAverageErrorSquared(points), 0.0f);
  }
  {
    Ray ray({0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f});
    std::vector<Vertex> points{{2.0f, 0.0f, 0.0f}, {-3.0f, 0.0f, 0.0f}};
    EXPECT_EQ(ray.getAverageErrorSquared(points), 0.0f);
  }
  {
    Ray ray({0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f});
    std::vector<Vertex> points{{2.0f, 1.0f, 0.0f}, {-3.0f, 0.0f, 1.0f}};
    EXPECT_GT(ray.getAverageErrorSquared(points), 0.0f);
  }
}

TEST(planeIntersection, Normals) {
  {
    Vertex common{1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionNormals(common, {1.0, 2.0f, 3.0f}, {3.0f, 1.0f, 2.0f}, {3.0f, 2.0f, 1.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{3.0f, -2.0f, 1.0f};
    auto result = getPlaneIntersectionNormals(common, {1.0, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{3.0f, -2.0f, -1.0f};
    auto result = getPlaneIntersectionNormals(common, {1.0, -2.0f, 3.0f}, {-1.0f, 2.0f, 3.0f}, {1.0f, 2.0f, -3.0f});
    EXPECT_TRUE(eq(common, result));
  }
}

Vertex getPlaneIntersectionProportionPoints(Vertex aCommon,
                                                  float aProportion1, Vertex aOne1,
                                                  float aProportion2, Vertex aOne2,
                                                  float aProportion3, Vertex aOne3) {
  Plane plane1 = Plane::createFrom1proportion2points(aProportion1, aOne1, aOne1 + 1.0f / aProportion1 * (aCommon - aOne1));
  Plane plane2 = Plane::createFrom1proportion2points(aProportion2, aOne2, aOne2 + 1.0f / aProportion2 * (aCommon - aOne2));
  Plane plane3 = Plane::createFrom1proportion2points(aProportion3, aOne3, aOne3 + 1.0f / aProportion3 * (aCommon - aOne3));
  return Plane::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, Proportion) {
  {
    Vertex common{0.0f, 0.0f, 0.0f};
    auto result = getPlaneIntersectionProportionPoints(common, 0.5f, {1.0f, 0.0f, 0.0f}, 0.5f, {0.0f, 1.0f, 0.0f}, 0.5f, {0.0f, 0.0f, 1.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{0.0f, 0.0f, 0.0f};
    auto result = getPlaneIntersectionProportionPoints(common, 0.5f, {1.0f, 0.0f, 0.0f}, 0.2f, {0.0f, 1.0f, 0.0f}, 0.1f, {0.0f, 0.0f, 1.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{-1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionProportionPoints(common, 0.1f, {10.0f, 10.0f, 0.0f}, 0.2f, {0.0f, 10.0f, 10.0f}, 0.3f, {10.0f, 0.0f, 10.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{-1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionProportionPoints(common, 0.1f, {-10.0f, 10.0f, 0.0f}, 0.2f, {0.0f, -10.0f, 10.0f}, 0.3f, {10.0f, 0.0f, 10.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{-1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionProportionPoints(common, 0.1f, {10.0f, 10.0f, 0.0f}, 0.2f, {0.0f, -10.0f, 10.0f}, 0.3f, {-10.0f, 0.0f, 10.0f});
    EXPECT_TRUE(eq(common, result));
  }
}

Vertex getPlaneIntersectionVertices(Vertex aCommon,
                                         Vertex aOther1,
                                         Vertex aOther2,
                                         Vertex aOther3) {
  Plane plane1 = Plane::createFrom3points(aOther1, aOther2, aCommon);
  Plane plane2 = Plane::createFrom3points(aOther2, aOther3, aCommon);
  Plane plane3 = Plane::createFrom3points(aOther1, aOther3, aCommon);
  return Plane::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, Vertices) {
  {
    Vertex common{1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionVertices(common, {10.0f, 0.0f, 0.0f}, {0.0f, 10.0f, 0.0f}, {0.0f, 0.0f, 10.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionVertices(common, {-10.0f, 0.0f, 0.0f}, {0.0f, 10.0f, 0.0f}, {0.0f, 0.0f, 10.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionVertices(common, {-10.0f, 0.0f, 0.0f}, {0.0f, -10.0f, 0.0f}, {0.0f, 0.0f, 10.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0f, 2.0f, 3.0f};
    auto result = getPlaneIntersectionVertices(common, {-10.0f, 0.0f, 0.0f}, {0.0f, -10.0f, 0.0f}, {0.0f, 0.0f, -10.0f});
    EXPECT_TRUE(eq(common, result));
  }
}

Vertex getPlaneIntersection1vector2points(Vertex aCommon,
                                                Vertex aOne1, Vector aDir1,
                                                Vertex aOne2, Vector aDir2,
                                                Vertex aOne3, Vector aDir3) {
  Plane plane1 = Plane::createFrom1vector2points(aDir1, aOne1, aCommon);
  Plane plane2 = Plane::createFrom1vector2points(aDir2, aOne2, aCommon);
  Plane plane3 = Plane::createFrom1vector2points(aDir3, aOne3, aCommon);
  return Plane::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, VectorPoints) {
  {
    Vertex common{1.0f, 2.0f, -3.0f};
    auto result = getPlaneIntersection1vector2points(common, {10.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f},
                                                             {0.0f, 10.0f, 0.0f}, {0.0f, 0.0f, 1.0f},
                                                             {0.0f, 0.0f, 10.0f}, {1.0f, 0.0f, 0.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0f, 2.0f, -3.0f};
    auto result = getPlaneIntersection1vector2points(common, {10.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 1.0f},
                                                             {0.0f, 10.0f, 0.0f}, {1.0f, 0.0f, -1.0f},
                                                             {0.0f, 0.0f, 10.0f}, {1.0f, 1.0f, 0.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0f, 2.0f, -3.0f};
    auto result = getPlaneIntersection1vector2points(common, {10.0f, 0.0f, 0.0f}, {-4.0f, 1.0f, 1.0f},
                                                             {0.0f, 10.0f, 0.0f}, {1.0f, -4.0f, -1.0f},
                                                             {0.0f, 0.0f, 10.0f}, {1.0f, 1.0f, -4.0f});
    EXPECT_TRUE(eq(common, result));
  }
}

Vertex getPlaneIntersection2vectors1point(Vertex aCommon,
                                                Vector aOne1, Vector aOther1,
                                                Vector aOne2, Vector aOther2,
                                                Vector aOne3, Vector aOther3) {
  Plane plane1 = Plane::createFrom2vectors1point(aOne1, aOther1, aCommon);
  Plane plane2 = Plane::createFrom2vectors1point(aOne2, aOther2, aCommon);
  Plane plane3 = Plane::createFrom2vectors1point(aOne3, aOther3, aCommon);
  return Plane::intersect(plane1, plane2, plane3);
}

TEST(planeIntersection, VectorsPoint) {
  {
    Vertex common{1.0f, 2.0f, -3.0f};
    auto result = getPlaneIntersection1vector2points(common, {10.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f},
                                                             {0.0f, 10.0f, 0.0f}, {0.0f, 0.0f, 1.0f},
                                                             {0.0f, 0.0f, 10.0f}, {1.0f, 0.0f, 0.0f});
    EXPECT_TRUE(eq(common, result));
  }
  {
    Vertex common{1.0f, 2.0f, -3.0f};
    auto result = getPlaneIntersection1vector2points(common, {10.0f, 1.0f, 0.0f}, {1.0f, 10.0f, 0.0f},
                                                             {0.0f, 10.0f, 1.0f}, {0.0f, 1.0f, 10.0f},
                                                             {1.0f, 0.0f, 10.0f}, {10.0f, 0.0f, 1.0f});
    EXPECT_TRUE(eq(common, result));
  }
}

TEST(planeIntersection, Ray) {
  {
    Ray ray(Vertex{1.0f, 2.0f, -3.0f}, Vector{1.0f, 1.0f, 1.0f});
    Plane plane = Plane::createFrom3points(Vertex{10.0f, 1.0f, 2.0f}, Vertex{11.0f, 11.1f, 2.0f}, Vertex{12.0f, 1.1f, 4.4f});
    auto result = plane.intersect(ray);
    EXPECT_TRUE(result.mValid);
  }
  {
    Ray ray(Vertex{1.0f, 2.0f, -3.0f}, Vector{-1.0f, 2.0f, 3.0f});
    Plane plane = Plane::createFrom3points(Vertex{10.0f, 1.0f, 2.0f}, Vertex{11.0f, 11.1f, 2.0f}, Vertex{12.0f, 1.1f, 4.4f});
    auto result = plane.intersect(ray);
    EXPECT_TRUE(result.mValid);
    EXPECT_LT(result.mDistance, 0.0f);
  }
  {
    Ray ray(Vertex{1.0f, 2.0f, -3.0f}, Vector{0.0f, 2.0f, 0.0f});
    Plane plane = Plane::createFrom3points(Vertex{10.0f, 1.0f, 2.0f}, Vertex{10.0f, 11.1f, 2.0f}, Vertex{10.0f, 1.1f, 4.4f});
    auto result = plane.intersect(ray);
    EXPECT_FALSE(result.mValid);
  }
  {
    Ray ray(Vertex{1.0f, 2.0f, -3.0f}, Vector{0.0f, 2.0f, 0.0f});
    Plane plane = Plane::createFrom3points(Vertex{10.0f, 10.0f, 2.0f}, Vertex{0.0f, 10.0f, 2.0f}, Vertex{10.0f, 10.0f, 10.4f});
    auto result = plane.intersect(ray);
    EXPECT_TRUE(result.mValid);
    EXPECT_LT((result.mPoint - Vertex{1.0f, 10.0f, -3.0f}).norm(), 0.00001f);
    EXPECT_GT(::abs(result.mCosIncidence), 0.9999f);
  }
}

TEST(planeProjection, Point) {
  {
    Vertex point{0.0f, 0.0f, 0.0f};
    Plane plane = Plane::createFrom3points({2.0f, 0.0f, 0.0f}, {0.0f, 2.0f, 0.0f}, {0.0f, 0.0f, 2.0f});
    auto projected = plane.project(point);
    EXPECT_TRUE(eq(projected, Vertex(0.666666f, 0.666666f, 0.666666f)));
  }
  {
    Vertex point{0.0f, 0.0f, 0.0f};
    Plane plane = Plane::createFrom3points({2.0f, 0.0f, 0.0f}, {2.0f, 1.0f, 0.0f}, {2.0f, 0.0f, 1.0f});
    auto projected = plane.project(point);
    EXPECT_TRUE(eq(projected, Vertex(2.0f, 0.0f, 0.0f)));
  }
  {
    Vertex point{1.0f, 2.0f, 3.0f};
    Plane plane = Plane::createFrom3points({3.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 3.0f}, {1.0f, 2.0f, 5.0f});
    auto projected = plane.project(point);
    EXPECT_TRUE(eq(projected, Vertex(1.666666f, 2.666666f, 3.666666f)));
  }
  {
    Vertex point{-1.0f, -2.0f, 3.0f};
    Plane plane = Plane::createFrom3points({1.0f, -2.0f, 3.0f}, {1.0f, -3.0f, 3.0f}, {1.0f, -2.0f, 4.0f});
    auto projected = plane.project(point);
    EXPECT_TRUE(eq(projected, Vertex(1.0f, -2.0f, 3.0f)));
  }
  {
    Vertex point{1.666666f, 2.666666f, 3.666666f};
    Plane plane = Plane::createFrom3points({3.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 3.0f}, {1.0f, 2.0f, 5.0f});
    auto projected = plane.project(point);
    EXPECT_TRUE(eq(projected, point));
  }
}

TEST(planeDistance, Point) {
  {
    Vertex point{0.0f, 0.0f, 0.0f};
    Plane plane = Plane::createFrom3points({2.0f, 0.0f, 0.0f}, {0.0f, 2.0f, 0.0f}, {0.0f, 0.0f, 2.0f});
    auto dist = ::abs(plane.distance(point));
    EXPECT_TRUE(eq(dist, 1.15468f));
  }
  {
    Vertex point{0.0f, 0.0f, 0.0f};
    Plane plane = Plane::createFrom3points({2.0f, 0.0f, 0.0f}, {2.0f, 1.0f, 0.0f}, {2.0f, 0.0f, 1.0f});
    auto dist = ::abs(plane.distance(point));
    EXPECT_TRUE(eq(dist, 2.0f));
  }
  {
    Vertex point{1.0f, 2.0f, 3.0f};
    Plane plane = Plane::createFrom3points({3.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 3.0f}, {1.0f, 2.0f, 5.0f});
    auto dist = ::abs(plane.distance(point));
    EXPECT_TRUE(eq(dist, 1.15468f));
  }
  {
    Vertex point{-1.0f, -2.0f, 3.0f};
    Plane plane = Plane::createFrom3points({1.0f, -2.0f, 3.0f}, {1.0f, -3.0f, 3.0f}, {1.0f, -2.0f, 4.0f});
    auto dist = ::abs(plane.distance(point));
    EXPECT_TRUE(eq(dist, 2.0f));
  }
  {
    Vertex point{1.666666f, 2.666666f, 3.666666f};
    Plane plane = Plane::createFrom3points({3.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 3.0f}, {1.0f, 2.0f, 5.0f});
    auto dist = ::abs(plane.distance(point));
    EXPECT_TRUE(eq(dist, 0.0f));
  }
}

TEST(toWhichSide, Points) {
  Vertex triangle0{3.0f, 2.0f, 5.0f};
  Vertex triangle1{1.0f, 4.0f, 5.0f};
  Vertex triangle2{6.0f, 5.0f, 5.0f};
  Vertex start = (triangle0 + triangle1 + triangle2) / 3.0f;
  {
    Vertex end = start + Vector{1.0f, 0.0f, 0.0f};
    auto converter = util::getBarycentricInverse(triangle0, triangle1, triangle2);
    EXPECT_EQ(util::toWhichSide(converter * start, converter * end), 2u);
  }
  {
    Vertex end = start + Vector{0.0f, 1.0f, 0.0f};
    auto converter = util::getBarycentricInverse(triangle0, triangle1, triangle2);
    EXPECT_EQ(util::toWhichSide(converter * start, converter * end), 1u);
  }
  {
    Vertex end = start + Vector{-1.0f, -1.0f, 0.0f};
    auto converter = util::getBarycentricInverse(triangle0, triangle1, triangle2);
    EXPECT_EQ(util::toWhichSide(converter * start, converter * end), 0u);
  }
}

TEST(polynomApprox, x2few) {
  PolynomApprox<2u> poly({1.0f, 2.0f}, {1.0f, 4.0f, 9.0f});
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(1.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(4.0f),  0.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), std::numeric_limits<float>::max()));
}

TEST(polynomApprox, x20) {
  PolynomApprox<2u> poly({1.0f, 2.0f, 3.0f}, {0.0f, 0.0f, 0.0f});
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(4.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(5.0f),  0.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2_21_27_30) {
  PolynomApprox<2u> poly({32.0f, 48.0f, 63.0f}, {21.0f, 27.0f, 30.0f});
  EXPECT_TRUE(eq(poly.eval(32.0f),  21.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(48.0f),  27.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(63.0f),  30.0f, 0.1f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.3f));
}

TEST(polynomApprox, x2) {
  PolynomApprox<2u> poly({1.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 9.0f});
  EXPECT_TRUE(poly.size() == 3u);
  EXPECT_TRUE(eq(poly[0u],  0.0f));
  EXPECT_TRUE(eq(poly[1u],  0.0f));
  EXPECT_TRUE(eq(poly[2u],  1.0f));
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(4.0f), 16.0f));
  EXPECT_TRUE(eq(poly.eval(5.0f), 25.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2_5) {
  PolynomApprox<2u> poly({1.0f, 2.0f, 3.0f}, {6.0f, 9.0f, 14.0f});
  EXPECT_TRUE(poly.size() == 3u);
  EXPECT_TRUE(eq(poly[0u],  5.0f));
  EXPECT_TRUE(eq(poly[1u],  0.0f));
  EXPECT_TRUE(eq(poly[2u],  1.0f));
  EXPECT_TRUE(eq(poly.eval(0.0f),  5.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(4.0f), 21.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(5.0f), 30.0f, 0.1f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2off) {
  PolynomApprox<2u> poly({1.0f, 2.0f, 3.0f, 4.0f}, {1.0f, 4.0f, 9.0f, 15.9f});
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(4.0f), 16.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(5.0f), 25.0f, 0.3f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.00059f));
}

TEST(polynomApprox, x2veryoff) {
  PolynomApprox<2u> poly({1.0f, 2.0f, 3.0f, 4.0f}, {1.0f, 3.0f, 13.0f, 13.0f});
  EXPECT_FALSE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_FALSE(eq(poly.eval(4.0f), 16.0f));
  EXPECT_FALSE(eq(poly.eval(5.0f), 25.0f));
  EXPECT_TRUE(poly.getRrmsError() > 0.1f);
}

TEST(polynomApprox, 2x3_3x2_4x_5) {
  PolynomApprox<3u> poly({0.0f, 1.0f, 2.0f, 3.0f}, {5.0f, 14.0f, 41.0f, 98.0f});
  EXPECT_TRUE(poly.size() == 4u);
  EXPECT_TRUE(eq(poly.eval(0.5f),  8.0f));
  EXPECT_TRUE(eq(poly.eval(1.5f), 24.5f));
  EXPECT_TRUE(eq(poly.eval(2.5f), 65.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
