#include "mesh.h"
#include "bezierMesh.h"

#include<deque>
#include<vector>
#include<iostream>

#include "gtest/gtest.h"

using Real = float;

constexpr Real cgEpsilon = 0.001f;

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
/*  Matrix<Real> matrix {
    { plane1.mNormal(0), plane1.mNormal(1), plane1.mNormal(2) },
    { plane2.mNormal(0), plane2.mNormal(1), plane2.mNormal(2) },
    { plane3.mNormal(0), plane3.mNormal(1), plane3.mNormal(2) }
  };
  std::cout << "p1: " << plane1 << "\np2: " << plane2 << "\np3: " << plane3 << '\n';
  std::cout << "m:\n" << matrix << '\n';
  Vector<Real> vector{ plane1.mConstant, plane2.mConstant, plane3.mConstant };
  std::cout << "v:\n" << vector << '\n';
  Matrix<Real> inverse = matrix.inverse();
  std::cout << "i:\n" << inverse << '\n';
  Vertex<Real> result = Plane<Real>::intersect(plane1, plane2, plane3);
  std::cout << "det:\n" << matrix.determinant() << "\n";
  std::cout << "res:\n" << result << "\n";
  std::cout << "error: " << ((matrix * result - vector).norm() / vector.norm()) << "\n\n";*/
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

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
