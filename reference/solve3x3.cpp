#include <Eigen/Dense>
#include <chrono>
#include <iostream>

/*
mult           0.0163621 s
solveInverse   0.0201562 s
solveFullPivLu 0.202951 s
solveColPivHh  0.547893 s
solveFullPivHh 0.621699 s
solveComplete  0.844406 s

*/

Eigen::Vector3f mult(Eigen::Matrix3f const &m, Eigen::Vector3f const &v) {
  return m * v;
}
 
Eigen::Vector3f solveInverse(Eigen::Matrix3f const &m, Eigen::Vector3f const &v) {
  return m.inverse() * v;
}
 
Eigen::Vector3f solveFullPivLu(Eigen::Matrix3f const &m, Eigen::Vector3f const &v) {
  return m.fullPivLu().solve(v);
} 

Eigen::Vector3f solveColPivHh(Eigen::Matrix3f const &m, Eigen::Vector3f const &v) {
  return m.colPivHouseholderQr().solve(v);
} 

Eigen::Vector3f solveFullPivHh(Eigen::Matrix3f const &m, Eigen::Vector3f const &v) {
  return m.fullPivHouseholderQr().solve(v);
} 

Eigen::Vector3f solveComplete(Eigen::Matrix3f const &m, Eigen::Vector3f const &v) {
  return m.completeOrthogonalDecomposition().solve(v);
} 

Eigen::Vector3f solveBdcsvd(Eigen::Matrix3f const &m, Eigen::Vector3f const &v) {
  return m.bdcSvd().solve(v);
} 

Eigen::Vector3f solveJacobiSvd(Eigen::Matrix3f const &m, Eigen::Vector3f const &v) {
  return m.jacobiSvd().solve(v);
} 

typedef Eigen::Vector3f func(Eigen::Matrix3f const &m, Eigen::Vector3f const &v);

int measure(char const *t, func &f) {
  Eigen::Matrix3f m;
  Eigen::Vector3f v,s;
  s << 0.0f, 0.0f, 0.0f;
  int n = 0;
  auto start = std::chrono::high_resolution_clock::now();
  for(int i = 0; i < 1000000; ++i) {
    m << ++n, ++n, ++n, ++n, ++n, ++n, ++n, ++n, ++n;
    v << ++n, ++n, ++n;
    n %= 71;
    s += f(m, v);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << t << ' ' << diff.count() << " s\n";
  return s(0);
}

int main() {
  return measure("mult", mult) +
  measure("solveInverse", solveInverse) +
  measure("solveFullPivLu", solveFullPivLu) +
  measure("solveColPivHh", solveColPivHh) +
  measure("solveFullPivHh", solveFullPivHh) +
  measure("solveComplete", solveComplete)/* +
  measure("solveBdcSvd", solveBdcsvd) +
  measure("solveJacobiSvd", solveJacobiSvd)*/;
}
