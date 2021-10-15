#include <Eigen/Dense>
#include <chrono>
#include <iostream>

Eigen::Vector3f mult(Eigen::Matrix3f const &m, Eigen::Vector3f const &v) {
  return m*v;
}
 
Eigen::Vector3f solve(Eigen::Matrix3f const &m, Eigen::Vector3f const &v) {
  return m.colPivHouseholderQr().solve(v);
} 

typedef Eigen::Vector3f func(Eigen::Matrix3f const &m, Eigen::Vector3f const &v);

int measure(char const *t, func &f) {
  Eigen::Matrix3f m;
  Eigen::Vector3f v,s;
  s << 0.0f, 0.0f, 0.0f;
  int n = 0;
  auto start = std::chrono::high_resolution_clock::now();
  for(int i = 0; i < 10000; ++i) {
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
  measure("solve", solve);
}
