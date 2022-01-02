#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_HOSTUTIL
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_HOSTUTIL

#include "util.h"
#include <vector>

class UniformHemisphere final {
private:
  float                                   const cmBeltWidth;
  std::vector<std::pair<float, uint32_t>> const mBelts;
  std::ranlux24_base                            mRandomGenerator;
  std::uniform_real_distribution<float>         mUniform1;
  std::uniform_real_distribution<float>         mUniformPi_2;
  std::uniform_real_distribution<float>         mUniform2pi;

  UniformHemisphere(uint32_t const aBelts);

  // Returns unit vector and the patch index.
  std::pair<Vector, uint32_t> getRandom() const;
};

#endif
