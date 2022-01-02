#ifndef CUDA_BEZIER_TRIANGLE_RAYTRACER_HOSTUTIL
#define CUDA_BEZIER_TRIANGLE_RAYTRACER_HOSTUTIL

#include "util.h"
#include <vector>
#include <random>

class UniformHemisphere final {
private:
  float const                             cmBeltWidth;
  std::vector<std::pair<float, uint32_t>>  mBelts;
  std::ranlux24_base                       mRandomGenerator;
  std::uniform_real_distribution<float>    mUniform1;
  std::uniform_real_distribution<float>    mUniformPi_2;
  std::uniform_real_distribution<float>    mUniform2pi;
  uint32_t                                 mPatchCount;

public:
  UniformHemisphere(uint32_t const aBelts);

  uint32_t getPatchCount() const { return mPatchCount; }

  // Returns unit vector and the patch index.
  std::pair<Vector, uint32_t> getRandom();
};

#endif
