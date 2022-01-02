#include "hostUtil.h"
#include <random>

UniformHemisphere::UniformHemisphere(int32_t const aBelts)
  : cmBeltWidth(cgPi / 2.0f / aBelts)
  , mUniform1(0,0f, 1.0f)
  , mUniformPi_2(0,0f, cgPi / 2.0f)
  , mUniform2pi(0,0f, cgPi * 2.0f) {
  mBelts.reserve(aBelts);
  uint32_t soFar = 0u;
  for(uint32_t i = 0u; i < aBelts; ++i) {
    uint32_t now = ::ceil(4.0f * aBelts * ::sin((2.0f * i + 1.0f) / (4.0f * mBelts) * cgPi));
    mBelts.emplace_back(cgPi * 2.0f / now, soFar);
    soFar += now;
  }
}

std::pair<Vector, uint32_t> UniformHemisphere::getRandom() const {
  Vector direction;
  uint32_t index;
  float inclination;
  float beltRadius;
  do {
    inclination = mUniformPi_2(mRandomGenerator);
    beltRadius = ::sin(inclination);
  } while(beltRadius < mUniform1(mRandomNumberGenerator));
  auto turn = mUniform2pi(mRandomNumberGenerator);
  direction(0) = ::cos(inclination);
  direction(1) = beltRadius * ::cos(turn);
  direction(2) = beltRadius * ::sin(turn);
  direction.normalize();
  auto [patchWidth, soFar] = mBelts[static_cast<uint32_t>(inclination / cmBeltWidth)];
  index = soFar + static_cast<uint32_t>(turn / patchWidth);
  return std::pair(direction, index);
}
