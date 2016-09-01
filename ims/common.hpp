#pragma once

#include <cstdint>
#include <vector>

namespace ims {

struct Position {
  uint32_t x;
  uint32_t y;
  uint32_t z;
};

struct Spectrum {
  std::vector<double> mzs;
  std::vector<float> intensities;

  Position coords;
};

struct Peak {
  Position coords;
  double mz;
  float intensity;
};

}  // namespace ims
