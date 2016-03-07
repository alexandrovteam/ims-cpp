#pragma once

#include "ims/image.hpp"
#include "ims/image_measures.hpp"

#include <vector>
#include <cstdint>
#include <memory>

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

class AbstractReader {
public:
  virtual bool readNextSpectrum(ims::Spectrum&) = 0;
  virtual uint32_t height() const = 0;
  virtual uint32_t width() const = 0;
  virtual ~AbstractReader() {}
};

typedef std::shared_ptr<AbstractReader> AbstractReaderPtr;

}
