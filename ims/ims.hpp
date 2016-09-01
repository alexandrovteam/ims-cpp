#pragma once

#include "ims/common.hpp"
#include "ims/image.hpp"
#include "ims/image_measures.hpp"

#include <memory>

namespace ims {

class AbstractReader {
 public:
  virtual bool readNextSpectrum(ims::Spectrum&) = 0;
  virtual uint32_t height() const = 0;
  virtual uint32_t width() const = 0;
  virtual ~AbstractReader() {}
};

typedef std::shared_ptr<AbstractReader> AbstractReaderPtr;
}
