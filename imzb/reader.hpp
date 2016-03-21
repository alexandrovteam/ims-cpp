#pragma once

#include "ims/ims.hpp"
#include "ims/image.hpp"
#include "imzb/index.hpp"

#include <string>
#include <fstream>
#include <vector>
#include <cstdint>

namespace imzb {

class ImzbReader {
  std::string fn_;
  std::ifstream in_;

  IndexPtr index_;
  size_t block_idx_;

  std::vector<char> buffer_;
  std::vector<ims::Peak> peaks_;
  size_t n_peaks_;
  size_t pos_;

  size_t decompressBlock(size_t block_idx,
                         std::ifstream& in,
                         std::vector<char>& inbuf,
                         std::vector<ims::Peak>& outbuf) const;

  bool readNextBlock();

  bool empty_;

public:
  ImzbReader(const std::string& filename);

  bool readNext(ims::Peak& peak);

  void reset();

  std::vector<ims::Peak> slice(double min_mz, double max_mz) const;

  ims::Image<float> image(double mz, double ppm) const;

  // for centroided data
  void readImage(double mz, double ppm, float* image) const;

  // for raw data
  void readCentroidedImage(double mz, double ppm, float* image) const;

  uint32_t height() const { return index_->header.mask.height; }
  uint32_t width() const { return index_->header.mask.width; }

  void close() { in_.close(); }
  const std::string& filename() const { return fn_; }
};

}
