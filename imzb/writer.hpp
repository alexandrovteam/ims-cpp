#pragma once

#include "ims/ims.hpp"
#include "imzb/fileheader.hpp"
#include "imzb/index.hpp"

#include <string>
#include <fstream>
#include <vector>

namespace imzb {

class ImzbWriter {
public:
  static const uint32_t BLOCK_SIZE = 4096;

private:
  std::ofstream out_;
  imzb::Index index_;

  std::vector<ims::Peak> in_buf_;
  std::vector<char> out_buf_;
  uint32_t filled_;
  uint64_t bytes_written_;

  void dump();

  std::string filename_;

public:
  ImzbWriter(const std::string& filename);

  void setMask(const imzb::Mask& mask) {
    index_.header.mask = mask;
  }

  void writePeak(const ims::Peak& peak);

  void close();

  ~ImzbWriter();
};

}
