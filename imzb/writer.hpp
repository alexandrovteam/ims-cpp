#pragma once

#include "ims/ims.hpp"
#include "imzb/fileheader.hpp"
#include "imzb/index.hpp"

#include <string>
#include <fstream>
#include <vector>

namespace imzb {

struct ImzbCompressionSettings {
  std::string compressor = "blosclz";
  uint8_t compression_level = 5;

  uint32_t block_size = 4096;
};

class ImzbWriter {
 private:
  std::ofstream out_;
  imzb::Index index_;

  std::vector<ims::Peak> in_buf_;
  std::vector<char> out_buf_;
  uint32_t filled_;
  uint64_t bytes_written_;

  void dump();

  std::string filename_;

  ImzbCompressionSettings c_;

  double min_mz_, max_mz_;

 public:
  ImzbWriter(const std::string& filename,
      const ImzbCompressionSettings& compression_settings = ImzbCompressionSettings());

  void setMask(const imzb::Mask& mask) { index_.header.mask = mask; }

  void writePeak(const ims::Peak& peak);

  void close();

  ~ImzbWriter();
};
}
