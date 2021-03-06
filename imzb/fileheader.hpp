#pragma once

#include "imzb/fileutils.hpp"

#include <cstdint>
#include <vector>
#include <fstream>
#include <cassert>

namespace imzb {

struct Mask {
  uint32_t height, width;
  std::vector<uint64_t> data;

  Mask() : height(0), width(0) {}

  Mask(uint32_t height, uint32_t width) { resize(height, width); }

  void resize(uint32_t height, uint32_t width) {
    this->height = height;
    this->width = width;
    data.resize(height * width / 64 + 1);
  }

  bool hasSpectrumAt(uint32_t x, uint32_t y) const {
    assert(x < height && y < width);
    size_t n = x * width + y;
    return bool(data[n / 64] & (1ULL << (n % 64)));
  }

  void set(uint32_t x, uint32_t y) {
    assert(x < height && y < width);
    size_t n = x * width + y;
    data[n / 64] |= (1ULL << (n % 64));
  }
};

struct FileHeader {
  uint32_t height() const { return mask.height; }
  uint32_t width() const { return mask.width; }
  Mask mask;

  uint32_t version;
  uint32_t block_size;

  double min_mz = 0;
  double max_mz = 0;

  void read(std::ifstream& stream) {
    uint64_t pos = stream.tellg();

    uint32_t header_len;
    binary_read(stream, header_len);
    binary_read(stream, mask.height);
    binary_read(stream, mask.width);
    mask.resize(mask.height, mask.width);
    stream.read(
        reinterpret_cast<char*>(&mask.data[0]), mask.data.size() * sizeof(uint64_t));

    if (stream.tellg() == pos + header_len) {
      version = 0;
      block_size = 4096;
    } else {
      binary_read(stream, version);
      if (version >= 1) binary_read(stream, block_size);
      if (version >= 2) {
        binary_read(stream, min_mz);
        binary_read(stream, max_mz);
      }
    }

    stream.seekg(pos + header_len, std::ios::beg);
  }

  void write(std::ofstream& stream) {
    uint32_t header_len = 3 * 4 + mask.data.size() * sizeof(uint64_t);
    header_len += sizeof(version) + sizeof(block_size);
    header_len += sizeof(min_mz) + sizeof(max_mz);

    binary_write(stream, header_len);
    binary_write(stream, mask.height);
    binary_write(stream, mask.width);
    stream.write(
        reinterpret_cast<char*>(&mask.data[0]), mask.data.size() * sizeof(uint64_t));

    binary_write(stream, version);
    binary_write(stream, block_size);
    binary_write(stream, min_mz);
    binary_write(stream, max_mz);
  }
};
}
