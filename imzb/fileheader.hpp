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

  Mask(): height(0), width(0) {}

  Mask(uint32_t height, uint32_t width) {
    resize(height, width);
  }

  void resize(uint32_t height, uint32_t width) {
    this->height = height;
    this->width = width;
    data.resize(height * width / 64 + 1);
  }

  bool hasSpectrumAt(uint32_t x, uint32_t y) {
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

  void read(std::ifstream& stream) {
    uint64_t pos = stream.tellg();

    uint32_t header_len;
    binary_read(stream, header_len);
    binary_read(stream, mask.height);
    binary_read(stream, mask.width);
    mask.resize(mask.height, mask.width);
    stream.read(reinterpret_cast<char*>(&mask.data[0]),
                mask.data.size() * sizeof(uint64_t));

    stream.seekg(pos + header_len, std::ios::beg);
  }

  void write(std::ofstream& stream) {
    uint32_t header_len = 3 * 4 + mask.data.size() * sizeof(uint64_t);
    binary_write(stream, header_len);
    binary_write(stream, mask.height);
    binary_write(stream, mask.width);
    stream.write(reinterpret_cast<char*>(&mask.data[0]),
                 mask.data.size() * sizeof(uint64_t));
  }
};

}
