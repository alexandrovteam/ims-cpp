#pragma once

#include "imzb/fileheader.hpp"
#include "imzb/fileutils.hpp"

#include <vector>
#include <memory>
#include <cstdint>
#include <algorithm>
#include <fstream>

namespace imzb {

struct Index {
  imzb::FileHeader header;

  std::vector<double> mzs;
  std::vector<uint64_t> offsets;

  size_t startBlock(double mz) const {
    size_t idx = std::lower_bound(mzs.begin(), mzs.end(), mz) - mzs.begin();
    return idx > 0 ? (idx - 1) : idx;
  }

  uint64_t startOffset(double mz) const {
    return offsets[startBlock(mz)];
  }

  size_t endBlock(double mz) const {
    return std::upper_bound(mzs.begin(), mzs.end(), mz) - mzs.begin();
  }

  uint64_t endOffset(double mz) const {
    return offsets[endBlock(mz)];
  }

  void read(std::ifstream& stream) {
    header.read(stream);
    double mz;
    uint64_t offset;
    while (binary_read(stream, mz)) {
      assert(binary_read(stream, offset));
      mzs.push_back(mz);
      offsets.push_back(offset);
    }
  }

  void write(std::ofstream& stream) {
    header.write(stream);
    for (size_t i = 0; i < mzs.size(); ++i) {
      binary_write(stream, mzs[i]);
      binary_write(stream, offsets[i]);
    }
  }
};

typedef std::shared_ptr<Index> IndexPtr;

}
