#include "imzb/writer.hpp"
#include "imzb/version.hpp"

#include "blosc.h"

#include <ios>
#include <cassert>
#include <stdexcept>
#include <limits>

imzb::ImzbWriter::ImzbWriter(
    const std::string& filename, const ImzbCompressionSettings& settings)
    : out_{filename, std::ios::binary},
      in_buf_(settings.block_size),
      out_buf_(settings.block_size * sizeof(ims::Peak) + BLOSC_MAX_OVERHEAD * 2),
      filled_(0),
      bytes_written_(0),
      filename_(filename),
      c_(settings),
      min_mz_(std::numeric_limits<double>::max()),
      max_mz_(std::numeric_limits<double>::min()) {
  if (c_.block_size < 512) throw std::runtime_error("block size is too small");
  if (blosc_compname_to_compcode(c_.compressor.c_str()) == -1)
    throw std::runtime_error("compressor '" + c_.compressor + "' is unavailable");
}

void imzb::ImzbWriter::dump() {
  if (filled_ == 0) return;
  int n = blosc_compress_ctx(int(c_.compression_level), 1, sizeof(ims::Peak),
      filled_ * sizeof(ims::Peak), &in_buf_[0], &out_buf_[0], out_buf_.size(),
      c_.compressor.c_str(), 0, 1);
  if (n <= 0) throw std::runtime_error("blosc compression error");
  filled_ = 0;
  out_.write(&out_buf_[0], n);
  index_.mzs.push_back(in_buf_[0].mz);

  min_mz_ = std::min(min_mz_, in_buf_.front().mz);
  max_mz_ = std::max(max_mz_, in_buf_.back().mz);

  index_.offsets.push_back(bytes_written_);
  bytes_written_ += n;
}

void imzb::ImzbWriter::writePeak(const ims::Peak& peak) {
  if (filled_ == in_buf_.size()) dump();
  in_buf_[filled_++] = peak;
}

void imzb::ImzbWriter::close() {
  dump();
  out_.close();

  std::ofstream out_idx(filename_ + ".idx", std::ios::binary);
  index_.header.version = imzb::VERSION;
  index_.header.block_size = c_.block_size;
  index_.header.min_mz = min_mz_;
  index_.header.max_mz = max_mz_;
  index_.write(out_idx);
  out_idx.close();
}

imzb::ImzbWriter::~ImzbWriter() {
  if (out_.is_open()) close();
}
