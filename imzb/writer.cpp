#include "imzb/writer.hpp"
#include "imzb/version.hpp"

#include "blosc.h"

#include <ios>
#include <cassert>
#include <stdexcept>

imzb::ImzbWriter::ImzbWriter(const std::string& filename,
                             uint32_t block_size,
                             const std::string& compressor,
                             uint8_t compression_level) :
    out_{filename, std::ios::binary},
    in_buf_(block_size),
    out_buf_(block_size * sizeof(ims::Peak) + BLOSC_MAX_OVERHEAD * 2),
    filled_(0), bytes_written_(0), filename_(filename),
    block_size_(block_size),
    compressor_(compressor),
    comp_level_(compression_level)
{
  if (block_size < 512)
    throw std::runtime_error("block size is too small");
  if (blosc_compname_to_compcode(compressor_.c_str()) == -1)
    throw std::runtime_error("compressor '" + compressor_ + "' is unavailable");
}

void imzb::ImzbWriter::dump()
{
  if (filled_ == 0)
    return;
  int n = blosc_compress_ctx(int(comp_level_), 1, sizeof(ims::Peak),
                             filled_ * sizeof(ims::Peak), &in_buf_[0],
                             &out_buf_[0], out_buf_.size(),
                             compressor_.c_str(), 0, 1);
  if (n <= 0)
    throw std::runtime_error("blosc compression error");
  filled_ = 0;
  out_.write(&out_buf_[0], n);
  index_.mzs.push_back(in_buf_[0].mz);
  index_.offsets.push_back(bytes_written_);
  bytes_written_ += n;
}

void imzb::ImzbWriter::writePeak(const ims::Peak& peak)
{
  if (filled_ == in_buf_.size())
    dump();
  in_buf_[filled_++] = peak;
}

void imzb::ImzbWriter::close()
{
  dump();
  out_.close();

  std::ofstream out_idx(filename_ + ".idx", std::ios::binary);
  index_.header.version = imzb::VERSION;
  index_.header.block_size = block_size_;
  index_.write(out_idx);
  out_idx.close();
}

imzb::ImzbWriter::~ImzbWriter()
{
  if (out_.is_open())
    close();
}
