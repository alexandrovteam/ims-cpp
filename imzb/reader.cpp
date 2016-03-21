#include "imzb/reader.hpp"

#include "blosc.h"

#include <ios>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <sstream>

imzb::ImzbReader::ImzbReader(const std::string& filename) :
    fn_(filename),
    in_(filename, std::ios::binary), block_idx_(0),
    n_peaks_(0), pos_(0), empty_(false)
{
  std::ifstream in_idx(filename + ".idx", std::ios::binary);
  if (!in_.is_open())
    throw std::runtime_error("couldn't open the file");
  if (!in_idx.is_open())
    throw std::runtime_error("couldn't open the index file");

  index_ = std::make_shared<Index>();
  index_->read(in_idx);

  in_.seekg(0, in_.end);
  index_->offsets.push_back(in_.tellg());

  in_.seekg(0, in_.beg);

  peaks_.resize(index_->header.block_size);
}

size_t imzb::ImzbReader::decompressBlock(size_t block_idx,
    std::ifstream& in,
    std::vector<char>& inbuf,
    std::vector<ims::Peak>& outbuf) const
{
  assert(block_idx + 1 < index_->offsets.size());
  uint64_t start = index_->offsets[block_idx];
  uint64_t end = index_->offsets[block_idx + 1];
  assert(start < end);
  inbuf.resize(end - start);
  in.seekg(start);
  in.read(&inbuf[0], end - start);
  int result = blosc_decompress_ctx(&inbuf[0], &outbuf[0],
      outbuf.size() * sizeof(ims::Peak), 1);
  assert(result > 0);
  return result / sizeof(ims::Peak);
}

bool imzb::ImzbReader::readNextBlock()
{
  if (block_idx_ == index_->mzs.size()) {
    n_peaks_ = 0;
    return false;
  }

  n_peaks_ = decompressBlock(block_idx_, in_, buffer_, peaks_);
  pos_ = 0;
  ++block_idx_;
  return true;
}

bool imzb::ImzbReader::readNext(ims::Peak& peak)
{
  if (empty_)
    return false;

  if (pos_ >= n_peaks_) {
    if (!readNextBlock()) {
      empty_ = true;
      return false;
    }
  }

  peak = peaks_[pos_];
  ++pos_;
  return true;
}

void imzb::ImzbReader::reset()
{
  in_.seekg(0, in_.beg);
  n_peaks_ = pos_ = block_idx_ = 0;
  empty_ = false;
}

std::vector<ims::Peak> imzb::ImzbReader::slice(double min_mz, double max_mz) const
{
  assert(min_mz < max_mz);
  std::vector<char> inbuf;
  std::vector<ims::Peak> result, outbuf{index_->header.block_size};
  size_t start_block = index_->startBlock(min_mz);
  size_t end_block = index_->endBlock(max_mz);
  std::ifstream in{fn_, std::ios::binary};

  for (size_t i = start_block; i < end_block; ++i) {
    size_t n = decompressBlock(i, in, inbuf, outbuf);
    auto beg = outbuf.cbegin(), end = outbuf.cbegin() + n;
    if (beg == end)
      continue;
    auto first_mz = beg->mz, last_mz = (end - 1)->mz;
    if (first_mz < min_mz)
      beg = std::lower_bound(beg, end, min_mz,
          [](const ims::Peak& p, double mz) { return p.mz < mz; });
    if (last_mz > max_mz)
      end = std::upper_bound(beg, end, max_mz,
          [](double mz, const ims::Peak& p) { return mz < p.mz; });
    result.insert(result.end(), beg, end);
  }
  return result;
}

ims::Image<float> imzb::ImzbReader::image(double mz, double ppm) const
{
  assert(ppm > 0);

  ims::Image<float> img(height(), width());
  readImage(mz, ppm, img.rawPtr());
  return img;
}

#define unlikely(x) __builtin_expect(!!(x), 0)

static void initImage(float* image, const imzb::Mask& mask) {
  for (size_t i = 0; i < mask.height; ++i)
    for (size_t j = 0; j < mask.width; ++j) {
      auto idx = ims::pixelIndex(i, j, mask.width);
      if (!mask.hasSpectrumAt(i, j))
        image[idx] = -1.0;
      else
        image[idx] = 0.0;
    }
}

void imzb::ImzbReader::readImage(double mz, double ppm, float* image) const
{
  initImage(image, index_->header.mask);

  auto peaks = slice(mz - mz * ppm * 1e-6, mz + mz * ppm * 1e-6);
  for (auto& peak: peaks) {
    if (unlikely(!index_->header.mask.hasSpectrumAt(peak.coords.x, peak.coords.y))) {
      std::stringstream ss;
      ss << "peak at x=" << peak.coords.x
         << ", y=" << peak.coords.y << ", m/z=" << peak.mz
         << " is outside the spatial mask";
      throw std::runtime_error(ss.str());
    }

    if (unlikely(peak.intensity < 0)) {
      std::stringstream ss;
      ss << "negative intensity peak at x=" << peak.coords.x
         << ", y=" << peak.coords.y << ", m/z=" << peak.mz;
      throw std::runtime_error(ss.str());
    }

    auto idx = ims::pixelIndex(peak.coords.x, peak.coords.y, width());
    image[idx] += peak.intensity;
  }
}

void imzb::ImzbReader::readCentroidedImage(double mz, double ppm, float* image) const
{
  initImage(image, index_->header.mask);

  auto peaks = slice(mz - mz * ppm * 1e-6, mz + mz * ppm * 1e-6);

  auto n = height() * width();
  auto float_inf = std::numeric_limits<float>::infinity();
  std::vector<double> m_prev(n), m_curr(n), m_next(n);
  std::vector<float> i_prev(n, float_inf), i_curr(n, float_inf), i_next(n, float_inf);

  for (auto& peak: peaks) {
    auto i = ims::pixelIndex(peak.coords.x, peak.coords.y, width());

    m_prev[i] = m_curr[i], i_prev[i] = i_curr[i];
    m_curr[i] = m_next[i], i_curr[i] = i_next[i];
    m_next[i] = peak.mz, i_next[i] = peak.intensity;

    if (i_prev[i] < i_curr[i] && i_curr[i] >= i_next[i]) {
      image[i] += i_curr[i];
    }
  }
}
