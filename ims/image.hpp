#pragma once

#include <valarray>
#include <algorithm>
#include <cassert>
#include <vector>

namespace ims {

inline size_t pixelIndex(size_t x, size_t y, size_t width) {
  return x * width + y;
}

template <typename T>
class Image {
  std::valarray<T> intensities_;
  size_t height_, width_;
  public:
  Image(size_t height, size_t width) :
    intensities_(height * width), height_(height), width_(width)
  {
  }

  size_t height() const { return height_; }
  size_t width() const { return width_; }

  const T& intensity(size_t x, size_t y) const {
    return intensities_[pixelIndex(x, y, width_)];
  }

  T& intensity(size_t x, size_t y) {
    return intensities_[pixelIndex(x, y, width_)];
  }

  T* rawPtr() { return &intensities_[0]; }

  size_t countEmptyPixels() const {
    size_t n = 0;
    for (auto& x: intensities_) if (x < -0.5) ++n;
    return n;
  }

  const std::valarray<T>& intensities() const { return intensities_; }

  std::pair<size_t, size_t> shape() const {
    return std::make_pair(height(), width());
  }

  void removeHotspots(double percentile, T* tmpbuf=nullptr) {
    assert(0 < percentile && percentile < 100);
    std::vector<T> tmp;
    if (tmpbuf == nullptr) {
      tmp.resize(height() * width());
      tmpbuf = &tmp[0];
    }
    T* pend = tmpbuf;
    for (T i: intensities_)
      if (i > 0)
        *pend++ = i;
    if (tmpbuf == pend)
      return; // empty image
    size_t k = size_t((pend - tmpbuf) * percentile) / 100;
    std::nth_element(tmpbuf, tmpbuf + k, pend);
    T threshold = tmpbuf[k];
    for (T& i: intensities_)
      if (i > threshold)
        i = threshold;
  }
};

typedef Image<float> ImageF;
typedef Image<double> ImageD;

}
