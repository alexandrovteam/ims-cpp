#include "ims/image.hpp"
#include "ims/image_measures.hpp"

#include <cstdint>
#include <limits>
#include <cstring>
#include <valarray>
#include <numeric>

namespace ims {
double isotopeImageCorrelation(const ims::ImageF* images, size_t n,
    const std::vector<double>& abundances, bool zero_mean_normalization) {
  assert(n <= abundances.size());
  if (n < 2) return 0.0;

  auto n_empty_pixels = images[0].countEmptyPixels();

  auto cov = [&](const ims::ImageF& i1, const ims::ImageF& i2) -> double {
    // skip division as it cancels out in corr. coef. calculation
    double result = std::inner_product(std::begin(i1.intensities()),
        std::end(i1.intensities()), std::begin(i2.intensities()), double(0.0));

    if (zero_mean_normalization) {
      auto total1 = i1.intensities().sum() + n_empty_pixels;
      auto total2 = i2.intensities().sum() + n_empty_pixels;
      result -= total1 * total2 / (i1.intensities().size() - n_empty_pixels);
    }

    return result - n_empty_pixels;
  };

  auto principle_peak_norm = std::sqrt(cov(images[0], images[0]));
  constexpr double eps = 1e-6;
  if (std::fabs(principle_peak_norm) < eps) return 0.0;

  std::valarray<float> correlations(n - 1);
  for (size_t i = 1; i < n; ++i) {
    assert(images[i].shape() == images[0].shape());
    auto peak_norm = std::sqrt(cov(images[i], images[i]));
    if (std::fabs(peak_norm) < eps)
      correlations[i - 1] = 0.0;
    else
      correlations[i - 1] = cov(images[0], images[i]) / principle_peak_norm / peak_norm;
  }

  double weighted_sum = 0.0, total_weight = 0.0;
  for (size_t i = 0; i < correlations.size(); i++) {
    weighted_sum += correlations[i] * abundances[i + 1];
    total_weight += abundances[i + 1];
  }

  return weighted_sum / total_weight;
}

double isotopeImageCorrelation(const std::vector<ims::ImageF>& images,
    const std::vector<double>& abundances, bool zero_mean_normalization) {
  return isotopeImageCorrelation(
      &images[0], images.size(), abundances, zero_mean_normalization);
}

double isotopePatternMatch(
    const ims::ImageF* images, size_t n, const std::vector<double>& abundances) {
  assert(n <= abundances.size());
  std::valarray<float> total_intensities(n), expected(n);
  double norm_squared1 = 0.0, norm_squared2 = 0.0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < images[0].intensities().size(); j++) {
      if (images[0].intensities()[j] <= 0) continue;
      total_intensities[i] += images[i].intensities()[j];
    }
    norm_squared1 += std::pow(total_intensities[i], 2);
    expected[i] = abundances[i];
    norm_squared2 += std::pow(expected[i], 2);
  }
  if (std::fabs(norm_squared1) < 1e-6) return 0.0;
  total_intensities /= std::sqrt(norm_squared1);
  expected /= std::sqrt(norm_squared2);
  double score = 1.0;
  for (size_t i = 0; i < n; i++)
    score -= std::abs(expected[i] - total_intensities[i]) / n;
  return score;
}

double isotopePatternMatch(
    const std::vector<ims::ImageF>& images, const std::vector<double>& abundances) {
  return isotopePatternMatch(&images[0], images.size(), abundances);
}

double measureOfChaos(const ims::ImageF& image, size_t n_levels) {
  assert(n_levels > 0);
  assert(n_levels <= 32);
  std::vector<uint32_t> tmp1, tmp2, tmp3;

  tmp1.assign(image.width() * image.height(), 0);
  tmp2.resize(image.width() * image.height());
  tmp3.resize(image.width() * image.height());

  std::vector<float> levels(n_levels + 1);
  const auto& intensities = image.intensities();
  auto max = *std::max_element(std::begin(intensities), std::end(intensities));

  // FIXME: n_levels - 1 is only for compatibility with pyIMS code,
  //        it would be more sensible to divide by n_levels
  for (size_t i = 0; i < n_levels; i++)
    levels[i] = float(max * i) / (n_levels - 1);
  levels[n_levels] = std::numeric_limits<float>::max();

  // set bits indicating if i-th pixel is in j-th level set
  for (size_t i = 0; i < intensities.size(); i++)
    for (size_t j = 0; intensities[i] > levels[j]; j++)
      tmp1[i] <<= 1, tmp1[i] |= 1;

  size_t w = image.width(), h = image.height();

#define idx(x, y) ((x)*w + (y))

  // dilate with 4-connectivity 3x3 structuring element
  std::copy(tmp1.begin(), tmp1.end(), tmp2.begin());
  for (size_t y = 0; y < w; y++) {
    for (size_t x = 0; x < h; x++) {
      auto mask = tmp1[idx(x, y)];
      if (y > 0) tmp2[idx(x, y - 1)] |= mask;
      if (y < w - 1) tmp2[idx(x, y + 1)] |= mask;
      if (x > 0) tmp2[idx(x - 1, y)] |= mask;
      if (x < h - 1) tmp2[idx(x + 1, y)] |= mask;
    }
  }

  // erode with 8-connectivity 3x3 structuring element
  std::copy(tmp2.begin(), tmp2.end(), tmp1.begin());
  constexpr int d[8][2] = {
      {-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};

  for (size_t y = 0; y < w; y++)
    for (size_t x = 0; x < h; x++)
      for (size_t i = 0; i < 8; i++) {
        int nx = x + d[i][0], ny = y + d[i][1];
        if (nx >= 0 && nx < h && ny >= 0 && ny < w) tmp1[idx(x, y)] &= tmp2[idx(nx, ny)];
      }

  auto& parent = tmp3;  // union-find structure

  auto findRoot = [&](uint32_t i) -> uint32_t {
    uint32_t root = i;
    while (parent[root] < root)
      root = parent[root];
    return root;
  };

#define setRoot(i, r)                            \
  do {                                           \
    uint32_t tmp, s = i;                         \
    while (parent[s] < s)                        \
      tmp = parent[s], parent[s] = (r), s = tmp; \
    parent[s] = (r);                             \
  } while (0);

  auto mergeLabels = [&](uint32_t i, uint32_t j) -> uint32_t {
    auto r = findRoot(i);
    if (i != j) {
      r = std::min(r, findRoot(j));
      setRoot(j, r);
    }
    setRoot(i, r);
    return r;
  };

  uint32_t l = 1;
  auto& labels = tmp2;  // pixel labels

  uint32_t mask = 0;

#define isSet(x, y) (tmp1[idx((x), (y))] & mask)

  std::vector<float> counts(n_levels);

  // count numbers of connected components with 4-connectivity for each level
  // set
  for (size_t i = 0; i < n_levels; i++) {
    mask = 1UL << i;

    parent[0] = 0;
    memset(&labels[0], 0, labels.size() * sizeof(labels[0]));
    l = 1;

    for (size_t y = 0; y < w; y++)
      for (size_t x = 0; x < h; x++) {
        if (!isSet(x, y)) { continue; }
        bool up = x > 0 && isSet(x - 1, y);
        bool left = y > 0 && isSet(x, y - 1);
        if (up && left)
          labels[idx(x, y)] = mergeLabels(labels[idx(x - 1, y)], labels[idx(x, y - 1)]);
        else if (up)
          labels[idx(x, y)] = labels[idx(x - 1, y)];
        else if (left)
          labels[idx(x, y)] = labels[idx(x, y - 1)];
        else
          labels[idx(x, y)] = parent[l] = l, l += 1;
      }

    uint32_t n_components = 1;
    for (size_t i = 1; i < l; i++) {
      if (parent[i] < i)
        parent[i] = parent[parent[i]];
      else
        parent[i] = n_components++;
    }

    counts[i] = n_components - 1;
  }

  double avg_n_objects = 0.0;
  for (auto x : counts)
    avg_n_objects += x;
  avg_n_objects /= counts.size();

  double not_null = 0.0;
  for (auto x : intensities)
    not_null += x > 0;

  // FIXME: works badly on noisy images with few pixels:
  //        dilation + erosion leave almost nothing,
  //        that leads to small number of objects divided by
  //        large number of pixels

  return 1.0 - avg_n_objects / not_null;
}

double medianMzShift(const std::vector<ims::Peak>& peaks, double mz) {
  if (peaks.size() == 0) return 0;
  std::vector<double> shifts(peaks.size());
  for (auto& p : peaks)
    shifts.push_back(p.mz - mz);
  std::sort(shifts.begin(), shifts.end());
  if (peaks.size() % 2 == 1)
    return shifts[peaks.size() / 2];
  else
    return (shifts[peaks.size() / 2] + shifts[peaks.size() / 2 - 1]) / 2;
}
}
