#include "ms/isotope_pattern.hpp"
#include "ms/periodic_table.hpp"

#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <array>
#include <functional>
#include <numeric>

namespace ms {

void IsotopePattern::addCharge(int charge) {
  for (auto& m: masses)
    m -= charge * ms::electronMass;
}

static bool isNormalized(const IsotopePattern& p) {
  if (p.abundances[0] != 1.0)
    return false;
  for (size_t i = 1; i < p.size(); i++)
    if (p.abundances[i - 1] < p.abundances[i])
      return false;
  return true;
}

IsotopePattern IsotopePattern::multiply(const IsotopePattern& other, double threshold) const {
  if (this->isUnit()) return other;
  if (other.isUnit()) return *this;

  const auto& p1 = *this, p2 = other;
  size_t n1 = p1.size(), n2 = p2.size();

  assert(isNormalized(p1));
  assert(isNormalized(p2));

  IsotopePattern result;

  for (size_t i = 0; i < n1; i++)
    for (size_t j = 0; j < n2; j++) {
      auto abundance = p1.abundances[i] * p2.abundances[j];
      if (abundance > threshold) {
        result.masses.push_back(p1.masses[i] + p2.masses[j]);
        result.abundances.push_back(abundance);
      } else {
        if (j == 0)
          break;
        n2 = j;
      }
    }
  result.normalize();
  return result;
}

static void sortIsotopePattern(std::vector<double>& masses, std::vector<double>& abundances,
                               std::function<bool(size_t, size_t)> cmp)
{
  auto n = masses.size();
  std::vector<size_t> index(masses.size());
  std::iota(index.begin(), index.end(), 0);
  std::sort(index.begin(), index.end(), cmp);
  std::vector<double> tmp(n);
  for (size_t i = 0; i < n; i++)
    tmp[i] = masses[index[i]];
  masses.swap(tmp);

  for (size_t i = 0; i < n; i++)
    tmp[i] = abundances[index[i]];
  abundances.swap(tmp);
}

void IsotopePattern::normalize() {
  sortIsotopePattern(masses, abundances, [&](size_t i, size_t j) {
    return abundances[i] > abundances[j];
  });
  auto top = abundances[0];
  if (top > 0)
    for (auto& item: abundances)
      item /= top;
}

IsotopePattern sortByMass(const ms::IsotopePattern& p) {
  IsotopePattern result = p;
  sortIsotopePattern(result.masses, result.abundances, [&](size_t i, size_t j) {
    return result.masses[i] < result.masses[j];
  });
  return result;
}

double IsotopePattern::envelope(double resolution, double mz, size_t width) const {
  double result = 0.0;

  double fwhm = masses[0] / resolution;
  double sigma = fwhm / fwhm_to_sigma;

  for (size_t k = 0; k < size(); ++k) {
    if (std::fabs(masses[k] - mz) > width * sigma)
      continue;
    result += abundances[k] * std::exp(-0.5 * std::pow((masses[k] - mz) / sigma, 2));
  }
  return result;
}

double EnvelopeGenerator::envelope(double mz) {
  if (empty_space_) return 0.0; // no isotopic peaks nearby
  double result = 0.0;
  int k = peak_index_, n = p_.size();
  while (k > 0 && mz - p_.masses[k - 1] < width_ * sigma_)
    --k;
  while (k < n && p_.masses[k] - mz <= width_ * sigma_) {
    double sd_dist = (p_.masses[k] - mz) / sigma_;
    result += p_.abundances[k] * std::exp(-0.5 * std::pow(sd_dist, 2));
    ++k;
  }
  return result;
}

double EnvelopeGenerator::operator()(double mz) {
  fwhm_ = mz / resolution_;
  sigma_ = fwhm_ / fwhm_to_sigma;

  if (last_mz_ > mz)
    throw std::runtime_error("input to EnvelopeGenerator must be sorted");

  if (mz > p_.masses[peak_index_] + width_ * sigma_)
    empty_space_ = true; // the last peak's influence is now considered negligible

  if (empty_space_ && peak_index_ + 1 < p_.size() &&
      mz >= p_.masses[peak_index_ + 1] - width_ * sigma_)
  {
    empty_space_ = false;
    ++peak_index_; // reached next peak
  }

  last_mz_ = mz;
  return envelope(mz);
}

constexpr size_t centroid_bins = 15;

IsotopePattern IsotopePattern::centroids(double resolution, double min_abundance, size_t points_per_fwhm) const {
  if (this->isUnit() || this->masses.size() == 0)
    return *this;
  if (points_per_fwhm < 5)
    throw std::logic_error("points_per_fwhm must be at least 5 for meaningful results");
  double fwhm = this->masses[0] / resolution;
  double sigma = fwhm / fwhm_to_sigma;
  const size_t width = 12;

  double step = fwhm / points_per_fwhm;
  double pad = step * centroid_bins / 2 + width * sigma;
  double min_mz = *std::min_element(this->masses.begin(), this->masses.end()) - pad;
  double max_mz = *std::max_element(this->masses.begin(), this->masses.end()) + pad;

  EnvelopeGenerator envelope(*this, resolution, width);

  std::array<double, centroid_bins> mz_window, int_window;
  size_t center = centroid_bins / 2, last_idx = centroid_bins - 1;
  mz_window[center] = min_mz + pad - width * sigma;
  for (size_t j = 0; j < centroid_bins; j++) {
    mz_window[j] = mz_window[center] + (int(j) - int(center)) * step;
    int_window[j] = envelope(mz_window[j]);
  }

  auto prev = [&](size_t idx) { return idx > 0 ? idx - 1 : centroid_bins - 1; };
  auto next = [&](size_t idx) { return idx == centroid_bins - 1 ? 0 : idx + 1; };
  auto rotateWindows = [&](int shift) -> void {
    auto abs_shift = (shift + int(centroid_bins)) % centroid_bins;
    std::rotate(mz_window.begin(),
                mz_window.begin() + abs_shift, mz_window.end());
    std::rotate(int_window.begin(),
                int_window.begin() + abs_shift, int_window.end());
  };

  ms::IsotopePattern result;
  for (;;) {
    center = next(center);

    auto next_mz = mz_window[last_idx] + step;
    if (next_mz > max_mz)
      break;

    last_idx = next(last_idx);
    mz_window[last_idx] = next_mz;
    int_window[last_idx] = envelope(next_mz);

    // check if it's a local maximum
    if (!(int_window[prev(center)] < int_window[center] &&
          int_window[center] >= int_window[next(center)]))
      continue;

    // skip low-intensity peaks
    if (int_window[center] < min_abundance)
      continue;

    double m, intensity;
    auto shift = int(centroid_bins - 1 - last_idx);
    rotateWindows(-shift);
    std::tie(m, intensity) = centroid(mz_window.begin(), mz_window.end(),
                                      int_window.begin());
    rotateWindows(shift);
    result.masses.push_back(m);
    result.abundances.push_back(intensity);
  }

  if (result.size() == 0)
    throw std::logic_error("the result contains no peaks, make min_abundance lower!");

  result.normalize();
  return result;
}


}
