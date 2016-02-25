#include "ms/isotope_pattern.hpp"
#include "ms/periodic_table.hpp"

#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <array>

namespace ms {

void IsotopePattern::addCharge(int charge) {
  for (auto& m: masses)
    m -= charge * ms::electronMass;
}

IsotopePattern IsotopePattern::multiply(const IsotopePattern& other, double threshold) const {
  if (this->isUnit()) return other;
  if (other.isUnit()) return *this;

  const auto& p1 = *this, p2 = other;
  size_t n1 = p1.size(), n2 = p2.size();
  IsotopePattern result;
  for (size_t i = 0; i < n1; i++)
    for (size_t j = 0; j < n2; j++) {
      auto abundance = p1.abundances[i] * p2.abundances[j];
      if (abundance > threshold) {
        result.masses.push_back(p1.masses[i] + p2.masses[j]);
        result.abundances.push_back(abundance);
      }
    }
  result.normalize();
  return result;
}

void IsotopePattern::normalize() {
  for (size_t i = 0; i < size(); i++)
    for (size_t j = i + 1; j < size(); j++)
      if (abundances[i] < abundances[j]) {
        std::swap(abundances[i], abundances[j]);
        std::swap(masses[i], masses[j]);
      }
  auto top = abundances[0];
  assert(top > 0);
  for (auto& item: abundances)
    item /= top;
}

IsotopePattern sortByMass(const ms::IsotopePattern& p) {
  IsotopePattern result = p;
  for (size_t i = 0; i < result.size(); i++)
    for (size_t j = i + 1; j < result.size(); j++)
      if (result.masses[i] > result.masses[j]) {
        std::swap(result.abundances[i], result.abundances[j]);
        std::swap(result.masses[i], result.masses[j]);
      }
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
  size_t center = centroid_bins / 2;
  mz_window[center] = min_mz + pad - width * sigma;
  for (size_t j = 0; j < centroid_bins; j++) {
    mz_window[j] = mz_window[center] + (int(j) - int(center)) * step;
    int_window[j] = envelope(mz_window[j]);
  }

  ms::IsotopePattern result;
  for (;;) {
    std::rotate(mz_window.begin(), mz_window.begin() + 1, mz_window.end());
    std::rotate(int_window.begin(), int_window.begin() + 1, int_window.end());
    mz_window.back() = mz_window[centroid_bins - 2] + step;
    if (mz_window.back() > max_mz)
      break;
    int_window.back() = envelope(mz_window.back());

    // check if it's a local maximum
    if (!(int_window[center - 1] < int_window[center] &&
          int_window[center] >= int_window[center + 1]))
      continue;

    // skip low-intensity peaks
    if (int_window[center] < min_abundance)
      continue;

    double m, intensity;
    std::tie(m, intensity) = centroid(mz_window.begin(), mz_window.end(),
                                      int_window.begin());
    result.masses.push_back(m);
    result.abundances.push_back(intensity);
  }

  if (result.size() == 0)
    throw std::logic_error("the result contains no peaks, make min_abundance lower!");

  result.normalize();
  return result;
}


}
