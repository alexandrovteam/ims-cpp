#include "ms/isotope_pattern.hpp"
#include "ms/periodic_table.hpp"

#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
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

constexpr double fwhm_to_sigma = 2.3548200450309493; // 2 \sqrt{2 \log 2}
constexpr int width = 12; // defines how many sigmas to each side from the peak we consider
constexpr size_t centroid_bins = 15;

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

double IsotopePattern::envelope(double resolution, double mz) const {
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

IsotopePattern IsotopePattern::centroids(double resolution, double min_abundance, size_t points_per_fwhm) const {
  if (this->isUnit() || this->masses.size() == 0)
    return *this;
  auto p = sortByMass(*this);
  assert(this->masses[0] > 0);
  if (points_per_fwhm < 5)
    throw std::logic_error("points_per_fwhm must be at least 5 for meaningful results");
  double fwhm = p.masses[0] / resolution;
  double step = fwhm / points_per_fwhm;
  double sigma = fwhm / fwhm_to_sigma;
  ms::IsotopePattern result;

  size_t peak_index = 0; // one of the two peaks closest to the current center
  bool empty = false; // true if there are no peaks within distance (width * sigma)

  auto envelope = [&](double mz) -> double {
    if (empty) return 0.0; // no isotopic peaks nearby
    double result = 0.0;
    int k = peak_index, n = p.size();
    while (k > 0 && mz - p.masses[k - 1] < width * sigma)
      --k;
    while (k < n) {
      if (p.masses[k] - mz > width * sigma)
        break;
      result += p.abundances[k] * std::exp(-0.5 * std::pow((p.masses[k] - mz) / sigma, 2));
      ++k;
    }
    return result;
  };

  std::array<double, centroid_bins> mz_window, int_window;
  size_t center = centroid_bins / 2;
  mz_window[center] = p.masses[0] - width * sigma;
  for (size_t j = 0; j < centroid_bins; j++) {
    mz_window[j] = mz_window[center] + (int(j) - int(center)) * step;
    int_window[j] = envelope(mz_window[j]);
  }

  for (;;) {
    std::rotate(mz_window.begin(), mz_window.begin() + 1, mz_window.end());
    std::rotate(int_window.begin(), int_window.begin() + 1, int_window.end());
    mz_window.back() = mz_window[centroid_bins - 2] + step;
    int_window.back() = envelope(mz_window.back());

    if (mz_window[center] > p.masses[peak_index] + width * sigma)
      empty = true;

    if (empty && peak_index + 1 == p.size())
      break;

    if (empty && peak_index + 1 < p.size() &&
        mz_window[center] >= p.masses[peak_index + 1] - width * sigma)
      empty = false, ++peak_index;

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
