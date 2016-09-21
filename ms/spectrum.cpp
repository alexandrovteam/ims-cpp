#include "ms/spectrum.hpp"
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

void Spectrum::addCharge(int charge) {
  for (auto& m : masses)
    m -= charge * ms::electronMass;
}

static bool isNormalized(const Spectrum& p) {
  if (p.intensities[0] != 1.0) return false;
  for (size_t i = 1; i < p.size(); i++)
    if (p.intensities[i - 1] < p.intensities[i]) return false;
  return true;
}

Spectrum Spectrum::convolve(const Spectrum& other, double threshold) const {
  if (this->isConvolutionUnit()) return other;
  if (other.isConvolutionUnit()) return *this;

  const auto& p1 = *this, p2 = other;
  size_t n1 = p1.size(), n2 = p2.size();

  assert(isNormalized(p1));
  assert(isNormalized(p2));

  Spectrum result;

  for (size_t i = 0; i < n1; i++)
    for (size_t j = 0; j < n2; j++) {
      auto abundance = p1.intensities[i] * p2.intensities[j];
      if (abundance > threshold) {
        result.masses.push_back(p1.masses[i] + p2.masses[j]);
        result.intensities.push_back(abundance);
      } else {
        if (j == 0) break;
        n2 = j;
      }
    }
  return result.normalize();
}

static void sortSpectrum(std::vector<double>& masses, std::vector<double>& intensities,
    std::function<bool(size_t, size_t)> cmp) {
  auto n = masses.size();
  std::vector<size_t> index(masses.size());
  std::iota(index.begin(), index.end(), 0);
  std::sort(index.begin(), index.end(), cmp);
  std::vector<double> tmp(n);
  for (size_t i = 0; i < n; i++)
    tmp[i] = masses[index[i]];
  masses.swap(tmp);

  for (size_t i = 0; i < n; i++)
    tmp[i] = intensities[index[i]];
  intensities.swap(tmp);
}

Spectrum& Spectrum::sortByIntensity(bool force) {
  if (peak_order_ == PeakOrder::intensity && !force) return *this;
  sortSpectrum(masses, intensities,
      [&](size_t i, size_t j) { return intensities[i] > intensities[j]; });
  peak_order_ = PeakOrder::intensity;
  return *this;
}

Spectrum& Spectrum::sortByMass(bool force) {
  if (peak_order_ == PeakOrder::mass && !force) return *this;
  sortSpectrum(
      masses, intensities, [&](size_t i, size_t j) { return masses[i] < masses[j]; });
  peak_order_ = PeakOrder::mass;
  return *this;
}

Spectrum& Spectrum::normalize() {
  sortByIntensity();
  auto top = intensities[0];
  if (top > 0)
    for (auto& item : intensities)
      item /= top;
  return *this;
}

static double sigmaAtResolvingPower(double mass, double resolving_power) {
  if (resolving_power <= 0) return NAN;
  auto fwhm = mass / resolving_power;
  return fwhm / fwhm_to_sigma;
}

static double calculateSigmaAt(double mass, const InstrumentProfile* instrument) {
  double resolving_power = instrument->resolvingPowerAt(mass);
  return sigmaAtResolvingPower(mass, resolving_power);
}

double Spectrum::envelope(
    const InstrumentProfile* instrument, double mass, size_t width) const {
  double result = 0.0;

  double sigma = calculateSigmaAt(mass, instrument);

  for (size_t k = 0; k < size(); ++k) {
    if (std::fabs(masses[k] - mass) > width * sigma) continue;
    result += intensities[k] * std::exp(-0.5 * std::pow((masses[k] - mass) / sigma, 2));
  }
  return result;
}

EnvelopeGenerator::EnvelopeGenerator(
    const Spectrum& p, const InstrumentProfile* instrument, size_t width)
    : p_(p.copy().sortByMass()),
      instrument_(instrument),
      width_(width),
      peak_index_(0),
      empty_space_(false),
      last_mz_(-std::numeric_limits<double>::min()) {
  assert(p.size() > 0 && p.masses[0] > 0);
  sigma_ = calculateSigmaAt(p_.masses[0], instrument_);
}

double EnvelopeGenerator::currentSigma() const {
  return sigma_;
}

double EnvelopeGenerator::envelope(double mz) {
  if (empty_space_) return 0.0;  // no isotopic peaks nearby
  double result = 0.0;
  int k = peak_index_, n = p_.size();
  while (k > 0 && mz - p_.masses[k - 1] < width_ * sigma_)
    --k;
  while (k < n && p_.masses[k] - mz <= width_ * sigma_) {
    double sd_dist = (p_.masses[k] - mz) / sigma_;
    result += p_.intensities[k] * std::exp(-0.5 * std::pow(sd_dist, 2));
    ++k;
  }
  return result;
}

double EnvelopeGenerator::operator()(double mz) {
  if (last_mz_ > mz)
    throw std::runtime_error("input to EnvelopeGenerator must be sorted");

  if (mz > p_.masses[peak_index_] + width_ * sigma_)
    empty_space_ = true;  // the last peak's influence is now considered negligible

  if (empty_space_ && peak_index_ + 1 < p_.size() &&
      mz >= p_.masses[peak_index_ + 1] - width_ * sigma_) {
    empty_space_ = false;
    ++peak_index_;  // reached next peak

    sigma_ = calculateSigmaAt(p_.masses[peak_index_], instrument_);
  }

  last_mz_ = mz;
  return envelope(mz);
}

Spectrum Spectrum::envelopeCentroids(const InstrumentProfile* instrument,
    double min_abundance, size_t points_per_fwhm, size_t centroid_bins) const {
  if (this->masses.size() <= 1) return *this;
  if (points_per_fwhm < 5)
    throw std::logic_error("points_per_fwhm must be at least 5 for meaningful results");
  if (centroid_bins < 3) throw std::logic_error("centroid_bins must be at least 3");

  const size_t width = 12;
  double min_mz = *std::min_element(this->masses.begin(), this->masses.end());
  double max_mz = *std::max_element(this->masses.begin(), this->masses.end());

  EnvelopeGenerator envelope(*this, instrument, width);

  double sigma = envelope.currentSigma();
  double step = (sigma * fwhm_to_sigma) / points_per_fwhm;

  int n_steps = 0;

  std::vector<double> mz_window(centroid_bins), int_window(centroid_bins);
  size_t center = centroid_bins / 2, last_idx = centroid_bins - 1;
  mz_window[center] = min_mz - width * sigma;
  for (size_t j = 0; j < centroid_bins; j++) {
    mz_window[j] = mz_window[center] + (int(j) - int(center)) * step;
    int_window[j] = envelope(mz_window[j]);
    ++n_steps;
  }

  auto prev = [&](size_t idx) { return idx > 0 ? idx - 1 : centroid_bins - 1; };
  auto next = [&](size_t idx) { return idx == centroid_bins - 1 ? 0 : idx + 1; };
  auto rotateWindows = [&](int shift) -> void {
    auto abs_shift = (shift + int(centroid_bins)) % centroid_bins;
    std::rotate(mz_window.begin(), mz_window.begin() + abs_shift, mz_window.end());
    std::rotate(int_window.begin(), int_window.begin() + abs_shift, int_window.end());
  };

  ms::Spectrum result;
  for (;;) {
    center = next(center);

    auto next_mz = mz_window[last_idx] + step;
    if (next_mz > max_mz + width * sigma) break;

    last_idx = next(last_idx);
    mz_window[last_idx] = next_mz;
    int_window[last_idx] = envelope(next_mz);

    sigma = envelope.currentSigma();
    step = (sigma * fwhm_to_sigma) / points_per_fwhm;
    ++n_steps;

    if (n_steps > 100000000)
      throw std::runtime_error("error in envelopeCentroids calculation");

    // check if it's a local maximum
    if (!(int_window[prev(center)] < int_window[center] &&
            int_window[center] >= int_window[next(center)]))
      continue;

    // skip low-intensity peaks
    if (int_window[center] < min_abundance) continue;

    double m, intensity;
    auto shift = int(centroid_bins - 1 - last_idx);
    rotateWindows(-shift);
    std::tie(m, intensity) =
        centroid(mz_window.begin(), mz_window.end(), int_window.begin());
    rotateWindows(shift);
    result.masses.push_back(m);
    result.intensities.push_back(intensity);
  }

  if (result.size() == 0)
    throw std::logic_error("the result contains no peaks, make min_abundance lower!");

  return result.sortByIntensity().normalize().removeIntensitiesBelow(min_abundance);
}
}
