#pragma once

#include <vector>
#include <cstdlib>
#include <initializer_list>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <tuple>
#include <cassert>
#include <cmath>
#include <functional>

namespace ms {

enum class PeakOrder {
  mass,       // ascending mass
  intensity,  // descending intensity
  none        // unsorted
};

struct Spectrum {
  std::vector<double> masses;
  std::vector<double> intensities;

  Spectrum() {}
  Spectrum(double mass) : masses{mass}, intensities{1.0} {}
  Spectrum(
      std::initializer_list<double> masses, std::initializer_list<double> intensities)
      : masses{masses}, intensities{intensities} {}

  Spectrum(std::vector<double> masses, std::vector<double> intensities)
      : masses{masses}, intensities{intensities} {}

  // shifts masses accordingly to the number of added/subtracted electrons
  void addCharge(int charge);

  // Scales all intensities by a factor
  Spectrum& operator*=(double mult) {
    for (auto& m : intensities)
      m *= mult;
    return *this;
  }

  // Adds another spectrum. For performance reasons the spectrum becomes unsorted,
  // so that after adding multiple spectra everything can be sorted at once.
  Spectrum& operator+=(const Spectrum& other) {
    std::copy(other.masses.begin(), other.masses.end(), std::back_inserter(this->masses));
    std::copy(other.intensities.begin(), other.intensities.end(),
        std::back_inserter(this->intensities));
    peak_order_ = PeakOrder::none;
    return *this;
  }

  Spectrum& charged(int charge) {
    addCharge(charge);
    return *this;
  }

  // Returns intensity-sorted spectrum with new_size peaks
  Spectrum& trimmed(size_t new_size) {
    sortByIntensity();
    if (new_size < size()) {
      masses.resize(new_size);
      intensities.resize(new_size);
    }
    return *this;
  }

  // Removes peaks with intensity less than min_intensity.
  // Keeps peak sorting order.
  Spectrum& removeIntensitiesBelow(double min_intensity) {
    std::vector<double> ms;
    std::vector<double> is;
    size_t n = size();
    for (size_t i = 0; i < n; i++)
      if (intensities[i] >= min_intensity)
        ms.push_back(masses[i]), is.push_back(intensities[i]);
    masses.swap(ms);
    intensities.swap(is);
    return *this;
  }

  Spectrum copy() const {
    Spectrum p{masses, intensities};
    p.peak_order_ = this->peak_order_;
    return p;
  }

  // Sorts peaks by decreasing intensity and scales max intensity to 1
  Spectrum& normalize();

  // Normalized intensity-sorted convolution of two spectra.
  // Here intensities are treated as probabilities, and masses as isotopes, so that
  // for each peak pair (m1, i1) from the first spectrum and (m2, i2) from the second
  // the peak of the convolution is computed as (m1 + m2, i1 * i2).
  ms::Spectrum convolve(const Spectrum& other, double threshold = 0.0) const;

  bool isConvolutionUnit() const { return masses.size() == 1 && masses[0] == 0.0; }

  // Intensity of a gaussian-smoothed spectrum at a given mass.
  //
  // WIdth parameter controls how many sigmas on each peak side to treat as non-zero,
  // the default value is very conservative corresponding to highest precision.
  double envelope(double resolving_power, double mass, size_t width = 12) const;

  // Normalized intensity-sorted list of centroids
  ms::Spectrum envelopeCentroids(double resolving_power, double min_abundance = 1e-4,
      size_t points_per_fwhm = 25, size_t centroid_bins = 15) const;

  // sorts by decreasing intensity in-place
  ms::Spectrum& sortByIntensity(bool force = false);

  // sorts by increasing mass in-place
  ms::Spectrum& sortByMass(bool force = false);

  // Number of peaks
  size_t size() const { return masses.size(); }

 private:
  PeakOrder peak_order_ = PeakOrder::none;

  // sorts by last used order
  void forcedSort();
};

template <typename MzIt, typename IntIt>
std::pair<double, double> centroid(
    const MzIt& mzs_begin, const MzIt& mzs_end, const IntIt& intensities_begin) {
  double mz = 0.0;
  double total_intensity = 0.0;
  double max_intensity = 0.0;
  auto centroid_bins = (mzs_end - mzs_begin);
  MzIt mzs_center = mzs_begin + centroid_bins / 2;
  MzIt mzs_it = mzs_center;
  IntIt int_it = intensities_begin + (mzs_it - mzs_begin);
  while (mzs_it != mzs_begin && *(int_it - 1) < *int_it)
    --mzs_it, --int_it;
  for (; mzs_it != mzs_end; ++mzs_it, ++int_it) {
    if (mzs_it > mzs_center && mzs_it + 1 != mzs_end && *(int_it + 1) > *int_it) break;
    mz += (*mzs_it) * (*int_it), total_intensity += *int_it;
    max_intensity = std::max(max_intensity, double(*int_it));
  }
  mz /= total_intensity;
  return std::make_pair(mz, max_intensity);
}

template <typename MzIt, typename IntIt>
Spectrum detectPeaks(const MzIt& mzs_begin, const MzIt& mzs_end,
    const IntIt& intensities_begin, int window_size) {
  Spectrum result;
  if (window_size < 1) throw std::runtime_error("window size must be positive");
  if (window_size > 1000) throw std::runtime_error("window size is too large!");
  auto mzs_it = mzs_begin + window_size / 2;
  auto int_it = intensities_begin + window_size / 2;
  auto end = mzs_end - (window_size - window_size / 2);

  for (; mzs_it < end; ++mzs_it, ++int_it) {
    bool is_local_maximum = (*(int_it - 1) < *int_it) && (*int_it >= *(int_it + 1));

    if (!is_local_maximum) continue;

    double m, intensity;
    std::tie(m, intensity) = centroid(mzs_it - window_size / 2,
        mzs_it - window_size / 2 + window_size, int_it - window_size / 2);

    result.masses.push_back(m);
    result.intensities.push_back(intensity);
  }

  if (result.size() > 0) result.normalize();

  return result;
}

constexpr static double fwhm_to_sigma = 2.3548200450309493;  // 2 \sqrt{2 \log 2}

double sigmaAtResolution(const Spectrum& p, double resolution);

class EnvelopeGenerator {
  Spectrum p_;
  double resolution_;
  size_t width_;

  size_t peak_index_;  // index of the currently processed peak
  bool empty_space_;   // true if there are no peaks within distance (width *
                       // sigma)

  double last_mz_;
  double fwhm_, sigma_;

  double envelope(double mz);

 public:
  EnvelopeGenerator(const Spectrum& p, double resolution, size_t width = 12)
      : p_(p.copy().sortByMass()),
        resolution_(resolution),
        width_(width),
        peak_index_(0),
        empty_space_(false),
        last_mz_(-std::numeric_limits<double>::min()) {
    assert(p.size() > 0 && p.masses[0] > 0);
  }

  double operator()(double mz);
};
}
