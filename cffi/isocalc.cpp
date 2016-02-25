#include "cffi/common.hpp"
#include "ms/isotope_pattern.hpp"
#include "ms/isocalc.hpp"

#include <stdexcept>
#include <new>

using namespace ms;
using namespace cffi;

extern "C" {

  IMS_EXTERN int isotope_pattern_size(IsotopePattern* p) { return p->size(); }

  IMS_EXTERN void isotope_pattern_masses(IsotopePattern* p, double* out) {
    for (size_t i = 0; i < p->size(); i++)
      out[i] = p->masses.at(i);
  }

  IMS_EXTERN void isotope_pattern_abundances(IsotopePattern* p, double* out) {
    for (size_t i = 0; i < p->size(); i++)
      out[i] = p->abundances.at(i);
  }

  IMS_EXTERN void isotope_pattern_add_charge(IsotopePattern* p, int charge) {
    p->addCharge(charge);
  }

  IMS_EXTERN void isotope_pattern_trim(IsotopePattern* p, unsigned n_peaks) {
    if (p->size() > n_peaks) {
      p->masses.resize(n_peaks);
      p->abundances.resize(n_peaks);
    }
  }

  IMS_EXTERN void isotope_pattern_free(IsotopePattern* p) { delete p; }

  IMS_EXTERN IsotopePattern* isotope_pattern_new(int n, double* masses, double* abundances) {
    auto p = new(std::nothrow) IsotopePattern();
    p->masses.assign(masses, masses + n);
    p->abundances.assign(abundances, abundances + n);
    return p;
  }

  IMS_EXTERN IsotopePattern* isotope_pattern_new_from_sf(const char* formula,
     double threshold, double fft_threshold)
  {
    return wrap_catch<IsotopePattern*>(nullptr, [&]() {
        return heapify(computeIsotopePattern(formula, threshold, fft_threshold));
    });
  }

  IMS_EXTERN IsotopePattern* isotope_pattern_new_from_raw(int n, double* masses, float* intensities, int window_size)
  {
    return wrap_catch<IsotopePattern*>(nullptr, [&]() {
        return heapify(detectPeaks(masses, masses + n, intensities, window_size));
    });
  }

  IMS_EXTERN IsotopePattern* isotope_pattern_copy(IsotopePattern* p) {
    return new(std::nothrow) IsotopePattern(*p);
  }

  IMS_EXTERN float isotope_pattern_envelope(IsotopePattern* p, double resolution, double mz)
  {
    return p->envelope(resolution, mz);
  }

  IMS_EXTERN void isotope_pattern_envelope_plot(IsotopePattern* p, double resolution,
                                                double* mzs, int n, float* out)
  {
    for (int i = 0; i < n; ++i)
      out[i] = p->envelope(resolution, mzs[i]);
  }

  IMS_EXTERN IsotopePattern* isotope_pattern_centroids(IsotopePattern* p, double resolution,
                                                       double min_abundance, int points_per_fwhm)
  {
    return wrap_catch<IsotopePattern*>(nullptr, [&]() {
      return heapify(p->centroids(resolution, min_abundance, points_per_fwhm));
    });
  }
}
