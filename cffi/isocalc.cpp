#include "cffi/common.hpp"
#include "ms/spectrum.hpp"
#include "ms/isocalc.hpp"

#include <stdexcept>
#include <new>

using namespace ms;
using namespace cffi;

extern "C" {

IMS_EXTERN int spectrum_size(Spectrum* s) {
  return s->size();
}

IMS_EXTERN void spectrum_masses(Spectrum* s, double* out) {
  for (size_t i = 0; i < s->size(); i++)
    out[i] = s->masses.at(i);
}

IMS_EXTERN void spectrum_abundances(Spectrum* s, double* out) {
  for (size_t i = 0; i < s->size(); i++)
    out[i] = s->abundances.at(i);
}

IMS_EXTERN void spectrum_add_charge(Spectrum* s, int charge) {
  s->addCharge(charge);
}

IMS_EXTERN void spectrum_multiply_inplace(Spectrum* s, double mult) {
  s->operator*=(mult);
}

IMS_EXTERN void spectrum_add_inplace(Spectrum* s1, Spectrum* s2) {
  (*s1) += *s2;
}

IMS_EXTERN void spectrum_trim(Spectrum* s, unsigned n_peaks) {
  if (s->size() > n_peaks) {
    s->masses.resize(n_peaks);
    s->abundances.resize(n_peaks);
  }
}

IMS_EXTERN void spectrum_free(Spectrum* s) {
  delete s;
}

IMS_EXTERN Spectrum* spectrum_new(int n, double* masses, double* abundances) {
  auto s = new (std::nothrow) Spectrum();
  s->masses.assign(masses, masses + n);
  s->abundances.assign(abundances, abundances + n);
  return s;
}

IMS_EXTERN Spectrum* spectrum_new_from_sf(
    const char* formula, double threshold, double fft_threshold) {
  return wrap_catch<Spectrum*>(nullptr, [&]() {
    return heapify(computeIsotopePattern(formula, threshold, fft_threshold));
  });
}

IMS_EXTERN Spectrum* spectrum_new_from_raw(
    int n, double* masses, float* intensities, int window_size) {
  return wrap_catch<Spectrum*>(nullptr, [&]() {
    return heapify(detectPeaks(masses, masses + n, intensities, window_size));
  });
}

IMS_EXTERN Spectrum* spectrum_copy(Spectrum* s) {
  return new (std::nothrow) Spectrum(*s);
}

IMS_EXTERN float spectrum_envelope(Spectrum* s, double resolution, double mz) {
  return s->envelope(resolution, mz);
}

IMS_EXTERN int spectrum_envelope_plot(
    Spectrum* s, double resolution, double* mzs, int n, float* out) {
  return wrap_catch<int>(-1, [&]() {
    ms::EnvelopeGenerator envelope(*s, resolution);
    for (int i = 0; i < n; ++i)
      out[i] = envelope(mzs[i]);
    return 0;
  });
}

IMS_EXTERN void spectrum_normalize(Spectrum* s) {
  s->normalize();
}

IMS_EXTERN Spectrum* spectrum_envelope_centroids(
    Spectrum* s, double resolution, double min_abundance, int points_per_fwhm) {
  return wrap_catch<Spectrum*>(nullptr, [&]() {
    return heapify(s->envelopeCentroids(resolution, min_abundance, points_per_fwhm));
  });
}
}
