#include "cffi/common.hpp"
#include "ms/spectrum.hpp"
#include "ms/isocalc.hpp"

#include <stdexcept>
#include <new>

using namespace ms;
using namespace cffi;

extern "C" {

IMS_EXTERN InstrumentProfile* instrument_profile_new(
    const char* type, double resolving_power, double at_mz) {
  std::string stype{type};
  if (stype == "orbitrap")
    return new (std::nothrow) OrbitrapProfile{resolving_power, at_mz};
  else if (stype == "fticr")
    return new (std::nothrow) FTICRProfile(resolving_power, at_mz);
  else if (stype == "tof")
    return new (std::nothrow) TOFProfile(resolving_power);
  else {
    setErrorMessage("Unknown instrument type: '" + stype + "'");
    return nullptr;
  }
}

IMS_EXTERN void instrument_profile_free(InstrumentProfile* instrument) {
  delete instrument;
}

IMS_EXTERN double instrument_resolving_power_at(
    InstrumentProfile* instrument, double mz) {
  return instrument->resolvingPowerAt(mz);
}

IMS_EXTERN int spectrum_size(Spectrum* s) {
  return s->size();
}

IMS_EXTERN void spectrum_masses(Spectrum* s, double* out) {
  for (size_t i = 0; i < s->size(); i++)
    out[i] = s->masses.at(i);
}

IMS_EXTERN void spectrum_intensities(Spectrum* s, double* out) {
  for (size_t i = 0; i < s->size(); i++)
    out[i] = s->intensities.at(i);
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
  s->trimmed(n_peaks);
}

IMS_EXTERN void spectrum_trim_intensity(Spectrum* s, double min_intensity) {
  s->removeIntensitiesBelow(min_intensity);
}

IMS_EXTERN void spectrum_sort_by_mass(Spectrum* s) {
  s->sortByMass();
}

IMS_EXTERN void spectrum_sort_by_intensity(Spectrum* s) {
  s->sortByIntensity();
}

IMS_EXTERN void spectrum_free(Spectrum* s) {
  delete s;
}

IMS_EXTERN Spectrum* spectrum_new(int n, double* masses, double* abundances) {
  auto s = new (std::nothrow) Spectrum();
  s->masses.assign(masses, masses + n);
  s->intensities.assign(abundances, abundances + n);
  return s;
}

IMS_EXTERN Spectrum* spectrum_new_from_sf(
    const char* formula, double desired_probability) {
  return wrap_catch<Spectrum*>(nullptr, [&]() {
    return heapify(computeIsotopePattern(formula, desired_probability));
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

IMS_EXTERN float spectrum_envelope(Spectrum* s, InstrumentProfile* instr, double mz) {
  return s->envelope(instr, mz);
}

IMS_EXTERN int spectrum_envelope_plot(
    Spectrum* s, InstrumentProfile* instr, double* mzs, int n, float* out) {
  return wrap_catch<int>(-1, [&]() {
    ms::EnvelopeGenerator envelope(*s, instr);
    for (int i = 0; i < n; ++i)
      out[i] = envelope(mzs[i]);
    return 0;
  });
}

IMS_EXTERN void spectrum_normalize(Spectrum* s) {
  s->normalize();
}

IMS_EXTERN Spectrum* spectrum_envelope_centroids(
    Spectrum* s, InstrumentProfile* instr, double min_abundance, int points_per_fwhm) {
  return wrap_catch<Spectrum*>(nullptr, [&]() {
    return heapify(s->envelopeCentroids(instr, min_abundance, points_per_fwhm));
  });
}
}
