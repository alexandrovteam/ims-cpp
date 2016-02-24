#include "utils/isotope_pattern_db.hpp"
#include "ms/isocalc.hpp"

#include "msgpack.hpp"

extern "C" {
#include "progressbar.h"
}

#include <cstdio>
#include <mutex>

namespace utils {

void IsotopePatternDB::computeIsotopePatterns(const utils::InstrumentProfile& instrument,
                                              size_t max_peaks)
{
  std::mutex map_mutex;

  progressbar* bar = nullptr;
  const int BAR_STEP = 50;
  if (use_progressbar_)
    bar = progressbar_new("", pairs_.size() / BAR_STEP);

#pragma omp parallel for
  for (size_t i = 0; i < pairs_.size(); i++) {
    const auto& f = pairs_[i].first;
    auto adduct = pairs_[i].second;
    int charge = adduct[0] == '-' ? -1 : +1;
    if (adduct[0] != '-' && adduct[0] != '+')
      adduct = "+" + adduct;
    auto full_formula = f + adduct;
    auto res = instrument.resolutionAt(ms::monoisotopicMass(full_formula));
    auto pattern = ms::computeIsotopePattern(full_formula)\
      .centroids(res).charged(charge).trimmed(max_peaks);
    auto key = std::make_pair(f, adduct);

    std::lock_guard<std::mutex> lock(map_mutex);
    patterns_[key]["mzs"] = pattern.masses;
    patterns_[key]["abundances"] = pattern.abundances;
    if (bar != nullptr && (i + 1) % BAR_STEP == 0)
      progressbar_inc(bar);
  }

  if (bar != nullptr) {
    progressbar_finish(bar);
    std::fflush(stdout);
  }
}

void IsotopePatternDB::save(const std::string& output_filename) {
  // TODO: write instrument profile parameters?
  msgpack::sbuffer sbuf;
  msgpack::pack(sbuf, patterns_);
  std::ofstream db(output_filename, std::ios::binary);
  db.write(sbuf.data(), sbuf.size());
}

void IsotopePatternDB::load(const std::string& input_filename)
{
  // load precomputed isotope patterns from msgpack file
  std::ifstream in(input_filename, std::ios::binary | std::ios::ate);
  size_t bufsize = in.tellg();
  msgpack::sbuffer sbuf(bufsize);
  in.seekg(0, std::ios::beg);
  in.read(sbuf.data(), bufsize);
  in.close();

  msgpack::unpacked unpacked;
  msgpack::unpack(&unpacked, sbuf.data(), bufsize);

  msgpack::object obj = unpacked.get();
  obj.convert(patterns_);
  for (auto& item: patterns_) pairs_.push_back(item.first);
}

}
