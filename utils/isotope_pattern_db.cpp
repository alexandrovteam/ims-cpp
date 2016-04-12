#include "utils/isotope_pattern_db.hpp"
#include "ms/isocalc.hpp"

#include "msgpack.hpp"

extern "C" {
#include "progressbar.h"
}

#include <cstdio>
#include <chrono>
#include <mutex>

namespace utils {

void IsotopePatternDB::computeIsotopePatterns(const utils::InstrumentProfile& instrument,
                                              size_t max_peaks)
{
  std::mutex map_mutex;

  progressbar* bar = nullptr;
  const int BAR_STEP = 100;
  if (use_progressbar_)
    bar = progressbar_new("", pairs_.size() / BAR_STEP);

  std::chrono::high_resolution_clock clock;

#pragma omp parallel for private(clock)
  for (size_t i = 0; i < pairs_.size(); i++) {
    const auto& f = pairs_[i].first;
    auto adduct = pairs_[i].second;
    int charge = adduct[0] == '-' ? -1 : +1;
    if (adduct[0] != '-' && adduct[0] != '+')
      adduct = "+" + adduct;
    auto full_formula = f + adduct;
    try {
      auto res = instrument.resolutionAt(ms::monoisotopicMass(full_formula));

      auto t1 = clock.now();
      auto pattern = ms::computeIsotopePattern(full_formula);

      auto t2 = clock.now();
      pattern = pattern.centroids(res).charged(charge).trimmed(max_peaks);

      auto t3 = clock.now();
      auto key = std::make_pair(f, adduct);

      std::lock_guard<std::mutex> lock(map_mutex);
      patterns_[key]["mzs"] = pattern.masses;
      patterns_[key]["abundances"] = pattern.abundances;
      patterns_[key]["time"] = std::vector<double>{
        std::chrono::duration<double, std::milli>(t2 - t1).count(),
        std::chrono::duration<double, std::milli>(t3 - t2).count()
      };
    } catch (sf_parser::NegativeTotalError) {
      // ignore formulas for which adduct subtraction is impossible
    }

    std::lock_guard<std::mutex> lock(map_mutex);
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
  if (db)
    db.write(sbuf.data(), sbuf.size());
  else
    throw std::runtime_error("can't open " + output_filename + " for writing");
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
