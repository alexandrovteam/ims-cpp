#pragma once

#include "ms/isotope_pattern.hpp"

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>

namespace utils {

class InstrumentProfile {
public:
  virtual double resolutionAt(double mz) const = 0;
  virtual ~InstrumentProfile() {}
};

class OrbitrapProfile : public InstrumentProfile {
  double resolution_200_;
public:
  OrbitrapProfile(double resolution_at_200) :
    resolution_200_(resolution_at_200)
  {}

  double resolutionAt(double mz) const {
    return resolution_200_ * (200.0 / mz);
  }

  OrbitrapProfile() {
  }
};

class IsotopePatternDB {
  typedef std::pair<std::string, std::string> SFAdductPair;
  std::vector<SFAdductPair> pairs_;

  // use a simple structure that doesn't require special msgpack methods
  std::map<SFAdductPair, std::map<std::string, std::vector<double>>> patterns_;

  bool use_progressbar_;

  static std::vector<SFAdductPair> makePairs(
    const std::vector<std::string>& sum_formulas,
    const std::vector<std::string>& adducts)
  {
    std::vector<SFAdductPair> pairs;
    for (auto& f: sum_formulas)
      for (auto& a: adducts)
        pairs.push_back(std::make_pair(f, a));
    return pairs;
  }

public:

  void save(const std::string& output_filename);
  void load(const std::string& input_filename);

  IsotopePatternDB(const std::vector<std::string>& sum_formulas,
                   const std::vector<std::string>& adducts)
    : IsotopePatternDB(makePairs(sum_formulas, adducts))
  {
  }

  IsotopePatternDB(const std::vector<std::pair<std::string, std::string>>& sf_adduct_pairs) :
    pairs_(sf_adduct_pairs), use_progressbar_(false)
  {
  }

  IsotopePatternDB(const std::string& dump_filename) {
    load(dump_filename);
  }

  void computeIsotopePatterns(const utils::InstrumentProfile& instrument,
                              size_t max_peaks=5);

  ms::IsotopePattern operator()(const std::string& formula, const std::string& adduct) const {
    auto& entry = patterns_.at(std::make_pair(formula, adduct));
    return ms::IsotopePattern{entry.at("mzs"), entry.at("abundances")};
  }

  void useProgressBar(bool use) {
    use_progressbar_ = use;
  } 

  size_t size() const {
    assert(patterns_.size() == pairs_.size());
    return patterns_.size();
  }

  const std::vector<std::pair<std::string, std::string>>& keys() const {
    return pairs_;
  }
};

}
