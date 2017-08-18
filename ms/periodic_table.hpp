#pragma once

#include "ms/spectrum.hpp"

#include <map>
#include <vector>
#include <string>

namespace ms {

constexpr double electronMass = 0.00054857990924;

struct Element;
extern const std::map<std::string, ms::Element> periodic_table;

struct Element {
  std::string abbr;
  unsigned number;
  ms::Spectrum isotope_pattern;  // intensity-sorted

private:
  std::vector<double> log_probs_;

public:

  Element(const std::string& abbr, unsigned atomicNumber, const ms::Spectrum& isotope_pattern) :
    abbr{abbr}, number{atomicNumber}, isotope_pattern{isotope_pattern},
    log_probs_(isotope_pattern.size())
  {
    this->isotope_pattern.normalize();

    double total_intensity = 0.0;
    for (double x: this->isotope_pattern.intensities)
      total_intensity += x;

    for (int i = 0; i < this->isotope_pattern.size(); i++)
      log_probs_[i] = std::log(this->isotope_pattern.intensities[i] / total_intensity);
  }

  // log-probability of i-th most common isotope (counting from zero)
  double logProbability(int i) const {
    return log_probs_[i];
  }

  static bool isKnown(const std::string& abbr) {
    return periodic_table.find(abbr) != periodic_table.end();
  }

  static Element getByName(const std::string& abbr) {
    auto el = periodic_table.find(abbr)->second;
    return el;
  }
};
}
