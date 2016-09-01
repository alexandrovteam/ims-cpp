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
  ms::Spectrum isotope_pattern;

  static bool isKnown(const std::string& abbr) {
    return periodic_table.find(abbr) != periodic_table.end();
  }

  static Element getByName(const std::string& abbr) {
    auto el = periodic_table.find(abbr)->second;
    el.isotope_pattern.normalize();
    return el;
  }
};
}
