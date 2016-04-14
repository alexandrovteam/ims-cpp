#pragma once

#include "ms/isotope_pattern.hpp"
#include "ms/periodic_table.hpp"

#include <string>
#include <sstream>
#include <map>
#include <cstdint>

namespace ms {
  typedef std::map<std::string, uint16_t> ElementCounter;

  double monoisotopicMass(const ElementCounter& counter);
  double monoisotopicMass(const std::string& formula);

  IsotopePattern computeIsotopePattern(const Element& element, size_t amount, double fft_threshold = 0.0);

  IsotopePattern computeIsotopePattern(const ElementCounter& element_counts,
      double threshold = 1e-4, double fft_threshold = 1e-8);

  IsotopePattern computeIsotopePattern(const std::string& sum_formula,
      double threshold = 1e-4, double fft_threshold = 1e-8);
}

namespace sf_parser {
  class ParseError : public std::exception {
    std::string s_;
    size_t pos_;
  public:
    ParseError(const std::string& str, size_t pos) : s_(str), pos_(pos) {}
    virtual const char* what() const noexcept {
      std::ostringstream ss;
      ss << s_ << " at position " << pos_;
      return ss.str().c_str();
    }
  };

  class NegativeTotalError : public std::exception {
    std::string element_;
    int total_;
  public:
    NegativeTotalError(const std::string& element, int total) : element_(element), total_(total) {}
    virtual const char* what() const noexcept {
      std::ostringstream ss;
      ss << "total number of " << element_ << " elements (" << total_ << ") is less than zero";
      return ss.str().c_str();
    }
  };

  ms::ElementCounter parseSumFormula(const std::string& formula);
}
