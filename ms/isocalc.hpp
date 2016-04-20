#pragma once

#include "ms/isotope_pattern.hpp"
#include "ms/periodic_table.hpp"

#include <string>
#include <sstream>
#include <map>
#include <cstdint>
#include <vector>

namespace ms {
  typedef std::map<std::string, uint16_t> ElementCounter;

  double monoisotopicMass(const ElementCounter& counter);
  double monoisotopicMass(const std::string& formula);

  IsotopePattern computeIsotopePattern(const Element& element, size_t amount, double fft_threshold = 0.0);

  IsotopePattern computeIsotopePattern(const ElementCounter& element_counts,
      double threshold = 1e-4, double fft_threshold = 1e-8);

  IsotopePattern computeIsotopePattern(const std::string& sum_formula,
      double threshold = 1e-4, double fft_threshold = 1e-8);

  namespace mass_search {

    struct ElementSettings {
      double monoisotopic_mass;
      size_t min_count, max_count;
      std::string name;

      ElementSettings(const std::string& name, // must be from the periodic table
                      size_t min_count = 0,
                      size_t max_count = 1000) :
        monoisotopic_mass(::ms::monoisotopicMass(name)),
        min_count{min_count}, max_count{max_count}, name{name}
      {
      }
    };

    typedef std::vector<size_t> ElementCounts;

    struct ExactMassSearch;

    struct Result {
      double mass;
      ElementCounts counts;
      const ExactMassSearch* const settings;

      Result(double, const ElementCounts&, const ExactMassSearch* const);

      std::string sumFormula() const;
    };

    /**
       Class for searching sum formulas that have a given monoisotopic mass.

       Example of usage:

       using ms::mass_search;
       ExactMassSearch search{512.436, 1.0,
           {{"C", 0, 100}, {"H", 0, 100}, {"O", 0, 100}, {"N", 0, 20}}
       };

       for (auto& r: search.run())
           std::cout << r.sumFormula() << std::endl;
     **/
    struct ExactMassSearch {
      ExactMassSearch(double mass, double ppm,
                      const std::vector<ElementSettings> elements,
                      size_t max_results=1000000);

      std::vector<Result> run() const;

      double mass;
      double ppm;
      std::vector<ElementSettings> elements;
      size_t max_results;
    };
  }
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
