#pragma once

#include "ms/spectrum.hpp"
#include "ms/periodic_table.hpp"

#include <string>
#include <sstream>
#include <map>
#include <cstdint>
#include <vector>
#include <memory>

namespace ms {
typedef std::map<std::string, int16_t> ElementCounter;

double monoisotopicMass(const ElementCounter& counter);
double monoisotopicMass(const std::string& formula);

Spectrum computeIsotopePattern(const ElementCounter& element_counts,
                               double desired_probability=0.99999);

Spectrum computeIsotopePatternThr(const ElementCounter& element_counts,
                                  double probability_threshold,
                                  bool relative=false);

Spectrum computeIsotopePattern(const std::string& sum_formula, double desired_probability=0.99999);

Spectrum computeIsotopePatternThr(const std::string& sum_formula,
                                  double probability_threshold=1e-5,
                                  bool relative=false);

namespace mass_search {

struct ElementSettings {
  double monoisotopic_mass;
  size_t min_count, max_count;
  std::string name;

  ElementSettings(const std::string& name,  // must be from the periodic table
      size_t min_count = 0, size_t max_count = 1000)
      : monoisotopic_mass(::ms::monoisotopicMass(name)),
        min_count{min_count},
        max_count{max_count},
        name{name} {}
};

typedef std::vector<size_t> ElementCounts;

struct ExactMassSearchSettings {
  double mass;
  double ppm;
  std::vector<ElementSettings> elements;
  size_t max_results;

  ExactMassSearchSettings(double mass, double ppm,
      const std::vector<ElementSettings>& elements, size_t max_results = 1000000)
      : mass{mass}, ppm{ppm}, elements{elements}, max_results{max_results} {}
};

struct Result {
  double mass;
  ElementCounts counts;
  std::shared_ptr<const ExactMassSearchSettings> settings;

  Result(double, const ElementCounts&,
      const std::shared_ptr<const ExactMassSearchSettings>&);

  Result() : mass{0.0}, counts{}, settings{nullptr} {};

  std::string sumFormula() const;
};

/**
   Class for searching sum formulas that have a given monoisotopic mass.

   Example of usage:

   using ms::mass_search;
   ExactMassSearch search{{512.436, 1.0,
       {{"C", 0, 100}, {"H", 0, 100}, {"O", 0, 100}, {"N", 0, 20}}
   }};

   for (auto& r: search.run())
       std::cout << r.sumFormula() << std::endl;
 **/
struct ExactMassSearch {
  ExactMassSearch(const ExactMassSearchSettings& settings);

  const std::shared_ptr<ExactMassSearchSettings> settings;

  std::vector<Result> run() const;
};

struct ExactMassSearchWithAdduct {
  ExactMassSearchWithAdduct(const ExactMassSearchSettings& settings,
      const std::vector<std::string>& possible_adducts, int charge);

  std::map<std::string, std::vector<Result>> run() const;

  ExactMassSearchSettings settings;

  std::vector<std::string> possible_adducts;
  int charge;
};
}
}

namespace sf_parser {
class ParseError : public std::exception {
  std::string msg_;

 public:
  ParseError(const std::string& str, size_t pos) {
    std::ostringstream ss;
    ss << str << " at position " << pos;
    msg_ = ss.str();
  }

  virtual const char* what() const noexcept { return msg_.c_str(); }
};

class NegativeTotalError : public std::exception {
  std::string msg_;

 public:
  NegativeTotalError(const std::string& element, int total) {
    std::ostringstream ss;
    ss << "total number of " << element << " elements (" << total
       << ") is less than zero";
    msg_ = ss.str();
  }

  virtual const char* what() const noexcept { return msg_.c_str(); }
};

ms::ElementCounter parseSumFormula(const std::string& formula);
}
