#include "ms/isocalc.hpp"

#include "fftw3.h"
#include <cassert>
#include <cstring>
#include <complex>
#include <vector>
#include <mutex>
#include <new>
#include <stdexcept>
#include <sstream>

namespace ms {
  namespace detail {
    __attribute ((__destructor__)) void fftwCleanup() { fftw_cleanup(); }

    class FftwArray {
      std::vector<int> dims_;
      std::complex<double>* data_;
      size_t n_;
      public:
      FftwArray(const std::vector<int>& dimensions) : dims_(dimensions) {
        n_ = 1;
        for (int x: dimensions)
          n_ *= x;
        if (n_ > 1e7)
          throw std::runtime_error("too many isotopic combinations");
        data_ = reinterpret_cast<std::complex<double>*>(fftw_alloc_complex(n_));
        if (data_ == nullptr)
          throw std::bad_alloc();
        std::memset(data_, 0, n_ * sizeof(fftw_complex));
      }

      fftw_complex* fftw_data() const { return reinterpret_cast<fftw_complex*>(data_); }
      std::complex<double>* data() const { return data_; }
      size_t size() const { return n_; }

      ~FftwArray() {
        fftw_free(data_);
      }
    };
  }

  double monoisotopicMass(const ElementCounter& counter) {
    double sum = 0.0;
    for (auto& item: counter)
      sum += item.second * ms::Element::getByName(item.first).isotope_pattern.masses[0];
    return sum;
  }

  double monoisotopicMass(const std::string& formula) {
    return monoisotopicMass(::sf_parser::parseSumFormula(formula));
  }

  IsotopePattern computeIsotopePattern(const ms::Element& element, size_t amount, double threshold) {
    size_t dim = element.isotope_pattern.size() - 1;
    std::vector<int> dimensions(static_cast<int>(dim), amount + 1);
    detail::FftwArray arr{dimensions};

    static std::mutex plan_mutex; // planning is not thread-safe

    // array setup
    const auto& iso = element.isotope_pattern;
    arr.data()[0] = iso.abundances[0];
    for (size_t i = 0, k = 1; i < dim; i++, k *= amount + 1)
      arr.data()[k] = iso.abundances[dim - i];

    plan_mutex.lock();
    // forward FFT
    auto fwd_plan = fftw_plan_dft(int(dim), &dimensions[0], arr.fftw_data(), arr.fftw_data(),
                                 FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(fwd_plan);

    // exponentiation
    for (size_t i = 0; i < arr.size(); ++i)
      arr.data()[i] = std::pow(arr.data()[i], amount);
    fftw_destroy_plan(fwd_plan);

    // inverse FFT
    auto bwd_plan = fftw_plan_dft(dim, &dimensions[0], arr.fftw_data(), arr.fftw_data(),
                                  FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(bwd_plan);
    for (size_t i = 0; i < arr.size(); ++i)
      arr.data()[i] /= double(arr.size()); // take care of FFTW normalization
    fftw_destroy_plan(bwd_plan);
    plan_mutex.unlock();

    ms::IsotopePattern isotope_pattern;

    std::vector<size_t> indices(dim, 0);
    for (size_t i = 0; i < arr.size(); ++i) {
      size_t k;
      size_t n = 0;
      double mass = 0.0;
      for (k = 0; k < dim; ++k) {
        mass += iso.masses[k + 1] * indices[k];
        n += indices[k];
      }

      if (n <= amount && arr.data()[i].real() >= threshold) {
        mass += (amount - n) * iso.masses[0];

        isotope_pattern.masses.push_back(mass);
        isotope_pattern.abundances.push_back(arr.data()[i].real());
      }

      if (i == arr.size() - 1) break;
      for (k = dim - 1; indices[k] == amount; --k);
      indices[k] += 1;
      while (++k < dim) indices[k] = 0;
    }

    isotope_pattern.normalize();
    return isotope_pattern;
  }

  IsotopePattern computeIsotopePattern(const ElementCounter& element_counts,
      double threshold, double fft_threshold)
  {
    assert(threshold < 1);
    ms::IsotopePattern result(0);

    for (auto& item: element_counts) {
      assert(ms::Element::isKnown(item.first));
      auto element = ms::Element::getByName(item.first);
      try {
        auto pattern = computeIsotopePattern(element, item.second, fft_threshold);
        result = result.multiply(pattern, threshold * threshold);
      } catch (std::runtime_error& e) {
        std::stringstream full_desc;
        full_desc << "failed to compute isotope pattern for "
                  << item.first << item.second << ": " << e.what();
        throw std::runtime_error(full_desc.str());
      }
    }

    return result.removeAbundancesBelow(threshold);
  }

  IsotopePattern computeIsotopePattern(const std::string& sum_formula,
      double threshold, double fft_threshold)
  {
    auto counts = sf_parser::parseSumFormula(sum_formula);
    return computeIsotopePattern(counts, threshold, fft_threshold);
  }
}

namespace sf_parser {
  typedef ms::ElementCounter ElementCounter;

  using ::sf_parser::ParseError;
  using ::sf_parser::NegativeTotalError;

  class SumFormulaParser {
    std::string s;
    size_t n;
    ElementCounter counter_;
  public:
    SumFormulaParser(const std::string& input) : s(input), n(0) {
      parseSumFormula(counter_);
    }

    const ElementCounter& elementCounts() const {
      return counter_;
    }

  private:
    bool eof() const { return n >= s.length(); }

    void checkAvailable() const {
      if (eof()) throw ParseError("unexpected end of input", n);
    }

    char nextChar() { checkAvailable(); return s[n++]; }

    char peek() const { checkAvailable(); return s[n]; }

    uint16_t parseOptionalNumber() {
      if (eof() || !std::isdigit(peek()))
        return 1;
      auto pos = n;
      uint16_t result = 0;
      while (!eof() && std::isdigit(peek())) {
        if (result >= 6553)
          throw ParseError("the number is too big", pos);
        result = result * 10 + (nextChar() - '0');
      }
      return result;
    }

    void parseElement(ElementCounter& counter) {
      std::string element;
      element.push_back(nextChar());
      auto pos = n;
      if (!std::isupper(element.back())) throw ParseError("expected an element", n);
      while (!eof() && std::islower(peek()))
        element.push_back(nextChar());
      uint16_t num = parseOptionalNumber();
      if (!ms::Element::isKnown(element))
        throw ParseError("unknown element " + element, pos);
      counter[element] += num;
    }

    void parseSimpleFragment(ElementCounter& counter) {
      while (!eof() && std::isupper(peek())) {
        parseElement(counter);
      }
    }

    void parseFragment(ElementCounter& counter) {
      checkAvailable();

      if (peek() == '(') {
        nextChar();
        ElementCounter tmp;
        parseFragment(tmp);
        if (eof() || nextChar() != ')')
          throw ParseError("expected closing parenthesis", n - 1);
        auto repeats = parseOptionalNumber();
        for (auto& item: tmp)
          counter[item.first] += item.second * repeats;
      } else if (std::isupper(peek())) {
        parseSimpleFragment(counter);
      } else {
        throw ParseError(std::string{"unexpected character '"} + peek() + "'", n);
      }
    }

    void parseMolecularComplex(ElementCounter& counter) {
      ElementCounter tmp;
      auto repeats = parseOptionalNumber();
      while (!eof()) {
        if (peek() == '.' || peek() == '-' || peek() == '+')
          break;
        parseFragment(tmp);
      }
      for (auto& item: tmp)
        counter[item.first] += repeats * item.second;
    }

    void parseSumFormula(ElementCounter& counter) {
      parseMolecularComplex(counter);

      while (!eof()) {
        if (peek() == '.') {
          nextChar();
          parseMolecularComplex(counter);
        } else {
          break;
        }
      }

      // FIXME: support multiple adducts
      if (!eof()) {
        ElementCounter adduct;
        char sign = nextChar();
        int mult;
        if (sign == '-') mult = -1;
        else if (sign == '+') mult = 1;
        else throw ParseError("expected +/-", n-1);
        parseMolecularComplex(adduct);
        for (auto& item: adduct) {
          auto total = int(counter[item.first]) + mult * int(item.second);
          if (total < 0)
            throw NegativeTotalError(item.first, total);
          if (total == 0)
            counter.erase(item.first);
          else
            counter[item.first] = total;
        }
      }
    }
  };

  ms::ElementCounter parseSumFormula(const std::string& formula) {
    sf_parser::SumFormulaParser parser{formula};
    return parser.elementCounts();
  }
}

namespace ms {
  namespace mass_search {

    std::string Result::sumFormula() const {
      std::stringstream ss;
      for (size_t i = 0; i < settings->elements.size(); i++) {
        if (counts[i] > 0) {
          ss << settings->elements[i].name;
          if (counts[i] > 1)
            ss << counts[i];
        }
      }
      return ss.str();
    }

    Result::Result(double mass, const ElementCounts& counts,
                   const ExactMassSearch* const settings) :
      mass{mass}, counts{counts}, settings{settings}
    {
    }

    /*
      value of N > 0 means that N'th element count is explored;
      value of N = 0 means checking and storing current result
    */
    template<size_t N>
    struct FormulaGenerator {
      static inline void run(std::vector<Result>& results,
                             const ExactMassSearch* const settings,
                             ElementCounts& current_counts,
                             double current_mass = 0.0)
      {
        const ElementSettings& element = settings->elements[N - 1];
        double element_mass = element.monoisotopic_mass;
        double max_mass = settings->mass * (1.0 + settings->ppm * 1e-6);
        size_t min_count = element.min_count;
        size_t max_count = std::min(element.max_count,
                                    size_t((max_mass - current_mass) / element_mass));

        auto old_value = current_counts[N - 1];

        for (size_t count = min_count; count <= max_count; count++) {
          double current_mass_i = current_mass + element_mass * count;
          current_counts[N - 1] = count;
          FormulaGenerator<N - 1>::run(results, settings, current_counts, current_mass_i);
        }

        current_counts[N - 1] = old_value;
      }
    };

    template<>
    struct FormulaGenerator<0> {
      static inline void run(std::vector<Result>& results,
                             const ExactMassSearch* const settings,
                             ElementCounts& current_counts,
                             double current_mass)
      {
        if (results.size() >= settings->max_results)
          return;

        auto ppm = (current_mass - settings->mass) / settings->mass * 1e6;
        if (std::fabs(ppm) > settings->ppm)
          return;

        results.emplace_back(current_mass, current_counts, settings);
      }
    };

    ExactMassSearch::ExactMassSearch(double mass, double ppm,
                                     const std::vector<ElementSettings> elements,
                                     size_t max_results) :
      mass{mass}, ppm{ppm}, elements{elements}, max_results{max_results}
    {
    }

    std::vector<Result> ExactMassSearch::run() const {
      std::vector<Result> results;
      ElementCounts counts(this->elements.size());

      switch (this->elements.size()) {
      case 0:  break;
      case 1:  FormulaGenerator< 1>::run(results, this, counts); break;
      case 2:  FormulaGenerator< 2>::run(results, this, counts); break;
      case 3:  FormulaGenerator< 3>::run(results, this, counts); break;
      case 4:  FormulaGenerator< 4>::run(results, this, counts); break;
      case 5:  FormulaGenerator< 5>::run(results, this, counts); break;
      case 6:  FormulaGenerator< 6>::run(results, this, counts); break;
      case 7:  FormulaGenerator< 7>::run(results, this, counts); break;
      case 8:  FormulaGenerator< 8>::run(results, this, counts); break;
      case 9:  FormulaGenerator< 9>::run(results, this, counts); break;
      case 10: FormulaGenerator<10>::run(results, this, counts); break;
      default:
        throw std::runtime_error("at most 10 distinct elements are supported");
      }

      return results;
    }
  }
}
