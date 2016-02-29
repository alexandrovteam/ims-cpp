#include "ms/isocalc.hpp"

#include "fftw3.h"
#include <cassert>
#include <cstring>
#include <complex>
#include <vector>
#include <mutex>
#include <new>

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
      auto pattern = computeIsotopePattern(element, item.second, fft_threshold);
      result = result.multiply(pattern, threshold * threshold);
    }

    return result;
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
