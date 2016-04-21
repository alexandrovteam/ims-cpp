#include "ms/isocalc.hpp"

#include "third-party/kissfft/kiss_fft.h"
#include "third-party/kissfft/tools/kiss_fftndr.h"

static_assert(sizeof(kiss_fft_scalar) == 8,
              "use -Dkiss_fft_scalar for double precision");

#include <cassert>
#include <cstring>
#include <complex>
#include <vector>
#include <mutex>
#include <new>
#include <stdexcept>
#include <sstream>
#include <memory>

namespace ms {
  namespace detail {

    class KissFftState {
      std::vector<int> dims_;
      bool inverse_;
      void* cfg_;
      bool is1d() const { return dims_.size() == 1; }
    public:
      KissFftState(const std::vector<int>& dimensions, bool inverse=false) :
        dims_(dimensions), inverse_(inverse)
      {
        if (is1d())
          cfg_ = kiss_fftr_alloc(dims_[0], inverse, nullptr, nullptr);
        else
          cfg_ = kiss_fftndr_alloc(&dims_[0], dims_.size(), inverse, nullptr, nullptr);
        if (cfg_ == nullptr)
          throw std::bad_alloc();
      }

      KissFftState(const KissFftState&) = delete;
      KissFftState operator=(const KissFftState&) = delete;
      ~KissFftState() {
        free(cfg_);
      }

      void runFFT(void* in, void* out) {
        if (is1d()) { // for some reason kiss_fftndr needs at least 2 dimensions to work
          if (!inverse_)
            kiss_fftr((kiss_fftr_cfg)cfg_, (kiss_fft_scalar*)in, (kiss_fft_cpx*)out);
          else
            kiss_fftri((kiss_fftr_cfg)cfg_, (kiss_fft_cpx*)in, (kiss_fft_scalar*)out);
        } else {
          if (!inverse_)
            kiss_fftndr((kiss_fftndr_cfg)cfg_, (kiss_fft_scalar*)in, (kiss_fft_cpx*)out);
          else
            kiss_fftndri((kiss_fftndr_cfg)cfg_, (kiss_fft_cpx*)in, (kiss_fft_scalar*)out);
        }
      }
    };

    static_assert(sizeof(kiss_fft_cpx) == sizeof(std::complex<kiss_fft_scalar>), "");

    class FftArray {
      std::vector<int> dims_;
      std::complex<kiss_fft_scalar>* data_;
      size_t n_;
    public:
      FftArray(const std::vector<int>& dimensions) :
        dims_(dimensions), data_(nullptr), n_(1)
      {
        for (int x: dims_) n_ *= x;

        if (n_ > 1e7)
          throw std::runtime_error("too many isotopic combinations");

        size_t n_bytes = n_ * sizeof(kiss_fft_cpx);

        data_ = reinterpret_cast<std::complex<kiss_fft_scalar>*>(KISS_FFT_MALLOC(n_bytes));
        if (data_ == nullptr)
          throw std::bad_alloc();
        std::memset(data_, 0, n_bytes);
      }

      void forwardFFT() {
        KissFftState state(dims_, false);
        state.runFFT(scalar_data(), complex_data());
      }

      void inverseFFT() {
        KissFftState state(dims_, true);
        state.runFFT(complex_data(), scalar_data());
      }

      kiss_fft_cpx* complex_data() const {
        return reinterpret_cast<kiss_fft_cpx*>(data_);
      }

      kiss_fft_scalar* scalar_data() const {
        return reinterpret_cast<kiss_fft_scalar*>(data_);
      }

      std::complex<kiss_fft_scalar>* data() const { return data_; }
      size_t size() const { return n_; }

      ~FftArray() {
        kiss_fftr_free(data_);
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
    if (dim == 0)
      return ms::IsotopePattern(element.isotope_pattern.masses[0] * amount);

    size_t edge_len = kiss_fftr_next_fast_size_real(amount + 1);

    std::vector<int> dimensions(int(dim), edge_len);
    detail::FftArray arr{dimensions};

    // array setup
    const auto& iso = element.isotope_pattern;
    arr.scalar_data()[0] = iso.abundances[0];
    for (size_t i = 0, k = 1; i < dim; i++, k *= edge_len)
      arr.scalar_data()[k] = iso.abundances[dim - i];

    // forward FFT
    arr.forwardFFT();

    // exponentiation
    for (size_t i = 0; i < arr.size(); ++i)
      arr.data()[i] = std::pow(arr.data()[i], amount);

    // inverse FFT
    arr.inverseFFT();
    for (size_t i = 0; i < arr.size(); ++i)
      arr.scalar_data()[i] /= double(arr.size()); // take care of FFT normalization

    ms::IsotopePattern isotope_pattern;

    std::vector<size_t> indices(dim, 0);
    for (size_t i = 0; i < arr.size(); ++i) {
      size_t k;
      size_t n = 0;
      double mass = 0.0;
      for (k = 0; k < indices.size(); ++k) {
        mass += iso.masses[k + 1] * indices[k];
        n += indices[k];
      }

      if (n <= amount && arr.scalar_data()[i] >= threshold) {
        mass += (amount - n) * iso.masses[0];

        isotope_pattern.masses.push_back(mass);
        isotope_pattern.abundances.push_back(arr.scalar_data()[i]);
      }

      if (i == arr.size() - 1) break;

      for (k = indices.size() - 1; indices[k] == edge_len - 1; --k);
      indices[k] += 1;
      while (++k < indices.size()) indices[k] = 0;
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
                   const std::shared_ptr<const ExactMassSearchSettings>& settings) :
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
                             const std::shared_ptr<const ExactMassSearchSettings>& settings,
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
                             const std::shared_ptr<const ExactMassSearchSettings>& settings,
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

    ExactMassSearch::ExactMassSearch(const ExactMassSearchSettings& settings) :
      settings(std::make_shared<ExactMassSearchSettings>(settings))
    {
    }

    std::vector<Result> ExactMassSearch::run() const {
      std::vector<Result> results;
      ElementCounts counts(this->settings->elements.size());
      auto p = settings;

      switch (this->settings->elements.size()) {
      case 0:  break;
      case 1:  FormulaGenerator< 1>::run(results, p, counts); break;
      case 2:  FormulaGenerator< 2>::run(results, p, counts); break;
      case 3:  FormulaGenerator< 3>::run(results, p, counts); break;
      case 4:  FormulaGenerator< 4>::run(results, p, counts); break;
      case 5:  FormulaGenerator< 5>::run(results, p, counts); break;
      case 6:  FormulaGenerator< 6>::run(results, p, counts); break;
      case 7:  FormulaGenerator< 7>::run(results, p, counts); break;
      case 8:  FormulaGenerator< 8>::run(results, p, counts); break;
      case 9:  FormulaGenerator< 9>::run(results, p, counts); break;
      case 10: FormulaGenerator<10>::run(results, p, counts); break;
      default:
        throw std::runtime_error("at most 10 distinct elements are supported");
      }

      return results;
    }

    ExactMassSearchWithAdduct::ExactMassSearchWithAdduct(const ExactMassSearchSettings& settings,
          const std::vector<std::string>& possible_adducts, int charge) :
      settings{settings}, possible_adducts{possible_adducts}, charge{charge}
    {
    }

    std::map<std::string, std::vector<Result>> ExactMassSearchWithAdduct::run() const {
      std::map<std::string, std::vector<Result>> results;
#pragma omp parallel for
      for (size_t i = 0; i < possible_adducts.size(); i++) {
        auto adduct = possible_adducts[i];
        ExactMassSearchSettings s = settings;
        s.mass -= -charge * ms::electronMass + ms::monoisotopicMass(adduct);
        ExactMassSearch search(s);
        auto r = search.run();
#pragma omp critical
        {
          results[adduct] = r;
        }
      }
      return results;
    }
  }
}
