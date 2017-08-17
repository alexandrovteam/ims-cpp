#include "ms/isocalc.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <new>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace ms {

struct LogFacTable {
  // just compute them all to avoid branch mispredictions;
  // nobody cares about 0.5MB these days
  double table[65536];

  LogFacTable() {
    for (int i = 0; i < 65536; i++) table[i] = std::lgamma(i + 1);
  }

  double operator[](size_t i) const { return table[i]; }
};

static LogFacTable log_fac_table;

double logFactorial(uint16_t n) {
  return log_fac_table[n];
}

class SingleElementConf {
public:
  // 10 is the maximum number of isotopes known (Sn element)
  using CountVec = std::array<uint16_t, 10>;

private:
  const ms::Element* element_;
  CountVec counts_;
  uint8_t n_isotopes_;

public:
  SingleElementConf(const ms::Element* element, const CountVec& counts={}) :
    element_{element},
    counts_(counts),
    n_isotopes_{static_cast<uint8_t>(element->isotope_pattern.size())}
  {
  }

  bool operator==(const SingleElementConf& other) const {
    assert(element_->abbr == other.element_->abbr);
    return counts_ == other.counts_;
  }

  bool operator<(const SingleElementConf& other) const {
    assert(element_->abbr == other.element_->abbr);
    return counts_ < other.counts_;
  }

  // multinomial distribution
  //
  // according to the paper, we don't try to perform log prob updates on neighboring moves
  // because that would lead to accumulating errors
  double getLogProbability() const {
    double result = 0.0;
    int total_count = 0;

    for (int i = 0; i < size(); i++) {
      uint16_t count = counts_[i];
      double log_prob = element_->logProbability(i);
      result += count * log_prob - logFactorial(count);
      total_count += count;
    }

    return result + logFactorial(total_count);
  }

  double getMass() const {
    double mass = 0.0;
    for (int i = 0; i < size(); i++)
      mass += counts_[i] * element_->isotope_pattern.masses[i];
    return mass;
  }

  // transition into a neighboring isotope configuration
  // by decrementing amount of i-th isotope and incrementing amount of j-th
  void neighborMove(size_t i, size_t j) {
    assert(counts_[i] > 0);
    --counts_[i];
    ++counts_[j];
  }

  uint16_t operator[](size_t idx) const {
    return counts_[idx];
  }

  size_t size() const {
    return n_isotopes_;
  }

  const ms::Element* element() const {
    return element_;
  }

  std::string toString() const {
    std::stringstream ss;
    ss << "SingleElementConf(" << this->element_->abbr << "{";
    for (int i = 0; i < this->size(); i++) {
      ss << counts_[i];
      if (i + 1 < this->size())
        ss << ", ";
    }
    ss << "}" << " | p = " << std::exp(getLogProbability()) << ")";
    return ss.str();
  }
};

}

namespace std {
template <>
struct hash<ms::SingleElementConf> {
  size_t operator()(const ms::SingleElementConf& conf) const {
    size_t seed = 0;
    for (int i = 0; i < conf.size(); i++)
      seed ^= conf[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};
}

namespace ms {

class KahanAdder {
  double sum_;
  double c_;
public:
  KahanAdder(double initial_sum=0.0) :
    sum_(initial_sum), c_(0.0)
  {}

  void add(double x) {
    double y = x - c_;
    double t = sum_ + y;
    c_ = (t - sum_) - y;
    sum_ = t;
  }

  double sum() const {
    return sum_;
  }
};

class SingleElementConfGenerator {
  using CountVec = SingleElementConf::CountVec;

  const ms::Element* element_;
  int atom_count_;

  struct CompareConf {
    const std::vector<SingleElementConf>* confs;
    bool operator()(size_t i, size_t j) {
      assert(i < confs->size());
      assert(j < confs->size());
      return (*confs)[i].getLogProbability() < (*confs)[j].getLogProbability();
    }
  };

  std::vector<SingleElementConf> configurations_;

  CompareConf compare_;

  // queue_ refers to configurations_ element indices
  std::priority_queue<size_t, std::vector<size_t>, CompareConf> queue_;

  // indices of configurations_ ordered by decreasing probability
  std::vector<size_t> order_;

  // configurations_ as a set to avoid duplicates
  std::unordered_set<SingleElementConf> visited_;

  // minimum log probability to consider (nobody cares about 1e-9 and below)
  double log_threshold_;

  void foreachNeighbor(SingleElementConf& conf,
                       const std::function<bool(const SingleElementConf&)>& continue_)
  {
    const int N = conf.size();
    for (int i = 0; i < N; i++) {
      if (conf[i] == 0)
        continue;
      for (int j = 0; j < N; j++) {
        if (i == j)
          continue;
        conf.neighborMove(i, j);
        bool done = !continue_(conf);
        conf.neighborMove(j, i);
        if (done)
          return;
      }
    }
  }

  SingleElementConf makeInitialGuess() {
    // TODO: make counts proportional to probabilities
    CountVec counts{};
    counts[0] = atom_count_;
    return SingleElementConf(element_, counts);
  }

  void pushConfiguration(const SingleElementConf& conf) {
    configurations_.push_back(conf);
    visited_.insert(conf);
    queue_.push(configurations_.size() - 1);
  }

public:
  SingleElementConfGenerator() : element_(nullptr) {}

  /*
    CAVEAT: since compare_ contains a pointer to configurations_,
    default copy constructor is NOT going to work as intended!
    Therefore only empty generators are allowed to be copied
    (which is required for std::vector operations)
  */
  SingleElementConfGenerator(const SingleElementConfGenerator& gen) {
    assert(gen.element_ == nullptr);
  }

  SingleElementConfGenerator(const ms::Element* element,
                             uint16_t atom_count,
                             double log_threshold=-20.0) :
    element_(element),
    atom_count_(atom_count),
    configurations_{},
    compare_{&configurations_},
    queue_{compare_},
    log_threshold_{log_threshold}
  {
    assert(atom_count > 0);
    SingleElementConf conf = makeInitialGuess();
    SingleElementConf next_conf(element);

    double curr_log_prob = conf.getLogProbability(), next_log_prob;
    bool found_more_probable_conf = true;
    while (found_more_probable_conf) {
      found_more_probable_conf = false;

      foreachNeighbor(conf, [&](const SingleElementConf& new_conf) -> bool {
        next_log_prob = new_conf.getLogProbability();
        if (next_log_prob > curr_log_prob) {
          next_conf = new_conf;
          found_more_probable_conf = true;
          return false;
        }
        return true;
      });

      if (found_more_probable_conf) {
        curr_log_prob = next_log_prob;
        conf = next_conf;
      }
    }

    pushConfiguration(conf);
  }

  // i-th most common configuration (0-based index)
  const SingleElementConf& getConfiguration(size_t i) const {
    return configurations_[order_[i]];
  }

  bool advance() {
    if (queue_.empty())
      return false;

    size_t conf_idx = queue_.top();
    queue_.pop();

    auto top_conf = configurations_[conf_idx];
    if (top_conf.getLogProbability() < log_threshold_)
      return false;

    order_.push_back(conf_idx);

    foreachNeighbor(top_conf, [&](const SingleElementConf& new_conf) -> bool {
      if (visited_.count(new_conf) == 0)
        pushConfiguration(new_conf);
      return true;
    });

    return true;
  }

  // number of most probable configurations computed so far
  size_t size() const {
    return order_.size();
  }
};

struct MultiElementConf {
  using SubgeneratorVec = std::vector<SingleElementConfGenerator>;

  const SubgeneratorVec* subgenerators;

private:
  double log_prob_;
  const size_t* sub_conf_indices_;

  void computeLogProbability() {
    log_prob_ = 0.0;
    for (size_t i = 0; i < size(); i++)
      log_prob_ += getSubConfiguration(i).getLogProbability();
  }

public:
  MultiElementConf(const SubgeneratorVec* subgenerators,
                   const size_t* indices) :
    subgenerators(subgenerators),
    sub_conf_indices_(indices)
  {
    computeLogProbability();
  }

  MultiElementConf() :
    subgenerators(nullptr),
    sub_conf_indices_(nullptr)
  {
  }

  const SingleElementConf& getSubConfiguration(size_t i) const {
    return (*subgenerators)[i].getConfiguration(sub_conf_indices_[i]);
  }

  double getLogProbability() const {
    return log_prob_;
  }

  double getMass() const {
    double mass = 0.0;
    for (size_t i = 0; i < size(); i++)
      mass += getSubConfiguration(i).getMass();
    return mass;
  }

  size_t operator[](size_t i) const {
    return sub_conf_indices_[i];
  }

  size_t size() const {
    return subgenerators->size();
  }

  void copyConfIndices(size_t* out) const {
    std::copy(sub_conf_indices_, sub_conf_indices_ + size(), out);
  }

  std::string toString() const {
    std::stringstream ss;
    ss << "MultiElementConf(";
    for (size_t i = 0; i < size(); i++) {
      if (i > 0)
        ss << ", ";
      const auto& sub_conf = getSubConfiguration(i);
      ss << sub_conf.element()->abbr << "{";
      for (size_t j = 0; j < sub_conf.size(); j++) {
        if (j > 0)
          ss << ", ";
        ss << sub_conf[j];
      }
      ss << "}";
    }
    ss << " | p = " << std::exp(getLogProbability()) << ")";
    return ss.str();
  }
};

class MultiElementConfGenerator {
  std::vector<ms::Element> elements_;
  MultiElementConf::SubgeneratorVec subgenerators_;
  KahanAdder total_prob_;

  double expansion_factor_ = 0.25;
  double log_prob_threshold_;

  using ConfigurationVec = std::vector<MultiElementConf>;

  ConfigurationVec current_layer_, next_layer_, accepted_;

  // allocation of many tiny vectors leads to terrible performance,
  // so we allocate configuration memory from large chunks
  std::vector<std::vector<size_t>> data_chunks_;
  static const size_t chunk_size_ = 1024;
  size_t used_in_last_chunk_ = 0;

  size_t* allocateConf() {
    size_t width = elements_.size();

    if (data_chunks_.empty() || used_in_last_chunk_ == chunk_size_) {
      data_chunks_.push_back(std::vector<size_t>(chunk_size_ * width));
      used_in_last_chunk_ = 0;
    }

    return &data_chunks_.back()[used_in_last_chunk_++ * width];
  }

public:
  MultiElementConfGenerator(const ms::ElementCounter& element_counts) {
    subgenerators_.reserve(element_counts.size());
    elements_.reserve(element_counts.size());

    current_layer_.reserve(128);
    next_layer_.reserve(128);
    accepted_.reserve(256);

    for (const auto& elem: element_counts) {
      elements_.push_back(ms::Element::getByName(elem.first));
    }

    int i = 0;
    for (const auto& elem: element_counts) {
      subgenerators_.emplace_back(&elements_[i++], elem.second);
      subgenerators_.back().advance();
    }

    size_t* indices = allocateConf();
    std::fill(indices, indices + elements_.size(), 0);
    current_layer_.emplace_back(&subgenerators_, indices);
    log_prob_threshold_ = current_layer_.back().getLogProbability();
  }

  // number of unique chemical elements
  size_t dim() const {
    return subgenerators_.size();
  }

  // generates the next layer and returns true if it's non-empty
  bool advance() {
    while (!current_layer_.empty()) {
      const auto conf = std::move(current_layer_.back());
      current_layer_.pop_back();

      double log_prob = conf.getLogProbability();
      assert(log_prob >= log_prob_threshold_);

      accepted_.push_back(conf);
      total_prob_.add(std::exp(log_prob));

      for (size_t i = 0; i < dim(); i++) {
        auto& subgen = subgenerators_[i];
        while (subgen.size() <= conf[i] + 1 && subgen.advance());
        if (subgen.size() <= conf[i] + 1)
          continue;

        size_t* indices = allocateConf();
        conf.copyConfIndices(indices);
        ++indices[i];
        MultiElementConf new_conf{&subgenerators_, indices};
        double new_prob = new_conf.getLogProbability();
        if (new_prob >= log_prob_threshold_)
          current_layer_.push_back(new_conf);
        else
          next_layer_.push_back(new_conf);

        // A trick to avoid visiting the same configuration twice:
        // each configuration (k_1, ..., k_n) is reached via the unique path
        // such that rightmost coordinates change first, i.e.
        // (0, ..., 0) -> (0, ..., 1) -> ... -> (0, ..., k_n) -> (0, ..., 1, k_n) ->
        // ... -> (0, ..., k_{n-1}, k_n) -> ... (1, ..., k_n) -> ... -> (k_1, ..., k_n)
        if (conf[i] != 0)
          break;
      }
    }

    if (next_layer_.empty())
      return false;

    current_layer_.swap(next_layer_);

    size_t num_kept = std::floor(current_layer_.size() * expansion_factor_);
    auto middle_it = current_layer_.begin() + num_kept;
    auto log_prob_compare = [](const MultiElementConf& c1, const MultiElementConf& c2) {
      return c1.getLogProbability() > c2.getLogProbability();
    };
    std::nth_element(current_layer_.begin(), middle_it, current_layer_.end(), log_prob_compare);
    log_prob_threshold_ = middle_it->getLogProbability();

    next_layer_.assign(middle_it + 1, current_layer_.end()); // all < new threshold
    current_layer_.resize(num_kept + 1); // all >= new threshold

    return true;
  }


  // unsorted top configurations whose probability adds up to collectedProbability()
  const std::vector<MultiElementConf>& configurations() const {
    return accepted_;
  }

  double collectedProbability() const {
    return total_prob_.sum();
  }
};

double monoisotopicMass(const ElementCounter& counter) {
  double sum = 0.0;
  for (auto& item : counter)
    sum += item.second * ms::Element::getByName(item.first).isotope_pattern.masses[0];
  return sum;
}

double monoisotopicMass(const std::string& formula) {
  return monoisotopicMass(::sf_parser::parseSumFormula(formula));
}

Spectrum computeIsotopePattern(const ElementCounter& element_counts, double desired_probability)
{
  assert(0 < desired_probability && desired_probability < 1);

  MultiElementConfGenerator gen{element_counts};
  while (gen.collectedProbability() < desired_probability && gen.advance());

  std::vector<double> masses, intensities;
  for (const auto& conf: gen.configurations()) {
    masses.push_back(conf.getMass());
    intensities.push_back(std::exp(conf.getLogProbability()));
  }

  ms::Spectrum isotope_pattern{masses, intensities};
  isotope_pattern.sortByIntensity();
  KahanAdder totalProb;
  size_t i = 0;
  while (totalProb.sum() < desired_probability)
    totalProb.add(isotope_pattern.intensities[i++]);
  isotope_pattern.masses.resize(i);
  isotope_pattern.intensities.resize(i);
  isotope_pattern.normalize();
  return isotope_pattern;
}

Spectrum computeIsotopePattern(
    const std::string& sum_formula, double desired_probability) {
  auto counts = sf_parser::parseSumFormula(sum_formula);
  return computeIsotopePattern(counts, desired_probability);
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

  const ElementCounter& elementCounts() const { return counter_; }

 private:
  bool eof() const { return n >= s.length(); }

  void checkAvailable() const {
    if (eof()) throw ParseError("unexpected end of input", n);
  }

  char nextChar() {
    checkAvailable();
    return s[n++];
  }

  char peek() const {
    checkAvailable();
    return s[n];
  }

  uint16_t parseOptionalNumber() {
    if (eof() || !std::isdigit(peek())) return 1;
    auto pos = n;
    uint16_t result = 0;
    while (!eof() && std::isdigit(peek())) {
      if (result >= 3000) throw ParseError("the number is too big", pos);
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
      for (auto& item : tmp)
        counter[item.first] += item.second * repeats;
    } else if (std::isupper(peek())) { parseSimpleFragment(counter); } else {
      throw ParseError(std::string{"unexpected character '"} + peek() + "'", n);
    }
  }

  void parseMolecularComplex(ElementCounter& counter) {
    ElementCounter tmp;
    auto repeats = parseOptionalNumber();
    parseFragment(tmp);
    while (!eof()) {
      if (peek() == '.' || peek() == '-' || peek() == '+') break;
      parseFragment(tmp);
    }
    for (auto& item : tmp)
      counter[item.first] += repeats * item.second;
  }

  void parseSumFormula(ElementCounter& counter) {
    parseMolecularComplex(counter);

    while (!eof()) {
      if (peek() == '.') {
        nextChar();
        parseMolecularComplex(counter);
      } else { break; }
    }

    while (!eof()) {
      ElementCounter adduct;
      char sign = nextChar();
      int mult;
      if (sign == '-')
        mult = -1;
      else if (sign == '+')
        mult = 1;
      else
        throw ParseError("expected +/-", n - 1);
      parseMolecularComplex(adduct);
      for (auto& item : adduct)
        counter[item.first] += mult * int(item.second);
    }

    for (auto it = counter.begin(); it != counter.end(); ) {
      if (it->second < 0) throw NegativeTotalError(it->first, it->second);
      if (it->second != 0)
        ++it;
      else
        it = counter.erase(it);
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
      if (counts[i] > 1) ss << counts[i];
    }
  }
  return ss.str();
}

Result::Result(double mass, const ElementCounts& counts,
    const std::shared_ptr<const ExactMassSearchSettings>& settings)
    : mass{mass}, counts{counts}, settings{settings} {
}

/*
  value of N > 0 means that N'th element count is explored;
  value of N = 0 means checking and storing current result
*/
template <size_t N>
struct FormulaGenerator {
  static inline void run(std::vector<Result>& results,
      const std::shared_ptr<const ExactMassSearchSettings>& settings,
      ElementCounts& current_counts, double current_mass = 0.0) {
    const ElementSettings& element = settings->elements[N - 1];
    double element_mass = element.monoisotopic_mass;
    double max_mass = settings->mass * (1.0 + settings->ppm * 1e-6);
    size_t min_count = element.min_count;
    size_t max_count =
        std::min(element.max_count, size_t((max_mass - current_mass) / element_mass));

    auto old_value = current_counts[N - 1];

    for (size_t count = min_count; count <= max_count; count++) {
      double current_mass_i = current_mass + element_mass * count;
      current_counts[N - 1] = count;
      FormulaGenerator<N - 1>::run(results, settings, current_counts, current_mass_i);
    }

    current_counts[N - 1] = old_value;
  }
};

template <>
struct FormulaGenerator<0> {
  static inline void run(std::vector<Result>& results,
      const std::shared_ptr<const ExactMassSearchSettings>& settings,
      ElementCounts& current_counts, double current_mass) {
    if (results.size() >= settings->max_results) return;

    auto ppm = (current_mass - settings->mass) / settings->mass * 1e6;
    if (std::fabs(ppm) > settings->ppm) return;

    results.emplace_back(current_mass, current_counts, settings);
  }
};

ExactMassSearch::ExactMassSearch(const ExactMassSearchSettings& settings)
    : settings(std::make_shared<ExactMassSearchSettings>(settings)) {
}

std::vector<Result> ExactMassSearch::run() const {
  std::vector<Result> results;
  ElementCounts counts(this->settings->elements.size());
  auto p = settings;

  switch (this->settings->elements.size()) {
    case 0:
      break;
    case 1:
      FormulaGenerator<1>::run(results, p, counts);
      break;
    case 2:
      FormulaGenerator<2>::run(results, p, counts);
      break;
    case 3:
      FormulaGenerator<3>::run(results, p, counts);
      break;
    case 4:
      FormulaGenerator<4>::run(results, p, counts);
      break;
    case 5:
      FormulaGenerator<5>::run(results, p, counts);
      break;
    case 6:
      FormulaGenerator<6>::run(results, p, counts);
      break;
    case 7:
      FormulaGenerator<7>::run(results, p, counts);
      break;
    case 8:
      FormulaGenerator<8>::run(results, p, counts);
      break;
    case 9:
      FormulaGenerator<9>::run(results, p, counts);
      break;
    case 10:
      FormulaGenerator<10>::run(results, p, counts);
      break;
    default:
      throw std::runtime_error("at most 10 distinct elements are supported");
  }

  return results;
}

ExactMassSearchWithAdduct::ExactMassSearchWithAdduct(
    const ExactMassSearchSettings& settings,
    const std::vector<std::string>& possible_adducts, int charge)
    : settings{settings}, possible_adducts{possible_adducts}, charge{charge} {
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
    { results[adduct] = r; }
  }
  return results;
}
}
}

int main() {
  for (int i = 0; i < 1000; i++) {
    ms::computeIsotopePattern("C100H200O50P5N10", 0.99999);
  }
}
