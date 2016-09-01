#include "ms/isocalc.hpp"
#include "utils/metrics.hpp"

#include "cxxopts.hpp"

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <functional>
#include <random>
#include <ctime>
#include <iomanip>

typedef utils::Metrics Metrics;

std::vector<Metrics> sortBy(
    const std::vector<Metrics>& metrics, std::function<double(const Metrics&)> f) {
  std::vector<double> values(metrics.size());
  std::vector<size_t> indices(values.size());
  for (size_t i = 0; i < metrics.size(); i++) {
    values[i] = f(metrics[i]);
    indices[i] = i;
  }
  std::sort(indices.begin(), indices.end(),
      [&](size_t n, size_t m) { return values[n] < values[m]; });

  std::vector<Metrics> result(metrics.size());
  size_t k = 0;
  for (auto i : indices)
    result[k++] = metrics[i];
  return result;
}

std::vector<Metrics> readResults(const std::string& filename) {
  std::vector<Metrics> metrics;
  std::ifstream in(filename);
  Metrics m;
  std::string line;

  std::getline(in, line);  // skip header
  if (line != Metrics::header())
    throw std::runtime_error("CSV header is incompatible with ims-fdr tool");

  while (Metrics::read(in, m))
    metrics.push_back(m);

  return sortBy(metrics, [](const Metrics& a) { return -a.msm(); });
}

static std::mt19937 rng;

std::vector<float> estimateFDR(const std::vector<Metrics>& target_metrics,
    const std::vector<Metrics>& decoy_metrics,
    size_t n_top = std::numeric_limits<size_t>::max(), bool verbose = false) {
  std::vector<float> fdr;
  size_t n = target_metrics.size();

  std::vector<size_t> decoy_indices(n);
  std::iota(decoy_indices.begin(), decoy_indices.end(), 0UL);
  // reservoir sampling
  for (size_t i = n; i < decoy_metrics.size(); i++) {
    std::uniform_int_distribution<size_t> gen(0, i);
    auto k = gen(rng);
    if (k < n) decoy_indices[k] = i;
  }

  std::sort(decoy_indices.begin(), decoy_indices.end(), [&](size_t n, size_t m) {
    return decoy_metrics[n].msm() > decoy_metrics[m].msm();
  });

  size_t n_decoy_hits = 0;
  size_t n_target_hits = 0;
  auto it = decoy_indices.begin();
  for (auto& m : target_metrics) {
    if (n_target_hits >= n_top) break;
    while (it != decoy_indices.end() && decoy_metrics[*it].msm() > m.msm()) {
      if (verbose)
        std::cout << "\x1b[31mDECOY  \x1b[39m" << decoy_metrics[*it] << std::endl;
      ++it;
      ++n_decoy_hits;
    }
    ++n_target_hits;
    fdr.push_back(float(n_decoy_hits) / float(n_target_hits));

    if (verbose)
      std::cout << "\x1b[32mTARGET \x1b[39m" << m << "," << fdr.back() << std::endl;
  }
  return fdr;
}

std::vector<float> estimateAverageFDR(const std::vector<Metrics>& target_metrics,
    const std::vector<Metrics>& decoy_metrics, size_t n_repetitions) {
  // FIXME: use median instead of mean
  std::vector<float> average_fdr(target_metrics.size());
  for (size_t j = 0; j < n_repetitions; j++) {
    auto fdr = estimateFDR(target_metrics, decoy_metrics);
    for (size_t i = 0; i < target_metrics.size(); ++i)
      average_fdr[i] += fdr[i] / float(n_repetitions);
  }
  return average_fdr;
}

int fdr_main(int argc, char** argv) {
  std::string target_csv_fn, decoy_csv_fn, adduct, output_filename;
  unsigned n_repeats;
  std::string groundtruth_filename;

  cxxopts::Options options("ims fdr", " <target_results.csv> <decoy_results.csv>");
  options.add_options()("out", "Output filename",
      cxxopts::value<std::string>(output_filename)->default_value("/dev/stdout"))(
      "repeats", "Number of random subsampling runs used for FDR estimation",
      cxxopts::value<unsigned>(n_repeats)->default_value("100"))("adduct",
      "Only use target results for a specific adduct",
      cxxopts::value<std::string>(adduct)->default_value(""))("groundtruth",
      "CSV file with ground-truth formula/adduct pairs",
      cxxopts::value<std::string>(groundtruth_filename)->default_value(""))(
      "help", "Print help");

  options.add_options("hidden")(
      "target_csv", "", cxxopts::value<std::string>(target_csv_fn))(
      "decoy_csv", "", cxxopts::value<std::string>(decoy_csv_fn))("debug", "");

  options.parse_positional(std::vector<std::string>{"target_csv", "decoy_csv"});

  options.parse(argc, argv);

  if (options.count("help") || target_csv_fn.empty() || decoy_csv_fn.empty()) {
    std::cout << options.help({""}) << std::endl;
    return 0;
  }

  auto target_metrics = readResults(target_csv_fn);
  auto decoy_metrics = readResults(decoy_csv_fn);

  if (!adduct.empty()) {
    auto it = std::remove_if(target_metrics.begin(), target_metrics.end(),
        [&](const Metrics& m) { return m.adduct != adduct; });
    target_metrics.erase(it, target_metrics.end());
  }

  rng.seed(std::random_device{}());

  if (options.count("debug")) {
    estimateFDR(target_metrics, decoy_metrics, 25, /* verbose= */ true);
    return 0;
  }

  auto fdr = estimateAverageFDR(target_metrics, decoy_metrics, n_repeats);
  std::ofstream out(output_filename);

  if (groundtruth_filename.empty()) {
    out << Metrics::header() << ",fdr" << std::endl;
    for (size_t j = 0; j < target_metrics.size(); j++)
      out << target_metrics[j] << "," << fdr[j] << std::endl;
  } else {
    std::set<std::pair<ms::ElementCounter, std::string>> groundtruth;

    auto is_correct = [&](const Metrics& metrics) -> bool {
      auto key = std::make_pair(sf_parser::parseSumFormula(metrics.sf), metrics.adduct);
      return groundtruth.find(key) != groundtruth.end();
    };

    std::ifstream gt(groundtruth_filename);
    std::string line;
    while (std::getline(gt, line)) {
      std::istringstream is(line);
      std::string sf, adduct;
      std::getline(is, sf, ',');
      std::getline(is, adduct, ',');
      groundtruth.insert(std::make_pair(sf_parser::parseSumFormula(sf), adduct));
    }

    out << Metrics::header() << ",est_fdr,true_fdr,correct" << std::endl;
    size_t true_hits = 0;

    for (size_t j = 0; j < target_metrics.size(); j++) {
      bool correct = is_correct(target_metrics[j]);
      true_hits += correct;
      double true_fdr = double(j + 1 - true_hits) / double(j + 1);
      out << target_metrics[j] << "," << fdr[j] << "," << true_fdr << "," << correct
          << std::endl;
    }
  }

  return 0;
}
