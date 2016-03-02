#include "utils/isotope_pattern_db.hpp"
#include "ms/isocalc.hpp"
#include "ms/periodic_table.hpp"
#include "cxxopts.hpp"

#include <iostream>
#include <regex>
#include <fstream>
#include <set>
#include <vector>
#include <string>

std::vector<std::string> decoy_adducts;

int isocalc_main(int argc, char** argv) {

  double resolution;
  unsigned max_peaks;

  std::string input_file, output_file;

  std::string adducts_str;

  cxxopts::Options options("ims isocalc",
  " <input.txt> <output.db>\n\t\t\twhere input contains one sum formula per line.");
  options.add_options()
    ("resolution",      "Resolving power at m/z=200",
     cxxopts::value<double>(resolution)->default_value("140000.0"))
    ("adducts",         "Comma-separated list of adducts",
     cxxopts::value<std::string>(adducts_str)->default_value("+H,+K,+Na"))
    ("max-peaks",       "Maximum number of peaks to store",
     cxxopts::value<unsigned>(max_peaks)->default_value("5"))
    ("generate-decoy",  "Generate a decoy database")
    ("help", "Print help");

  options.add_options("hidden")
    ("input-file",      "List of sum formulas, one per line",
     cxxopts::value<std::string>(input_file))
    ("output-file",     "Output file",
     cxxopts::value<std::string>(output_file));

  options.parse_positional(std::vector<std::string>{"input-file", "output-file"});
  options.parse(argc, argv);

  if (options.count("help") || input_file.empty() || output_file.empty()) {
    std::cout << options.help({""}) << std::endl;
    return 0;
  }

  std::regex regex{","};
  std::set<std::string> adducts_set{
    std::sregex_token_iterator(adducts_str.begin(), adducts_str.end(), regex, -1),
    std::sregex_token_iterator()};
  std::set<std::string> tmp;
  for (auto& adduct: adducts_set)
    if (adduct[0] != '+' && adduct[0] != '-')
      tmp.insert("+" + adduct);
    else tmp.insert(adduct);
  adducts_set.swap(tmp);
  tmp.clear();
  bool positive_mode = std::all_of(adducts_set.begin(), adducts_set.end(),
                                   [](const std::string& s) { return s[0] == '+'; });
  bool negative_mode = std::all_of(adducts_set.begin(), adducts_set.end(),
                                   [](const std::string& s) { return s[0] == '-'; });
  if (!positive_mode && !negative_mode) {
    std::cerr << "all adducts must have the same sign" << std::endl;
    return -2;
  }
  std::vector<std::string> target_adducts{adducts_set.begin(), adducts_set.end()};

  std::cout << "Resolution @ m/z=200: " << resolution << std::endl;
  std::cout << "Generating isotope pattern database for adducts ";
  for (auto& a: target_adducts) {
    std::cout << a;
    if (a != target_adducts.back()) std::cout << ", ";
  }
  std::cout << "..." << std::endl;

  std::ifstream in(input_file);

  std::vector<std::string> sum_formulas;

  std::string f;
  while (!in.eof()) {
    std::getline(in, f);
    if (ms::monoisotopicMass(f) > 100 && ms::monoisotopicMass(f) < 2000) {
      sum_formulas.push_back(f);
    }
  }

  utils::IsotopePatternDB db{sum_formulas, target_adducts};
  utils::OrbitrapProfile orbitrap{resolution};
  std::ios_base::sync_with_stdio(true);
  db.useProgressBar(true);
  db.computeIsotopePatterns(orbitrap, max_peaks);
  db.save(output_file);

  std::cout << "Isotope patterns have been saved to " << output_file << std::endl;

  if (options.count("generate-decoy")) {
    auto sign = positive_mode ? "+" : "-";
    for (const auto& item: ms::periodic_table)
      if (adducts_set.find(sign + item.first) == adducts_set.end())
        decoy_adducts.push_back(sign + item.first);
    std::cout << "Generating decoy database..." << std::endl;
    utils::IsotopePatternDB decoy_db{sum_formulas, decoy_adducts};
    utils::OrbitrapProfile orbitrap{resolution};
    std::ios_base::sync_with_stdio(true);
    decoy_db.useProgressBar(true);
    decoy_db.computeIsotopePatterns(orbitrap, max_peaks);
    decoy_db.save(output_file + ".decoy");

    std::cout << "Isotope patterns have been saved to " << output_file << ".decoy" << std::endl;
  }

  return 0;
}
