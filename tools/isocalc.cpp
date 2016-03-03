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

void printAdducts(const std::vector<std::string>& adducts) {
  for (auto& a: adducts) {
    std::cout << a;
    if (a != adducts.back())
      std::cout << ", ";
  }
}

std::vector<std::string> parseAdducts(const std::string& adducts_str) {
  std::regex regex{","};
  std::set<std::string> tmp{
    std::sregex_token_iterator(adducts_str.begin(), adducts_str.end(), regex, -1),
      std::sregex_token_iterator()};
  std::set<std::string> adducts;
  for (auto& a: tmp) {
    bool has_sign = !a.empty() && (a[0] == '+' || a[0] == '-');
    try {
      sf_parser::parseSumFormula(a);
    } catch (sf_parser::NegativeTotalError&) {
      // ignore this error type for adducts
    } catch (sf_parser::ParseError& e) {
      // die if it raises a parsing error
      std::cerr << "'" << a << "' is not a valid adduct!" << std::endl;
      throw e;
    }
    adducts.insert(has_sign ? a : "+" + a);
  }
  return std::vector<std::string>{adducts.begin(), adducts.end()};
}

int isocalc_main(int argc, char** argv) {

  double resolution;
  unsigned max_peaks;

  std::string input_file, output_file;

  std::string adducts_str, decoy_adducts_str;

  cxxopts::Options options("ims isocalc",
  " <input.txt> <output.db>\n\t\t\twhere input contains one sum formula per line.");
  options.add_options()
    ("resolution",      "Resolving power at m/z=200",
     cxxopts::value<double>(resolution)->default_value("140000.0"))
    ("adducts",         "Comma-separated list of adducts",
     cxxopts::value<std::string>(adducts_str)->default_value("+H,+K,+Na"))
    ("max-peaks",       "Maximum number of peaks to store",
     cxxopts::value<unsigned>(max_peaks)->default_value("5"))
    ("decoy-adducts",  "Adducts for generating a decoy database (stored in a separate file named <output.db.decoy>)",
     cxxopts::value<std::string>(decoy_adducts_str)->default_value(""))
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

  auto target_adducts = parseAdducts(adducts_str);
  auto decoy_adducts = parseAdducts(decoy_adducts_str);

  std::cout << "Resolution @ m/z=200: " << resolution << std::endl;
  std::cout << "Generating isotope pattern database for adducts ";
  printAdducts(target_adducts);
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

  if (!decoy_adducts.empty()) {
    std::cout << "Generating decoy database using adducts ";
    printAdducts(decoy_adducts);
    std::cout << "..." << std::endl;
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
