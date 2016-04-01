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
#include <memory>

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
    if (a.empty())
      continue;
    bool has_sign = a[0] == '+' || a[0] == '-';
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

void saveIsotopeDB(utils::IsotopePatternDB& db, std::string output_fn) {
  db.save(output_fn);
  std::cout << "Isotope patterns have been saved to " << output_fn << std::endl;
}

int isocalc_main(int argc, char** argv) {

  double resolution;
  unsigned max_peaks;

  std::string input_file, output_file;
  std::string adducts_str, decoy_adducts_str;
  std::string instrument_type;

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
    ("instrument",     "Instrument type (orbitrap|fticr)",
     cxxopts::value<std::string>(instrument_type)->default_value("orbitrap"))
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

  using InstrumentProfilePtr = std::unique_ptr<utils::InstrumentProfile>;
  InstrumentProfilePtr instrument;
  if (instrument_type == "orbitrap")
    instrument = InstrumentProfilePtr(new utils::OrbitrapProfile{resolution});
  else if (instrument_type == "fticr")
    instrument = InstrumentProfilePtr(new utils::FTICRProfile(resolution));
  else {
    std::cerr << "unknown instrument type: " << instrument_type << std::endl;
    return -1;
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
  std::ios_base::sync_with_stdio(true);
  db.useProgressBar(true);
  db.computeIsotopePatterns(*instrument, max_peaks);
  saveIsotopeDB(db, output_file);

  if (!decoy_adducts.empty()) {
    std::cout << "Generating decoy database using adducts ";
    printAdducts(decoy_adducts);
    std::cout << "..." << std::endl;
    utils::IsotopePatternDB decoy_db{sum_formulas, decoy_adducts};
    std::ios_base::sync_with_stdio(true);
    decoy_db.useProgressBar(true);
    decoy_db.computeIsotopePatterns(*instrument, max_peaks);
    saveIsotopeDB(decoy_db, output_file + ".decoy");
  }

  return 0;
}
