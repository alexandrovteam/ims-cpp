#include "ms/isocalc.hpp"
#include "imzb/reader.hpp"

#include "utils/isotope_pattern_db.hpp"

#include "cxxopts.hpp"
extern "C" {
#include "progressbar.h"
}

#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <string>
#include <mutex>

struct Metrics {
  double img_corr;
  double iso_corr;
  double moc;
};

int detect_main(int argc, char** argv) {
  std::string db_filename, imzb_filename;
  std::string output_filename;
  double ppm;
  bool remove_hotspots;

  cxxopts::Options options("ims detect", " <isotope_patterns.db> <input.imzb>");
  options.add_options()
    ("ppm", "m/z-window half-width in ppm",
     cxxopts::value<double>(ppm)->default_value("3.0"))
    ("out", "",
     cxxopts::value<std::string>(output_filename)->default_value("/dev/stdout"))
    ("remove-hotspots", "apply hotspot removal prior to computing metrics",
     cxxopts::value<bool>(remove_hotspots)->default_value("true"))
    ("help", "Print help");

  options.add_options("hidden")
    ("db_filename", "", cxxopts::value<std::string>(db_filename))
    ("input_file", "", cxxopts::value<std::string>(imzb_filename));

  options.parse_positional(std::vector<std::string>{"db_filename", "input_file"});

  options.parse(argc, argv);

  if (options.count("help") || imzb_filename.empty() || db_filename.empty()) {
    std::cout << options.help({""}) << std::endl;
    return 0;
  }

  utils::IsotopePatternDB db(db_filename);
  imzb::ImzbReader reader(imzb_filename);
  std::ofstream out(output_filename);

  auto& keys = db.keys();
  std::vector<Metrics> metrics(keys.size());

  std::vector<ims::ImageF> images;
  std::vector<float> hotspot_removal_buf;

  const int BAR_STEP = 50;
  auto bar = progressbar_new("Computing scores", keys.size() / BAR_STEP);
  std::mutex bar_mutex;

#pragma omp parallel for private(images, hotspot_removal_buf) shared(bar, bar_mutex)
  for (size_t i = 0; i < keys.size(); i++) {
    std::string f = keys[i].first;
    std::string adduct = keys[i].second;

    auto p = db(f, adduct);
    if (p.size() == 0)
      continue;

    if (images.size() < p.size()) {
      for (size_t i = images.size(); i < p.size(); i++)
        images.emplace_back(reader.height(), reader.width());
    }

    if (remove_hotspots && hotspot_removal_buf.empty())
      hotspot_removal_buf.resize(reader.width() * reader.height());

    for (size_t j = 0; j < p.size(); ++j) {
      reader.readImage(p.masses[j], ppm, images[j].rawPtr());
      if (remove_hotspots)
        images[j].removeHotspots(99.0, &hotspot_removal_buf[0]);
    }

    auto img_corr = ims::isotopeImageCorrelation(&images[0], p.size(), p);
    auto iso_corr = ims::isotopePatternMatch(&images[0], p.size(), p);
    auto moc = ims::measureOfChaos(images[0], 30);
    metrics[i] = Metrics{img_corr, iso_corr, moc};

    std::lock_guard<std::mutex> lock(bar_mutex);
    if ((i + 1) % BAR_STEP == 0)
      progressbar_inc(bar);
  }

  progressbar_finish(bar);

  const auto sep = ",";

  out << "formula,adduct,img,iso,moc\n";

  for (size_t i = 0; i < keys.size(); ++i) {
    std::string f = keys[i].first;
    std::string adduct = keys[i].second;
    out << f << sep << adduct << sep <<
      metrics[i].img_corr << sep << metrics[i].iso_corr << sep << metrics[i].moc << "\n";
  }

  return 0;
}
