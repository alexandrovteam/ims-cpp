#pragma once

#include <sstream>
#include <vector>
#include <string>

namespace utils {

bool readDouble(std::istringstream& in, double& num);

struct Metrics {
  std::string sf;
  std::string adduct;
  double img_corr;
  double iso_corr;
  double chaos;
  double median_mz_shift;
  double num_pixels;

  double msm() const { return img_corr * iso_corr * chaos; }

  static std::string header() { return "formula,adduct,img,iso,moc,mz_shift,n_pixels"; }

  static bool read(std::istream& in, Metrics& m) {
    static std::string line;
    if (!std::getline(in, line)) return false;

    std::istringstream ss{line};
    if (!std::getline(ss, m.sf, ',')) return false;
    if (!std::getline(ss, m.adduct, ',')) return false;
    if (!readDouble(ss, m.img_corr)) return false;
    if (!readDouble(ss, m.iso_corr)) return false;
    if (!readDouble(ss, m.chaos)) return false;
    if (!readDouble(ss, m.median_mz_shift)) return false;
    if (!readDouble(ss, m.num_pixels)) return false;
    return true;
  }
};

std::ostream& operator<<(std::ostream& os, const Metrics& m);
}
