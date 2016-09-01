#include "utils/metrics.hpp"

#include <cmath>
#include <iomanip>

namespace utils {

bool readDouble(std::istringstream& in, double& num) {
  static std::string tmp;
  if (!std::getline(in, tmp, ',')) return false;
  num = std::stod(tmp);
  if (std::isnan(num)) num = 0.0;
  return true;
}

std::ostream& operator<<(std::ostream& os, const Metrics& m) {
  os << m.sf << "," << m.adduct << "," << std::setprecision(3) << m.img_corr << ","
     << m.iso_corr << "," << m.chaos << "," << m.median_mz_shift << "," << m.num_pixels;
  return os;
}
}
