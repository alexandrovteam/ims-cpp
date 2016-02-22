#pragma once

#include <vector>
#include <cstdlib>
#include <initializer_list>
#include <cassert>

namespace ms {
  struct IsotopePattern {
    std::vector<double> masses;
    std::vector<double> abundances;

    IsotopePattern() {}
    IsotopePattern(double mass) : masses{mass}, abundances{1.0} {}
    IsotopePattern(std::initializer_list<double> masses,
                   std::initializer_list<double> abundances) : masses{masses}, abundances{abundances}
    {}

    IsotopePattern(std::vector<double> masses,
                   std::vector<double> abundances) : masses{masses}, abundances{abundances}
    {}

    // shifts masses accordingly to the number of added/subtracted electrons
    void addCharge(int charge);

    IsotopePattern& charged(int charge) {
      addCharge(charge);
      return *this;
    }

    IsotopePattern& trimmed(size_t new_size) {
      if (new_size < size()) {
        masses.resize(new_size);
        abundances.resize(new_size);
      }
      return *this;
    }

    IsotopePattern copy() const {
      IsotopePattern p{masses, abundances};
      return p;
    }

    bool isUnit() const { return masses.size() == 1 && masses[0] == 0.0; }

    ms::IsotopePattern multiply(const IsotopePattern& other, double threshold = 0.0) const;

    double envelope(double resolution, double mz) const;

    ms::IsotopePattern centroids(double resolution,
                                 double min_abundance = 1e-4,
                                 size_t points_per_fwhm = 25) const;

    // sorts by decreasing intensity and scales max intensity to 1
    void normalize();

    size_t size() const { return masses.size(); }
  };
}
