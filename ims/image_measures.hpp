#pragma once

#include "ims/image.hpp"
#include "ms/isotope_pattern.hpp"

#include <vector>
#include <valarray>
#include <cassert>
#include <iostream>

namespace ims {

  double isotopeImageCorrelation(const ims::ImageF* images, size_t n,
                                 const ms::IsotopePattern& pattern);

  double isotopeImageCorrelation(const std::vector<ims::ImageF>& images,
                                 const ms::IsotopePattern& pattern);

  double isotopePatternMatch(const ims::ImageF* images, size_t n,
                             const ms::IsotopePattern& pattern);

  double isotopePatternMatch(const std::vector<ims::ImageF>& images,
      const ms::IsotopePattern& pattern);

  double measureOfChaos(const ims::ImageF& image, size_t n_levels=32);
}
