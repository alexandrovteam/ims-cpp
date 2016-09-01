#pragma once

#include "ims/image.hpp"
#include "ims/common.hpp"

#include <vector>
#include <valarray>
#include <cassert>

namespace ims {

double isotopeImageCorrelation(const ims::ImageF* images, size_t n,
    const std::vector<double>& abundances, bool zero_mean_normalization = false);

double isotopeImageCorrelation(const std::vector<ims::ImageF>& images,
    const std::vector<double>& abundances, bool zero_mean_normalization = false);

double isotopePatternMatch(
    const ims::ImageF* images, size_t n, const std::vector<double>& pattern);

double isotopePatternMatch(
    const std::vector<ims::ImageF>& images, const std::vector<double>& pattern);

double measureOfChaos(const ims::ImageF& image, size_t n_levels = 32);

double medianMzShift(const std::vector<ims::Peak>& peaks, double mz);
}
