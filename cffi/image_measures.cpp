#include "cffi/common.hpp"
#include "ims/image_measures.hpp"

#include <cmath>

using namespace ims;

template <typename T>
ims::ImageF makeImage(T* image, int width, int height) {
  ims::ImageF img(height, width);
  auto p = img.rawPtr();
  for (size_t i = 0, n = height * width; i < n; i++)
    p[i] = std::isnan(image[i]) ? -1.0 : image[i];
  return img;
}

template <typename T>
double measure_of_chaos(T* image, int width, int height, int n_levels) {
  return measureOfChaos(makeImage(image, width, height), n_levels);
}

typedef double (*MultiF)(const std::vector<ims::ImageF>&,
                         const std::vector<double>&);

template <typename T, MultiF f>
double datacube_score(T** images, int n, int width, int height,
                      double* isotope_abundances)
{
  std::vector<ims::ImageF> img_vec;
  for (size_t i = 0; i < n; i++)
    img_vec.push_back(makeImage(images[i], width, height));
  std::vector<double> abundances(isotope_abundances, isotope_abundances + n);
  return f(img_vec, abundances);
}

extern "C" {
  IMS_EXTERN double measure_of_chaos_f(float* image, int width, int height, int n_levels) {
    return measure_of_chaos<float>(image, width, height, n_levels);
  }

  IMS_EXTERN double measure_of_chaos_d(double* image, int width, int height, int n_levels) {
    return measure_of_chaos<double>(image, width, height, n_levels);
  }

  IMS_EXTERN double iso_img_correlation_f(float** images, int n, int width, int height, double* abundances)
  {
    return datacube_score<float, isotopeImageCorrelation>(images, n, width, height, abundances);
  }

  IMS_EXTERN double iso_img_correlation_d(double** images, int n, int width, int height, double* abundances)
  {
    return datacube_score<double, isotopeImageCorrelation>(images, n, width, height, abundances);
  }

  IMS_EXTERN double pattern_match_f(float** images, int n, int width, int height, double* abundances)
  {
    return datacube_score<float, isotopePatternMatch>(images, n, width, height, abundances);
  }

  IMS_EXTERN double pattern_match_d(double** images, int n, int width, int height, double* abundances)
  {
    return datacube_score<double, isotopePatternMatch>(images, n, width, height, abundances);
  }
}
