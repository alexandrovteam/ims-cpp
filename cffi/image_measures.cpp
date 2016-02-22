#include "cffi/common.hpp"
#include "ims/image_measures.hpp"
#include "ms/isotope_pattern.hpp"

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
												 const ms::IsotopePattern&);

template <typename T, MultiF f>
double datacube_score(T** images, int n, int width, int height,
											void* isotope_pattern)
{
	std::vector<ims::ImageF> img_vec;
	for (size_t i = 0; i < n; i++)
		img_vec.push_back(makeImage(images[i], width, height));
	auto p = static_cast<ms::IsotopePattern*>(isotope_pattern);
	return f(img_vec, *p);
}

extern "C" {
	DLL_PUBLIC double measure_of_chaos_f(float* image, int width, int height, int n_levels) {
		return measure_of_chaos<float>(image, width, height, n_levels);
	}

	DLL_PUBLIC double measure_of_chaos_d(double* image, int width, int height, int n_levels) {
		return measure_of_chaos<double>(image, width, height, n_levels);
	}

	DLL_PUBLIC double iso_img_correlation_f(float** images, int n, int width, int height, void* isotope_pattern)
	{
		return datacube_score<float, isotopeImageCorrelation>(images, n, width, height, isotope_pattern);
	}

	DLL_PUBLIC double iso_img_correlation_d(double** images, int n, int width, int height, void* isotope_pattern)
	{
		return datacube_score<double, isotopeImageCorrelation>(images, n, width, height, isotope_pattern);
	}

	DLL_PUBLIC double pattern_match_f(float** images, int n, int width, int height, void* isotope_pattern)
	{
		return datacube_score<float, isotopePatternMatch>(images, n, width, height, isotope_pattern);
	}

	DLL_PUBLIC double pattern_match_d(double** images, int n, int width, int height, void* isotope_pattern)
	{
		return datacube_score<double, isotopePatternMatch>(images, n, width, height, isotope_pattern);
	}
}
