char* ims_strerror();

typedef void* ImzbReader;
ImzbReader imzb_reader_new(char*);
void imzb_reader_free(ImzbReader);
int imzb_reader_height(ImzbReader);
int imzb_reader_width(ImzbReader);
double imzb_reader_min_mz(ImzbReader);
double imzb_reader_max_mz(ImzbReader);
int imzb_reader_image(ImzbReader, double mz, double ppm, float* out);
int imzb_reader_centroided_image(ImzbReader, double mz, double ppm, float* out);

typedef void* IsotopePattern;
IsotopePattern isotope_pattern_new(int n, double* masses, double* abundances);
IsotopePattern isotope_pattern_new_from_sf(char* formula, double thr, double fft_thr);
IsotopePattern isotope_pattern_new_from_raw(int n, double* masses, double* intensities, int window_size);
IsotopePattern isotope_pattern_copy(IsotopePattern);
float isotope_pattern_envelope(IsotopePattern, double resolution, double mz);
int isotope_pattern_envelope_plot(IsotopePattern, double resolution, double* mzs, int n, float* out);
IsotopePattern isotope_pattern_centroids(IsotopePattern, double resolution, double min_abundance,
                                         int points_per_fwhm);
int isotope_pattern_size(IsotopePattern);
void isotope_pattern_masses(IsotopePattern, double*);
void isotope_pattern_abundances(IsotopePattern, double*);
void isotope_pattern_add_charge(IsotopePattern, int);
void isotope_pattern_trim(IsotopePattern, unsigned);
void isotope_pattern_free(IsotopePattern);

double measure_of_chaos_f(float* image, int width, int height, int n_levels);
double measure_of_chaos_d(double* image, int width, int height, int n_levels);

double iso_img_correlation_f(float** images, int n, int width, int height, double* isotope_abundances);
double iso_img_correlation_d(double** images, int n, int width, int height, double* isotope_abundances);

double pattern_match_f(float** images, int n, int width, int height, double* isotope_abundances);
double pattern_match_d(double** images, int n, int width, int height, double* isotope_abundances);
