char* ims_strerror();

typedef struct {
  uint32_t x;
  uint32_t y;
  uint32_t z;
} Position;

typedef struct {
  Position coords;
  double mz;
  float intensity;
} Peak;

typedef struct {
  double left, right;
  uint64_t count;
  uint64_t core_count;
  double sum, sumsq;
  double intensity;
} MzBin;

typedef void* ImzbReader;
ImzbReader imzb_reader_new(char*);
void imzb_reader_free(ImzbReader);
int imzb_reader_height(ImzbReader);
int imzb_reader_width(ImzbReader);
double imzb_reader_min_mz(ImzbReader);
double imzb_reader_max_mz(ImzbReader);
int imzb_reader_image(ImzbReader, double mz, double ppm, float* out);
int imzb_reader_centroided_image(ImzbReader, double mz, double ppm, float* out);
int imzb_reader_slice(ImzbReader, double min_mz, double max_mz, Peak** out);
int imzb_reader_dbscan(ImzbReader, int minPts, double eps, MzBin** out);
int imzb_reader_dbscan2(ImzbReader, int minPts, double eps, double min_mz, double max_mz, MzBin** out);

typedef void* InstrumentProfile;
InstrumentProfile instrument_profile_new(const char* type, double resolving_power, double at_mz);
void instrument_profile_free(InstrumentProfile);
double instrument_resolving_power_at(InstrumentProfile, double mz);

typedef void* Spectrum;
Spectrum spectrum_new(int n, double* masses, double* intensities);
Spectrum spectrum_new_from_sf(char* formula, double thr, double fft_thr);
Spectrum spectrum_new_from_raw(int n, double* masses, float* intensities, int window_size);
Spectrum spectrum_copy(Spectrum);
float spectrum_envelope(Spectrum, InstrumentProfile, double mz);
int spectrum_envelope_plot(Spectrum, InstrumentProfile, double* mzs, int n, float* out);
Spectrum spectrum_envelope_centroids(Spectrum, InstrumentProfile, double min_abundance,
                                     int points_per_fwhm);
int spectrum_size(Spectrum);
void spectrum_masses(Spectrum, double*);
void spectrum_intensities(Spectrum, double*);
void spectrum_multiply_inplace(Spectrum, double);
void spectrum_add_inplace(Spectrum, Spectrum);
void spectrum_add_charge(Spectrum, int);
void spectrum_normalize(Spectrum);
void spectrum_sort_by_mass(Spectrum);
void spectrum_sort_by_intensity(Spectrum);
Spectrum spectrum_convolve(Spectrum, Spectrum);
void spectrum_trim(Spectrum, unsigned);
void spectrum_trim_intensity(Spectrum, double);
void spectrum_free(Spectrum);

double measure_of_chaos_f(float* image, int width, int height, int n_levels);
double measure_of_chaos_d(double* image, int width, int height, int n_levels);

double iso_img_correlation_f(float** images, int n, int width, int height, double* isotope_abundances);
double iso_img_correlation_d(double** images, int n, int width, int height, double* isotope_abundances);

double pattern_match_f(float** images, int n, int width, int height, double* isotope_abundances);
double pattern_match_d(double** images, int n, int width, int height, double* isotope_abundances);

void free(void*);
