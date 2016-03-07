#pragma once

#include "ims/ims.hpp"

#include "hdf5.h"
#include "hdf5_hl.h"

#include <vector>
#include <set>
#include <string>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <stdexcept>

namespace scils {

class H5Reader : public ims::AbstractReader {
  hid_t file_;

  std::vector<double> x_coords_, y_coords_, z_coords_;
  size_t nrow_, ncol_;
  size_t curr_idx_;
  size_t grid_size_;

  std::vector<float> intensities_buf_;
  std::vector<double> mzs_buf_;

  void readCoordinates() {
    hsize_t ndims[2] = {0};
    std::string dataset = version() >= 5 ? "/Registrations/0/Coordinates" : "Coordinates";
    H5LTget_dataset_info(file_, dataset.c_str(), ndims, NULL, NULL);
    assert(ndims[0] == 3);
    auto n = ndims[1];
    std::vector<double> coords(3 * n);

    H5LTread_dataset_double(file_, dataset.c_str(), &coords[0]);

    x_coords_.assign(coords.begin(), coords.begin() + n);
    y_coords_.assign(coords.begin() + n, coords.begin() + 2 * n);
    z_coords_.assign(coords.begin() + 2 * n, coords.begin() + 3 * n);
  }

  size_t normalizeGrid(std::vector<double>& values) {
    double sum_diff = 0.0;
    size_t n = 0;
    std::set<double> v(values.begin(), values.end());
    auto it = v.begin();
    for (;;) {
      auto curr_value = *it;
      if (++it == v.end())
        break;
      auto diff = *it - curr_value;
      if (diff > 1e-4)
        sum_diff += diff, ++n;
    }
    auto step = n > 0 ? sum_diff / n : 1.0;
    for (auto& val: values)
      val = (val - *v.begin()) / step;

    return 1 + *std::max_element(values.begin(), values.end());
  }

public:
  H5Reader(const std::string& filename) :
    file_(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT))
  {
    readCoordinates();
    ncol_ = normalizeGrid(x_coords_);
    nrow_ = normalizeGrid(y_coords_);
    normalizeGrid(z_coords_);
    curr_idx_ = 0;
  }

  ~H5Reader() {
    H5Fclose(file_);
  }

  int version() {
    int version;
    if (H5LTread_dataset_int(file_, "/Version", &version) < 0)
      throw std::runtime_error("couldn't read the format version");
    return version;
  }

  uint32_t width() const {
    return ncol_;
  }

  uint32_t height() const {
    return nrow_;
  }

  bool readNextSpectrum(ims::Spectrum& spectrum) {
    if (curr_idx_ >= x_coords_.size())
      return false;

    auto dataset = "Spots/" + std::to_string(curr_idx_) + "/InitialMeasurement/Intensities";

    if (intensities_buf_.empty()) {
      hsize_t ndims[1] = {0};
      H5LTget_dataset_info(file_, dataset.c_str(), ndims, NULL, NULL);
      grid_size_ = ndims[0];
    }

    spectrum.intensities.resize(grid_size_);
    H5LTread_dataset_float(file_, dataset.c_str(), &spectrum.intensities[0]);

    if (mzs_buf_.empty()) {
      mzs_buf_.resize(grid_size_);
      dataset = "SamplePositions/GlobalMassAxis/SamplePositions";
      H5LTread_dataset_double(file_, dataset.c_str(), &mzs_buf_[0]);
    }

    spectrum.mzs = mzs_buf_;
    spectrum.coords = ims::Position{uint32_t(nrow_ - 1 - y_coords_[curr_idx_]),
                                    uint32_t(x_coords_[curr_idx_]),
                                    uint32_t(z_coords_[curr_idx_])};

    ++curr_idx_;
    return true;
  }
};

}
