#include "imzb/dbscan.hpp"

namespace imzb {

DBScan::DBScan(int minPts, double eps) : min_pts_(minPts), eps_(eps) {}

const std::vector<MzBin>& DBScan::bins() const { return this->bins_; }

void DBScan::put(const ims::Peak& peak) {
  ++n_points_total_;
  buf_.push_back(peak.mz);
  peak_buf_.push_back(peak);

  if (buf_.size() < min_pts_) {
    core_pos_ = -1;
    curr_bin_.left = -1;
    curr_bin_.right = -1;
    curr_bin_.count = -1;
    return;
  }

  size_t new_core_pos = buf_.size() - min_pts_;
  double new_core = buf_[new_core_pos];
  auto left_it = std::lower_bound(buf_.begin(), buf_.begin() + new_core_pos,
                                  new_core - eps_);
  auto right_it = std::upper_bound(buf_.begin() + new_core_pos, buf_.end(),
                                   new_core + eps_) - 1;
  int count = right_it - left_it + 1;
  bool is_core = count >= min_pts_;
  int right_pos = right_it - buf_.begin();

  if (is_core && core_pos_ >= 0 &&
      (new_core - buf_[core_pos_] <= eps_ || *left_it < curr_bin_.right)) {
    core_pos_ = new_core_pos;
    if (*right_it > curr_bin_.right) {
      curr_bin_.right = *right_it;
      curr_bin_.count += right_pos - curr_right_pos_;
      for (int i = curr_right_pos_ + 1; i <= right_pos; i++) {
        curr_bin_.sum += buf_[i];
        curr_bin_.sumsq += buf_[i] * buf_[i];
        curr_bin_.intensity += peak_buf_[i].intensity;
      }
      curr_right_pos_ = right_pos;
      ++curr_bin_.core_count;
    }
  } else if (is_core) {
    MzBin new_bin;
    new_bin.count = count;
    new_bin.left = *left_it;
    new_bin.right = *right_it;
    new_bin.sum = new_bin.sumsq = new_bin.intensity = 0.0;
    new_bin.core_count = 1;

    for (auto it = left_it; it <= right_it; ++it) {
      new_bin.sum += *it;
      new_bin.sumsq += (*it) * (*it);
      new_bin.intensity += peak_buf_[it - buf_.begin()].intensity;
    }

    if (core_pos_ < 0) {
      core_pos_ = new_core_pos;
      curr_right_pos_ = right_pos;
    } else {
      bins_.push_back(curr_bin_);
      n_points_in_clusters_ += curr_bin_.count;
      std::vector<double> trimmed(left_it, buf_.end());
      std::vector<ims::Peak> peak_trimmed(
          peak_buf_.begin() + (left_it - buf_.begin()), peak_buf_.end());
      int shift = left_it - buf_.begin();
      buf_.swap(trimmed);
      peak_buf_.swap(peak_trimmed);
      core_pos_ = new_core_pos - shift;
      curr_right_pos_ = right_pos - shift;
    }

    curr_bin_ = new_bin;
  }
}

DBScan dbscan(ImzbReader* reader, uint32_t minPts, double eps) {
  DBScan dbscan(minPts, eps);

  ims::Peak peak;
  while (reader->readNext(peak))
    dbscan.put(peak);

  return dbscan;
}
}
