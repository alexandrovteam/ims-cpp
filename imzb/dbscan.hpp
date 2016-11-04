#include "imzb/reader.hpp"
#include "ims/common.hpp"

#include <vector>
#include <algorithm>
#include <memory>
#include <functional>

namespace imzb {

// [left, right] interval
struct MzBin {
  double left, right;
  uint64_t count;
  uint64_t core_count;
  double intensity;
  double median_mz;
  double mean, sd;
};

// DBScan for sorted one-dimensional data
class DBScan {
  int min_pts_;
  std::function<double(double)> eps_;
  int core_pos_;
  double curr_eps_ = -1;
  std::vector<MzBin> bins_;

  std::vector<double> buf_;
  std::vector<ims::Peak> peak_buf_;
  size_t curr_right_pos_;
  MzBin curr_bin_;
  double curr_sum_ = 0.0, curr_sumsq_ = 0.0;

 public:
  size_t n_points_total_ = 0;
  size_t n_points_in_clusters_ = 0;
  size_t eps_update_interval_ = 50;

  DBScan(int minPts, double eps);
  DBScan(int minPts, std::function<double(double)> eps);

  void put(const ims::Peak& peak);

  const std::vector<MzBin>& bins() const;
};

DBScan dbscan(ImzbReader* reader, uint32_t minPts, double eps);
DBScan dbscan(ImzbReader* reader, uint32_t minPts, double eps,
              double min_mz, double max_mz);

DBScan dbscan(ImzbReader* reader, uint32_t minPts, std::function<double(double)> eps);
DBScan dbscan(ImzbReader* reader, uint32_t minPts, std::function<double(double)> eps,
              double min_mz, double max_mz);

}
