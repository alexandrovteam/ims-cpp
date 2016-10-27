#include "imzb/reader.hpp"
#include "ims/common.hpp"

#include <vector>
#include <algorithm>
#include <memory>

namespace imzb {

// [left, right] interval
struct MzBin {
  double left, right;
  uint64_t count;
  uint64_t core_count;
  double sum, sumsq;
  double intensity;
};

// DBScan for sorted one-dimensional data
class DBScan {
  int min_pts_;
  double eps_;
  int core_pos_;
  std::vector<MzBin> bins_;

  std::vector<double> buf_;
  std::vector<ims::Peak> peak_buf_;
  size_t curr_right_pos_;
  MzBin curr_bin_;

 public:
  size_t n_points_total_ = 0;
  size_t n_points_in_clusters_ = 0;
  DBScan(int minPts, double eps);

  void put(const ims::Peak& peak);

  const std::vector<MzBin>& bins() const;
};

DBScan dbscan(ImzbReader* reader, uint32_t minPts, double eps);
}
