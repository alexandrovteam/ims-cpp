#pragma once

#include <cmath>

namespace ms {

class InstrumentProfile {
 public:
  virtual double resolvingPowerAt(double mz) const = 0;
  virtual ~InstrumentProfile() {}
};

class OrbitrapProfile : public InstrumentProfile {
  double res_power_;
  double at_;

 public:
  OrbitrapProfile(double resolving_power, double at = 200)
      : res_power_(resolving_power), at_(at) {}

  double resolvingPowerAt(double mz) const { return res_power_ * sqrt(at_ / mz); }
};

class FTICRProfile : public InstrumentProfile {
  double res_power_;
  double at_;

 public:
  FTICRProfile(double resolving_power, double at = 200)
      : res_power_(resolving_power), at_(at) {}

  double resolvingPowerAt(double mz) const { return res_power_ * at_ / mz; }
};

class TOFProfile : public InstrumentProfile {
  double res_power_;

 public:
  TOFProfile(double resolving_power) : res_power_(resolving_power) {}

  double resolvingPowerAt(double mz) const { return res_power_; }
};
}
