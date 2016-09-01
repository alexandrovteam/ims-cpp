#pragma once

#ifdef IMS_PROFILE_ISOCALC
#include <fstream>
#include <chrono>
#include <string>
#include <sstream>
#include <iostream>

namespace utils {

// not thread-safe
class ProfilingTimer {
  std::ofstream out_;

  typedef std::chrono::high_resolution_clock Clock;
  Clock clock_;
  std::stringstream message_;
  std::chrono::time_point<Clock> start_;
  bool started_;

  ProfilingTimer(const std::string& output_file) : started_(false) {
    out_.open(output_file, std::ofstream::out | std::ofstream::app);
  }

  ~ProfilingTimer() { out_.close(); }

 public:
  inline static ProfilingTimer& instance() {
    static ProfilingTimer timer("ims_isocalc_profile.log");
    return timer;
  }

  void start() {
    start_ = clock_.now();
    started_ = true;
  }

  void reset() {
    std::stringstream().swap(message_);
    started_ = false;
  }

  template <typename T>
  ProfilingTimer& operator<<(const T& obj) {
    message_ << obj;
    return *this;
  }

  void stop() {
    if (!started_) return;

    started_ = false;
    auto end_ = clock_.now();
    auto diff = end_ - start_;
    auto milli = std::chrono::duration<double, std::milli>(diff).count();

    out_ << "[ " << message_.str() << " ] ~ " << milli << std::endl;
  }
};
}
#else  // no-op logger
namespace utils {
struct ProfilingTimer {
  inline static ProfilingTimer& instance() {
    static ProfilingTimer timer;
    return timer;
  }

  void start() {}
  void stop() {}
  void reset() {}

  template <typename T>
  ProfilingTimer& operator<<(const T&) {
    return *this;
  }
};
}
#endif
