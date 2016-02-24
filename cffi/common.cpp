#include "cffi/common.hpp"

#include <string>
#include <mutex>

namespace cffi {
  static std::string last_error;

  static std::mutex m_;

  void setErrorMessage(const std::string& e) {
    std::lock_guard<std::mutex> lock_guard(m_);
    last_error = e;
  }
}

extern "C" {

  IMS_EXTERN const char* ims_strerror() {
    return cffi::last_error.c_str();
  }
}
