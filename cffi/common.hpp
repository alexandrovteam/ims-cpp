#pragma once

#if defined _WIN32
  #define IMS_EXTERN __declspec(dllexport)
#else
  #define IMS_EXTERN __attribute__ ((visibility ("default")))
#endif

extern "C" {
  IMS_EXTERN const char* ims_strerror();
}

#include <string>
#include <exception>
#include <functional>

namespace cffi {
  void setErrorMessage(const std::string& e);

  template <typename R>
  R wrap_catch(R on_error, std::function<R()>&& setter) {
    try {
      return setter();
    } catch (std::exception& e) {
      setErrorMessage(e.what());
      return on_error;
    }
  }

  template <typename T>
  T* heapify(const T& value) {
    return new T(value);
  }
}
