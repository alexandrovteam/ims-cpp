#pragma once

#if defined _WIN32
  #define DLL_PUBLIC __declspec(dllexport)
#else
  #define DLL_PUBLIC __attribute__ ((visibility ("default")))
#endif

extern "C" {
  DLL_PUBLIC const char* ims_strerror();
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
