#include <string>

namespace utils {

  bool endsWith(const std::string& s, const std::string& tail) {
    if (s.length() < tail.length())
      return false;
    return s.substr(s.length() - tail.length()) == tail;
  }

}
