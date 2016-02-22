#pragma once

#include <fstream>

template <typename T>
bool binary_read(std::ifstream& stream, T& value) {
  return bool(stream.read(reinterpret_cast<char*>(&value), sizeof(value)));
}

template <typename T>
void binary_write(std::ofstream& stream, T& value) {
  stream.write(reinterpret_cast<char*>(&value), sizeof(value));
}
