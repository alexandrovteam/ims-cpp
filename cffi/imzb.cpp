#include "cffi/common.hpp"
#include "imzb/reader.hpp"

using namespace cffi;

extern "C" {
IMS_EXTERN imzb::ImzbReader* imzb_reader_new(const char* filename) {
  return wrap_catch<imzb::ImzbReader*>(
      nullptr, [&]() { return new imzb::ImzbReader(filename); });
}

IMS_EXTERN void imzb_reader_free(imzb::ImzbReader* reader) {
  delete reader;
}

IMS_EXTERN int imzb_reader_height(imzb::ImzbReader* reader) {
  return reader->height();
}

IMS_EXTERN int imzb_reader_width(imzb::ImzbReader* reader) {
  return reader->width();
}

IMS_EXTERN int imzb_reader_image(
    imzb::ImzbReader* reader, double mz, double ppm, float* outbuf) {
  return wrap_catch<int>(-1, [&]() -> int {
    reader->readImage(mz, ppm, outbuf);
    return 0;
  });
}

IMS_EXTERN int imzb_reader_centroided_image(
    imzb::ImzbReader* reader, double mz, double ppm, float* outbuf) {
  return wrap_catch<int>(-1, [&]() -> int {
    reader->readCentroidedImage(mz, ppm, outbuf);
    return 0;
  });
}

IMS_EXTERN double imzb_reader_min_mz(imzb::ImzbReader* reader) {
  return reader->minMz();
}

IMS_EXTERN double imzb_reader_max_mz(imzb::ImzbReader* reader) {
  return reader->maxMz();
}
}
