#include "cffi/common.hpp"
#include "imzb/reader.hpp"
#include "imzb/dbscan.hpp"

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

IMS_EXTERN int imzb_reader_slice(imzb::ImzbReader* reader,
                                 double min_mz, double max_mz, ims::Peak** out) {
  return wrap_catch<int>(-1, [&]() -> int {
    auto slice = reader->slice(min_mz, max_mz);
    *out = new ims::Peak[slice.size()];
    std::copy(slice.begin(), slice.end(), *out);
    return slice.size();
  });
}

static int copyBins(const imzb::DBScan& dbscan, imzb::MzBin** dst) {
  const auto& src = dbscan.bins();
  *dst = new imzb::MzBin[src.size()];
  std::copy(src.begin(), src.end(), *dst);
  return src.size();
}

IMS_EXTERN int imzb_reader_dbscan(imzb::ImzbReader* reader, int minPts, double eps,
                                  imzb::MzBin** out) {
  return wrap_catch<int>(-1, [&]() -> int {
    return copyBins(imzb::dbscan(reader, uint32_t(minPts), eps), out);
  });
}

IMS_EXTERN int imzb_reader_dbscan2(imzb::ImzbReader* reader, int minPts, double eps,
                                   double min_mz, double max_mz, imzb::MzBin** out) {
  return wrap_catch<int>(-1, [&]() -> int {
    return copyBins(imzb::dbscan(reader, uint32_t(minPts), eps, min_mz, max_mz), out);
  });
}
IMS_EXTERN int imzb_reader_dbscan3(imzb::ImzbReader* reader, int minPts, double (*eps)(double),
                                   imzb::MzBin** out) {
  return wrap_catch<int>(-1, [&]() -> int {
    return copyBins(imzb::dbscan(reader, uint32_t(minPts), eps), out);
  });
}

IMS_EXTERN int imzb_reader_dbscan4(imzb::ImzbReader* reader, int minPts, double (*eps)(double),
                                   double min_mz, double max_mz, imzb::MzBin** out) {
  return wrap_catch<int>(-1, [&]() -> int {
    return copyBins(imzb::dbscan(reader, uint32_t(minPts), eps, min_mz, max_mz), out);
  });
}

}
