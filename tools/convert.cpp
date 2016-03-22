#include "imzml/reader.hpp"
#include "utils/string.hpp"

#include "cxxopts.hpp"

#ifdef SCILS_H5
#include "scils/h5reader.hpp"
#endif

#include "imzb/reader.hpp"
#include "imzb/writer.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <memory>
#include <queue>
#include <cstdio>

class Sorter {
  std::string fn_;
  std::vector<ims::Peak> buffer_;
  std::vector<std::string> tmp_filenames_;
  size_t filled_;
  bool closed_;
  uint32_t block_size_;
  std::string compressor_;
  uint8_t comp_level_;

  void dump() {
    std::sort(buffer_.begin(), buffer_.begin() + filled_,
        [](const ims::Peak& a, const ims::Peak& b) { return a.mz < b.mz; }
    );
    std::stringstream tmp_fn;
    tmp_fn << fn_ << "." << tmp_filenames_.size();
    tmp_filenames_.push_back(tmp_fn.str());
    imzb::ImzbWriter writer(tmp_filenames_.back(), block_size_, "blosclz");
    writer.setMask(mask_);
    for (size_t i = 0; i < filled_; ++i)
      writer.writePeak(buffer_[i]);
    writer.close();

    filled_ = 0;
  }

  struct PeakAndFile {
    ims::Peak peak;
    size_t file_index;

    bool operator<(const PeakAndFile& other) const {
      return peak.mz > other.peak.mz;
    }
  };

  void merge() {
    std::cout << "merging files..." << std::endl;
    std::priority_queue<PeakAndFile> queue;

    std::vector<std::shared_ptr<imzb::ImzbReader>> readers;
    imzb::ImzbWriter writer(fn_, block_size_, compressor_, comp_level_);
    writer.setMask(mask_);

    ims::Peak peak;
    for (const auto& tmp_fn: tmp_filenames_) {
      readers.push_back(std::make_shared<imzb::ImzbReader>(tmp_fn));
      if (readers.back()->readNext(peak))
        queue.push(PeakAndFile{peak, readers.size() - 1});
    }

    size_t n = 0;
    while (!queue.empty()) {
      auto item = queue.top();
      writer.writePeak(item.peak);
      ++n;
      queue.pop();

      if (readers[item.file_index]->readNext(peak))
        queue.push(PeakAndFile{peak, item.file_index});
    }

    writer.close();

    std::cout << "removing temporary files" << std::endl;
    for (const auto& r: readers) {
      auto fn = r->filename();
      auto idx_fn = fn + ".idx";
      r->close();
      std::remove(fn.c_str());
      std::remove(idx_fn.c_str());
    }
    std::cout << "done!" << std::endl;
  }

  const imzb::Mask& mask_;
public:
  Sorter(const std::string& filename, const imzb::Mask& mask,
         size_t buffer_size, uint32_t block_size,
         const std::string& compressor,
         uint8_t compression_level) :
    fn_(filename), buffer_(buffer_size), filled_(0), closed_(false),
    block_size_(block_size), compressor_(compressor), comp_level_(compression_level),
    mask_(mask)
  {
    std::cout << "dumping chunks sorted by m/z..." << std::endl;
  }

  void addPeak(const ims::Peak& peak) {
    if (filled_ == buffer_.size())
      dump();
    buffer_[filled_++] = peak;
  }

  void close() {
    dump();
    merge();
    closed_ = true;
  }

  ~Sorter() {
    if (!closed_)
      close();
  }
};

int convert_main(int argc, char** argv) {
  std::string input_filename, output_filename, compressor;
  uint32_t block_size;
  size_t buffer_size = 10000000;
  int compression_level;

  cxxopts::Options options("ims convert", " <input.imzML> <output.imzb>");
  options.add_options()
    ("block-size", "maximum number of records in a compressed block; larger values lead to slower m/z queries but smaller file size",
     cxxopts::value<uint32_t>(block_size)->default_value("4096"))
    ("compressor", "blosc compressor to be used",
     cxxopts::value<std::string>(compressor)->default_value("blosclz"))
    ("compression-level", "compression level (0-9)",
     cxxopts::value<int>(compression_level)->default_value("5"))
    ("help", "Print help");

  options.add_options("hidden")
    ("in", "", cxxopts::value<std::string>(input_filename))
    ("out", "", cxxopts::value<std::string>(output_filename));

  options.parse_positional(std::vector<std::string>{"in", "out"});

  options.parse(argc, argv);

  if (options.count("help") || input_filename.empty() || output_filename.empty()) {
    std::cout << options.help({""}) << std::endl;
    return 0;
  }

  ims::AbstractReaderPtr reader;
  if (utils::endsWith(input_filename, ".imzML"))
    reader = std::make_shared<imzml::ImzmlReader>(input_filename);
#ifdef SCILS_H5
  else if (utils::endsWith(input_filename, ".h5"))
    reader = std::make_shared<scils::H5Reader>(input_filename);
#endif
  else {
    std::cerr << "unsupported file extension" << std::endl;
    return -3;
  }
  imzb::Mask mask{reader->height(), reader->width()};

  Sorter sorter(output_filename, mask, buffer_size, block_size,
                compressor, compression_level);

  ims::Spectrum sp;
  while (reader->readNextSpectrum(sp)) {
    if (sp.coords.x >= reader->height() || sp.coords.y >= reader->width()) {
      std::cerr << "WARNING: skipped spectrum with invalid coordinates ("
                << int32_t(sp.coords.x) << ", " << int32_t(sp.coords.y) << ")" << std::endl;
      continue;
    }
    mask.set(sp.coords.x, sp.coords.y);
    for (size_t i = 0; i < sp.mzs.size(); ++i) {
      // skip invalid (NaN) and zero peaks
      if (sp.mzs[i] > 0 && sp.intensities[i] > 0)
        sorter.addPeak(ims::Peak{sp.coords, sp.mzs[i], sp.intensities[i]});
    }
  }
  sorter.close();
  return 0;
}
