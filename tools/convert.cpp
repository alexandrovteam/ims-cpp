#include "utils/convert.hpp"

#include "cxxopts.hpp"

#include <string>
#include <vector>

int convert_main(int argc, char** argv) {
  std::string input_filename, output_filename;
  size_t buffer_size = 10000000;

  int compression_level;
  imzb::ImzbCompressionSettings compression_settings;

  cxxopts::Options options("ims convert", " <input.imzML> <output.imzb>");
  options.add_options()("block-size",
      "maximum number of records in a compressed block; larger "
      "values lead to slower m/z queries but smaller file size",
      cxxopts::value<uint32_t>(compression_settings.block_size)->default_value("4096"))(
      "compressor", "blosc compressor to be used",
      cxxopts::value<std::string>(compression_settings.compressor)
          ->default_value("blosclz"))("compression-level", "compression level (0-9)",
      cxxopts::value<int>(compression_level)->default_value("5"))("help", "Print help");

  options.add_options("hidden")("in", "", cxxopts::value<std::string>(input_filename))(
      "out", "", cxxopts::value<std::string>(output_filename));

  options.parse_positional(std::vector<std::string>{"in", "out"});

  options.parse(argc, argv);
  compression_settings.compression_level = compression_level;

  if (options.count("help") || input_filename.empty() || output_filename.empty()) {
    std::cout << options.help({""}) << std::endl;
    return 0;
  }

  try {
    convertToImzb(input_filename, output_filename, buffer_size, compression_settings);
  } catch (std::runtime_error& e) {
    std::cerr << e.what() << std::endl;
    return -3;
  }
}
