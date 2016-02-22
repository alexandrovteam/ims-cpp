#include "cffi/common.hpp"
#include "imzb/reader.hpp"

using namespace cffi;

extern "C" {
	DLL_PUBLIC imzb::ImzbReader* imzb_reader_new(const char* filename) {
		return wrap_catch<imzb::ImzbReader*>(nullptr, [&]() {
			return new imzb::ImzbReader(filename);
		});
	}

	DLL_PUBLIC void imzb_reader_free(imzb::ImzbReader* reader) {
		delete reader;
	}

	DLL_PUBLIC int imzb_reader_height(imzb::ImzbReader* reader) {
		return reader->height();
	}

	DLL_PUBLIC int imzb_reader_width(imzb::ImzbReader* reader) {
		return reader->width();
	}

	DLL_PUBLIC int imzb_reader_image(imzb::ImzbReader* reader, double mz, double ppm, float* outbuf) {
		return wrap_catch<int>(-1, [&]() -> int {
				reader->readImage(mz, ppm, outbuf);
				return 0;
		});
	}
}
