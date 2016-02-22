#include "imzb/writer.hpp"

#include "blosc.h"

#include <ios>
#include <cassert>

imzb::ImzbWriter::ImzbWriter(const std::string& filename) :
		out_{filename, std::ios::binary},
		in_buf_(BLOCK_SIZE),
		out_buf_(BLOCK_SIZE * sizeof(ims::Peak) + BLOSC_MAX_OVERHEAD * 2),
		filled_(0), bytes_written_(0), filename_(filename)
{
}

void imzb::ImzbWriter::dump()
{
	if (filled_ == 0)
		return;
	int n = blosc_compress_ctx(5, 1, sizeof(ims::Peak),
														 filled_ * sizeof(ims::Peak), &in_buf_[0],
														 &out_buf_[0], out_buf_.size(),
														 "blosclz", 0, 1);
	assert(n > 0);
	filled_ = 0;
	out_.write(&out_buf_[0], n);
	index_.mzs.push_back(in_buf_[0].mz);
	index_.offsets.push_back(bytes_written_);
	bytes_written_ += n;
}

void imzb::ImzbWriter::writePeak(const ims::Peak& peak)
{
	if (filled_ == in_buf_.size())
		dump();
	in_buf_[filled_++] = peak;
}

void imzb::ImzbWriter::close()
{
	dump();
	out_.close();

	std::ofstream out_idx(filename_ + ".idx", std::ios::binary);
	index_.write(out_idx);
	out_idx.close();
}

imzb::ImzbWriter::~ImzbWriter()
{
	if (out_.is_open())
		close();
}


