#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "imzml/reader.hpp"

namespace imzml {

LibXmlString getAttribute(xmlTextReaderPtr reader, const char* attribute) {
	return xmlTextReaderGetAttribute(reader, (xmlChar*)attribute);
}

void Metadata::processNode(xmlTextReaderPtr reader) {
	auto accession = getAttribute(reader, "accession");
	if (accession == nullptr)
		return;

	auto it = supported_accessions_.find(accession);
	if (it == supported_accessions_.end())
		return;

	auto value = getAttribute(reader, "value");
	dict_[it->second] = value;
}

const std::map<std::string, std::string> Metadata::supported_accessions_ =
{
	{ "IMS:1000042", "max count of pixels x" },
	{ "IMS:1000043", "max count of pixels y" },
	{ "IMS:1000044", "max dimension x" },
	{ "IMS:1000045", "max dimension y" },
	{ "IMS:1000046", "pixel size x" },
	{ "IMS:1000047", "pixel size y" },
	{ "IMS:1000835", "matrix solution concentration" },
	{  "MS:1000843", "wavelength" },
	{  "MS:1000844", "focus diameter x" },
	{  "MS:1000845", "focus diameter y" },
	{  "MS:1000846", "pulse energy" },
	{  "MS:1000847", "pulse duration" },
	{  "MS:1000848", "attenuation" }
};

void ImzmlReader::readMetadata() {
	std::string current_param_group;

	ret = xmlTextReaderRead(xml_reader_);
	while (ret == 1) {
		LibXmlString name = xmlTextReaderName(xml_reader_);
		if (name == "cvParam") {
			metadata_.processNode(xml_reader_);

			auto accession = getAttribute("accession");
			if (accession == "MS:1000514")
				mz_group_name_ = current_param_group;
			else if (accession == "MS:1000515")
				intensity_group_name_ = current_param_group;
		} else if (name == "spectrum" && isNodeStart()) {
			have_next_ = true;
			break;
		} else if (name == "referenceableParamGroup") {
			current_param_group = getAttribute("id");
		}

		ret = xmlTextReaderRead(xml_reader_);
	}
}

LibXmlString ImzmlReader::getAttribute(const char* attribute) {
	return imzml::getAttribute(xml_reader_, attribute);
}

ImzmlReader::ImzmlReader(const std::string& filename) : filename_{filename}, have_next_{false} {
	xml_reader_ = xmlNewTextReaderFilename(filename.c_str());
	ibd_filename_ = filename_.substr(0, filename_.length() - 5) + "ibd";
	readMetadata();
}

bool ImzmlReader::readNextSpectrum(ims::Spectrum& spectrum) {
	if (!have_next_)
		return false;

	have_next_ = false;

	spectrum.coords = ims::Position{0, 0, 0};
	spectrum.mzs.clear();
	spectrum.intensities.clear();

	ExternalArray* array;

	ret = xmlTextReaderRead(xml_reader_);
	while (ret == 1) {
		LibXmlString name = xmlTextReaderName(xml_reader_);

		if (name == "cvParam") {
			auto accession = getAttribute("accession");
			if      (accession == "IMS:1000050") readIntValue(spectrum.coords.x);
			else if (accession == "IMS:1000051") readIntValue(spectrum.coords.y);
			else if (accession == "IMS:1000052") readIntValue(spectrum.coords.z);
			else if (accession == "IMS:1000102") readIntValue(array->file_offset);
			else if (accession == "IMS:1000103") readIntValue(array->length);
			else if (accession == "IMS:1000104") readIntValue(array->encoded_length);
		} else if (name == "spectrum" && isNodeStart()) {
			have_next_ = true;
			break;
		} else if (name == "referenceableParamGroupRef") {
			auto ref = getAttribute("ref");
			if (ref == mz_group_name_)
				array = &mzs_;
			else if (ref == intensity_group_name_)
				array = &intensities_;
		}
		ret = xmlTextReaderRead(xml_reader_);
	}

	readExternalArray<double>(mzs_, spectrum.mzs);
	readExternalArray<float>(intensities_, spectrum.intensities);

	return true;
}

ImzmlReader::~ImzmlReader() {
	xmlFreeTextReader(xml_reader_);
}

const std::map<std::string, std::string>& ImzmlReader::dict() const {
	return metadata_.dict();
}

uint32_t ImzmlReader::height() const {
	return 1 + atoi(dict().find("max count of pixels x")->second.c_str());
}

uint32_t ImzmlReader::width() const {
	return 1 + atoi(dict().find("max count of pixels y")->second.c_str());
}

}
