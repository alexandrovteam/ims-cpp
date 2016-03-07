#pragma once

#include "ims/ims.hpp"

#include <libxml2/libxml/xmlreader.h>
#include <cstdint>
#include <map>
#include <string>
#include <vector>
#include <fstream>

namespace imzml{

class Metadata {
  static const std::map<std::string, std::string> supported_accessions_;

  std::map<std::string, std::string> dict_;
public:
  Metadata() {}

  void processNode(xmlTextReaderPtr reader);

  const std::map<std::string, std::string>& dict() const {
    return dict_;
  }
};

struct LibXmlString {
  xmlChar* str;

  LibXmlString(xmlChar* str) : str(str) {
  }

  bool operator==(const char* other) {
    return !xmlStrcmp(str, (xmlChar*)other);
  }

  bool operator==(const std::string& other) {
    return this->operator==(other.c_str());
  }

  bool operator!=(const char* other) {
    return !this->operator==(other);
  }

  bool operator!=(const std::string& other) {
    return !this->operator==(other);
  }

  operator std::string() const {
    return std::string((const char*)str);
  }

  ~LibXmlString() {
    xmlFree(str);
  }
};

LibXmlString getAttribute(xmlTextReaderPtr, const char*);

class ImzmlReader final : public ims::AbstractReader {
  Metadata metadata_;
  std::string filename_;
  std::string ibd_filename_;
  xmlTextReaderPtr xml_reader_;
  int ret;
  bool have_next_;

  std::string mz_group_name_, intensity_group_name_;

  struct ExternalArray {
    uint64_t file_offset;
    uint32_t length;
    uint32_t encoded_length;
  };

  ExternalArray mzs_;
  ExternalArray intensities_;

  LibXmlString getAttribute(const char* attribute);

  template <typename T>
  T getIntValue() {
    auto value = getAttribute("value");
    return static_cast<T>(std::atoll((const char*)value.str));
  }

  template <typename T>
  void readIntValue(T& value) {
    value = getIntValue<T>();
  }

  bool isNodeStart() const {
    return xmlTextReaderNodeType(xml_reader_) == XML_READER_TYPE_ELEMENT;
  }

  template <typename T>
  void readExternalArray(const ExternalArray& array, std::vector<T>& buffer) {
    if (array.encoded_length / array.length == 4)
      readExternalArrayHelper<float, T>(array, buffer);
    else if (array.encoded_length / array.length == 8)
      readExternalArrayHelper<double, T>(array, buffer);
  }

  template <typename T, typename U>
  void readExternalArrayHelper(const ExternalArray& array, std::vector<U>& buffer) {
    std::ifstream in(ibd_filename_);
    in.seekg(array.file_offset);
    // FIXME: endianness?
    std::vector<T> vec(array.length);
    in.read(reinterpret_cast<char*>(&vec[0]), array.encoded_length);

    buffer.resize(array.length);
    std::copy(vec.begin(), vec.end(), buffer.begin());
  }

  void readMetadata();

public:
  ImzmlReader(const std::string& filename);

  bool readNextSpectrum(ims::Spectrum& spectrum);

  const std::map<std::string, std::string>& dict() const;

  uint32_t height() const;
  uint32_t width() const;

  ~ImzmlReader();
};

}
