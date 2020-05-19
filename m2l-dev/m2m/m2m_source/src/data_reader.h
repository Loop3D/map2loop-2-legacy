/******************************************************************************
* A class for reading input data.
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#ifndef data_reader_h
#define data_reader_h

#include <string>
#include <map>

#include "converter_types.h"

namespace ConverterLib {

// Reading objects (multi-polygons or faults) from file.
// Returns the number of objects read.
int ReadDataObj(const std::string &filename, const std::string& keyword,
                const std::map<std::string, std::string>& constNames,
                Objects &objects,
                const std::vector<int>& idsToRead,
                const std::vector<std::string>& idsToReadString = std::vector<std::string>());

// Reads field value in a given column (removes quotes if any).
std::string ReadField(const std::string &__line, size_t column, const char delimiter = '\t');

// Reads header field names.
std::map<std::string, size_t> ReadHeader(const std::string &__line, const char delimiter = '\t');

// Reads integer value in a given column (removes quotes if any).
int ReadIntField(const std::string &__line, size_t column, const char delimiter = '\t');

};

#endif
