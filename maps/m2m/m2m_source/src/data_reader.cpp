/******************************************************************************
 * A class for reading input data.
 *
 * Author: Vitaliy Ogarko, vogarko@gmail.com
 *******************************************************************************/

#include <fstream>
#include <iostream>

#include "clipper.hpp"
#include "converter_types.h"
#include "data_reader.h"
#include "parameters_utils.h"

namespace ConverterLib {

// Reads coordinates of an object from a string in WKT format.
static Paths ReadCoordinates(const std::string &__line,
                             const std::string keyword) {
  std::string line = __line;

  std::size_t pos = line.find(keyword);
  if (pos == std::string::npos) {
    std::cout << "Error: not found keyword " << keyword << "\n"
              << "line =" << line << std::endl;
    exit(0);
  }

  Paths paths;

  // Extract line part followed by "keyword".
  line = line.substr(pos + keyword.length());

  Path path;

  while (true) {
    // Remove leading brackets and spaces.
    while (line[0] == '(' || line[0] == ' ') {
      line.erase(0, 1);
    }

    std::istringstream iss(line);
    double x, y;
    cInt ix, iy;

    // Reading next pair of coordinates (separated by a space).
    if (!(iss >> x >> y)) {
      std::cout << "Error in reading a pair of coordinates!" << std::endl;
      break;
    }

    // Converting coordinates to integers (since the Clipper lib works with
    // integers).
    ix = ConverterUtils::CoordToInteger(x);
    iy = ConverterUtils::CoordToInteger(y);

    path.push_back(IntPoint(ix, iy));

    // Getting positions of the following comma and closing bracket.
    std::size_t pos_comma = line.find(',');
    std::size_t pos_close_bracket = line.find(')');

    if (pos_comma < pos_close_bracket) {
      // There is one more pair of coordinates (separated by a comma).
      // Remove from the line a processed part.
      line.erase(0, pos_comma + 1);
    } else {
      // End of polygon, i.e., found closing bracket ")".
      // Adding next path.
      paths.push_back(path);
      path.clear();

      // Search if there are more polygons at this line.
      std::size_t pos_open_bracket = line.find('(');
      if (pos_open_bracket != std::string::npos) {
        // Move to the beginning of the next polygon.
        line.erase(0, pos_open_bracket + 1);

        // Will be adding a new path for this object.
      } else {
        // No more paths.
        break;
      }
    }
  }
  return paths;
}

static std::string RemoveQuotations(const std::string &line) {
  std::string res = line;
  res.erase(std::remove(res.begin(), res.end(), '\"'), res.end());
  return res;
}

static int String2Int(const std::string str) {
  int res;
  std::istringstream(str) >> res;
  return res;
}

static double String2Double(const std::string str) {
  double res;
  std::istringstream(str) >> res;
  return res;
}

std::string ReadField(const std::string &__line, size_t column,
                      const char delimiter) {
  std::string line = __line;

  for (size_t i = 1; i < column; i++) {
    std::size_t pos = line.find(delimiter);
    line = line.substr(pos + 1);
  }

  std::size_t pos = line.find(delimiter);
  std::string field = line.substr(0, pos);
  // Removing quotes if any.
  field = RemoveQuotations(field);

  return field;
}

std::map<std::string, size_t> ReadHeader(const std::string &__line,
                                         const char delimiter) {
  std::map<std::string, size_t> headerFields;
  std::string line = __line;

  std::size_t pos = line.find(delimiter);
  size_t column = 1;
  while (pos != std::string::npos) {
    std::string field = line.substr(0, pos);
    line = line.substr(pos + 1);

    // Store the field.
    headerFields[field] = column;

    column++;
    pos = line.find(delimiter);
  }
  // Store the last field.
  std::string field = line;
  headerFields[field] = column;

  return headerFields;
}

int ReadIntField(const std::string &__line, size_t column,
                 const char delimiter) {
  return String2Int(ReadField(__line, column, delimiter));
}

// Trim a string from left.
static inline std::string ltrim(std::string s, const char *t = " \t\n\r\f\v") {
  s.erase(0, s.find_first_not_of(t));
  return s;
}

// Trim a string from right.
static inline std::string rtrim(std::string s, const char *t = " \t\n\r\f\v") {
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

// Trim a string from left and right.
static inline std::string trim(std::string s, const char *t = " \t\n\r\f\v") {
  return ltrim(rtrim(s, t), t);
}

static size_t getFieldColumn(const std::map<std::string, size_t> &headerFields,
                             const std::string &fieldName) {
  const std::map<std::string, size_t>::const_iterator it =
      headerFields.find(fieldName);
  if (it != headerFields.end()) {
    return it->second;
  } else {
    std::cout << "Unknown field name: " << fieldName << std::endl;
    exit(0);
  }
}

int ReadDataObj(const std::string &filename, const std::string &keyword,
                const std::map<std::string, std::string> &constNames,
                Objects &objects, const std::vector<int> &idsToRead,
                const std::vector<std::string> &idsToReadString) {
  std::ifstream myfile(filename.c_str());

  if (!myfile.good()) {
    std::cout << "Error in opening the file " << filename << "\n";
    exit(0);
  } else {
    std::cout << "Reading data from the file " << filename << "\n";
  }

  std::string line;

  // Reading the header.
  std::getline(myfile, line);
  std::map<std::string, size_t> headerFields = ReadHeader(line);

  // Loop over the file lines.
  while (std::getline(myfile, line).good()) {
    // Skip lines without the keyword.
    std::size_t pos = line.find(keyword);
    if (pos == std::string::npos)
      continue;

    std::string fieldName;

    // Extracting ID.
    int id;
    if (keyword == "LINESTRING") {
      fieldName = getConstName(constNames, "FIELD_FAULT_ID");
      id = ReadIntField(line, getFieldColumn(headerFields, fieldName));

    } else if (keyword == "POLYGON") {
      fieldName = getConstName(constNames, "FIELD_POLYGON_ID");
      id = ReadIntField(line, getFieldColumn(headerFields, fieldName));

    } else if (keyword == "POINT") {
      // Currently deposits have not integer id, so assign id here.
      id = objects.size();
    } else {
      std::cout << "Unknown object type: " << keyword << "\n";
      exit(0);
    }

    // Processing only IDs form the input list, if it is not empty.
    if (!idsToRead.empty() &&
        std::find(idsToRead.begin(), idsToRead.end(), id) == idsToRead.end()) {
      continue;
    }

    if (keyword == "LINESTRING") {
      // Ignore faults with features of this type.
      fieldName = getConstName(constNames, "FIELD_FAULT_FEATURE");
      std::string feature =
          trim(ReadField(line, getFieldColumn(headerFields, fieldName)));

      if (feature == getConstName(constNames, "FAULT_AXIAL_FEATURE_NAME")) {
        continue;
      }
    }

    if (keyword == "POINT") {
      // Ignore deposits of this site type.
      fieldName = getConstName(constNames, "FIELD_SITE_TYPE");
      std::string site_type =
          trim(ReadField(line, getFieldColumn(headerFields, fieldName)));

      if (site_type == getConstName(constNames, "IGNORE_DEPOSITS_SITE_TYPE")) {
        continue;
      }
    }

    Object object;
    object.id = id;

    if (keyword == "POINT") {
      fieldName = getConstName(constNames, "FIELD_SITE_CODE");
      object.site_code =
          trim(ReadField(line, getFieldColumn(headerFields, fieldName)));

      fieldName = getConstName(constNames, "FIELD_SITE_COMMO");
      object.site_commo =
          trim(ReadField(line, getFieldColumn(headerFields, fieldName)));

      // Processing only codes form the input list, if it is not empty.
      if (!idsToReadString.empty() &&
          std::find(idsToReadString.begin(), idsToReadString.end(),
                    object.site_code) == idsToReadString.end()) {
        continue;
      }

    } else if (keyword == "POLYGON") {
      fieldName = getConstName(constNames, "FIELD_POLYGON_LEVEL1_NAME");
      object.name =
          trim(ReadField(line, getFieldColumn(headerFields, fieldName)));

      fieldName = getConstName(constNames, "FIELD_POLYGON_LEVEL2_NAME");
      object.group =
          trim(ReadField(line, getFieldColumn(headerFields, fieldName)));

      fieldName = getConstName(constNames, "FIELD_POLYGON_MIN_AGE");
      object.min_age = String2Double(
          ReadField(line, getFieldColumn(headerFields, fieldName)));

      fieldName = getConstName(constNames, "FIELD_POLYGON_MAX_AGE");
      object.max_age = String2Double(
          ReadField(line, getFieldColumn(headerFields, fieldName)));

      fieldName = getConstName(constNames, "FIELD_POLYGON_CODE");
      object.code = ReadField(line, getFieldColumn(headerFields, fieldName));

      fieldName = getConstName(constNames, "FIELD_POLYGON_DESCRIPTION");
      object.description =
          ReadField(line, getFieldColumn(headerFields, fieldName));

      fieldName = getConstName(constNames, "FIELD_POLYGON_ROCKTYPE1");
      object.rocktype1 =
          ReadField(line, getFieldColumn(headerFields, fieldName));

      fieldName = getConstName(constNames, "FIELD_POLYGON_ROCKTYPE2");
      object.rocktype2 =
          ReadField(line, getFieldColumn(headerFields, fieldName));
    }

    // Reading coordinates.
    fieldName = getConstName(constNames, "FIELD_COORDINATES");
    line = ReadField(line, getFieldColumn(headerFields, fieldName));
    object.paths = ReadCoordinates(line, keyword);

    // Adding next object.
    objects.push_back(object);
  }
  myfile.close();

  std::cout << "Objects read: " << objects.size() << std::endl;

  return (int)(objects.size());
}
//========================================================================================

}; // namespace ConverterLib
