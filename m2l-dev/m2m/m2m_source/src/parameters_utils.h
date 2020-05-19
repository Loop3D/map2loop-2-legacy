/******************************************************************************
* Parameter utils.
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#ifndef parameters_utils_h
#define parameters_utils_h

#include <map>
#include <string>

namespace ConverterLib {

std::string getConstName(const std::map<std::string, std::string>& constNames,
                         const std::string& key);

}

#endif

