/******************************************************************************
* Parameter utils.
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#include <stdexcept>

#include "parameters_utils.h"

namespace ConverterLib {

std::string getConstName(const std::map<std::string, std::string>& constNames,
                                const std::string& key) {
    const std::map<std::string, std::string>::const_iterator it = constNames.find(key);
    if (it != constNames.end()) {
        return it->second;
    } else {
        throw std::runtime_error("Unknown constant name key: " + key);
    }
}

}


