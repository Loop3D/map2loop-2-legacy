/******************************************************************************
* A class for gml graphs generation.
*
* Author: Vitaliy Ogarko (2016).
*******************************************************************************/

#ifndef graph_h
#define graph_h

#include <string>

#include "converter_types.h"

namespace ConverterLib {

class GraphWriter
{
public:
  //! Writes the graph of units connections into a file in GML format.
  void WriteGraph(const std::string &file_graph, const UnitContacts &contacts,
                  const std::string &label_type, bool exclude_sills,
                  int edge_width_type, const std::vector<double> &edge_width_categories,
                  int edge_direction_type,
                  const std::string& depositName,
                  const std::string &comments = "");

private:

};

}

#endif
