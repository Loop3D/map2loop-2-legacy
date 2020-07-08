/******************************************************************************
 * Parameter reader.
 *
 * Author: Vitaliy Ogarko, vogarko@gmail.com
 *******************************************************************************/

#ifndef parameters_reader_h
#define parameters_reader_h

#include <map>
#include <string>
#include <vector>

namespace ConverterLib {

class Parameters {
public:
  // Field names in CSV data files, and some other literal constants.
  std::map<std::string, std::string> constNames;

  // Path to the output data folder.
  std::string path_output;

  // Converter lib input parameters.
  std::string path_geology;
  std::string path_faults;
  std::string path_points;

  double clipping_window[4];

  std::vector<double> graph_edge_width_categories;
  int graph_edge_direction_type;

  int partial_graph_polygon_id;
  int partial_graph_depth;

  // Size of the map subregion size for calculating Jaccard indexes.
  double subregion_size_x;
  double subregion_size_y;

  // Minimum length fraction in the mixed contact for the strat/faults graphs.
  double minFractionInMixedContact;

  // Constants for IntersectContactWithFault routine.
  double angleEpsilon;
  double distanceEpsilon;
  // Distance buffer from a point to segment (to find if a fault stops on //
  // another fault).
  double faultFaultDistanceBuffer;

  // Distance buffer from a point (deposit) to contact (to find if a point is
  // sitting on the contact).
  double pointToContactDistanceBuffer;

  // Distance buffer for IntersectPolygons (for finding contacts) - for bad
  // maps.
  double intersectPolygonsDistanceBuffer;

  // Deposit names for which we are adding info on the graph.
  std::vector<std::string> depositListForGraphInfo;

  // Reading input parameters from a file.
  // Returns parameter lines that were read.
  std::string Read(const std::string &filename);
  std::string directRead(const std::string output, const std::string geology,
                         const std::string faults, const std::string points);
};

std::string getConstName(const std::map<std::string, std::string> &constNames,
                         const std::string &key);

} // namespace ConverterLib

#endif
