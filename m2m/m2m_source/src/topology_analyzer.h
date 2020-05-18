/******************************************************************************
* A class for analyzing the map topology (graphs).
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#ifndef topology_analyzer_h
#define topology_analyzer_h

#include "converter.h"
#include "intersector.h"

namespace ConverterLib {

class TopologyAnalyzer
{
public:
  // Constructor.
  TopologyAnalyzer(AABB window = AABB(0));

  void Initialize(const Parameters &par);

  void AnalyzeLocalTopology(const Parameters &par);
  void AnalyzeGlobalTopology(const Parameters &par, const std::string &graph_comments);

private:
  //! Object for generating topology graphs.
  Converter converter;
  //! Clipping window.
  AABB window;

  // Writes unit vs fault intersection list.
  static void WriteUnitFaultIntersectionList(const std::string &fileName,
                const UnitFaultIntersectionList &unitFaultIntersectionList);

  // Writes fault vs fault intersection list.
  static void WriteFaultIntersectionList(const std::string &fileName,
                const FaultIntersectionList &faultIntersectionList);

  // Writes points vs polygons & contacts topology table.
  static void WritePointPolygonIntersectionList(const std::string &fileName,
          const PointPolygonIntersectionList &pointPolygonIntersectionList);
};

}

#endif
