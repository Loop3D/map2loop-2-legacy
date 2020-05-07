/******************************************************************************
* A class for analyzing the local map topology (graphs).
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#ifndef local_topology_analyzer_h
#define local_topology_analyzer_h

#include "converter.h"
#include "parameters_reader.h"

namespace ConverterLib {

typedef UnitContacts Graph;
// 2D container of local graphs (on a grid).
typedef std::vector< std::vector< Graph> > Graphs;
// 2D container of local graph world coordinates.
typedef std::vector< std::vector < std::pair <double, double> > > GraphCoords;

class TopologyIndexes
{
public:
  double ind[4];
};

class Vec2D
{
public:
  double x, y;
  double GetLength() { return sqrt(x * x + y * y); }
};

class TopologyIndexesGradient
{
public:
  Vec2D grad[4];
};

class LocalTopologyAnalyzer
{
public:
  // Constructor.
  LocalTopologyAnalyzer(const AABB &window, int Nx, int Ny);

  //! Build graphs for subregions (clipped rectangles) of the map.
  void BuildLocalGraphs(const Parameters &par, Converter &converter);

  //! Calculate topology indexes (Jaccard index and its modifications) for all local map regions.
  //! Store the output data in files.
  void BuildTopologyGrid(const std::string &file_path) const;

private:
  //! Clipping window.
  AABB window;
  //! Topology grid (of local graphs) size.
  int Nx, Ny;
  //! Container for storing local graphs.
  Graphs local_graphs;
  //! Store local graph world coordinates (of the map subregion centers).
  GraphCoords local_graph_coords;
};

}

#endif
