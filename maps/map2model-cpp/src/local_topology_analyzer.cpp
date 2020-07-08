/******************************************************************************
* A class for analyzing the local map topology (graphs).
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#include <iostream>
#include <vector>
#include <exception>
#include <fstream>

#include "graph.h"
#include "graph_utils.h"
#include "local_topology_analyzer.h"

namespace ConverterLib {

LocalTopologyAnalyzer::LocalTopologyAnalyzer(const AABB &window, int Nx, int Ny):
                                             window(window),
                                             Nx(Nx),
                                             Ny(Ny)
{
  if (Nx <= 0 || Ny <= 0) {
    throw std::runtime_error("Error: wrong topology grid dimensions!");
  }

  // Allocate the graphs grid.
  local_graphs.resize(Nx);
  for (int i = 0; i < Nx; i++) {
    local_graphs[i].resize(Ny);
  }
  // Allocate the graph coordinates grid.
  local_graph_coords.resize(Nx);
  for (int i = 0; i < Nx; i++) {
    local_graph_coords[i].resize(Ny);
  }
}

void LocalTopologyAnalyzer::BuildTopologyGrid(const std::string &file_path) const
{
  std::cout << "Building the topology grid (Jaccard indexes). Nx, Ny = " << Nx << ", " << Ny << std::endl;

  // Container to store averaged Jaccard indexes for every local map region.
  std::vector< std::vector< TopologyIndexes > > topo_indexes_avg(Nx, std::vector<TopologyIndexes>(Ny));

  std::vector<double> indexes;
  std::vector<double> topo_indexes;
  topo_indexes.resize(4);

  for (int coord_A_x = 0; coord_A_x < Nx; coord_A_x++)
    for (int coord_A_y = 0; coord_A_y < Ny; coord_A_y++)
  {
    // Initialize values.
    for (int i = 0; i < 4; i++)
      topo_indexes_avg[coord_A_x][coord_A_y].ind[i] = 0.0;

    double weight_sum = 0.;

    // Loops over neighboring cells around graph A.
    for (int coord_B_x = coord_A_x - 1; coord_B_x <= coord_A_x + 1; coord_B_x++)
      for (int coord_B_y = coord_A_y - 1; coord_B_y <= coord_A_y + 1; coord_B_y++)
    {
      // Prevent from going out of the map for subregions on the border,
      // and exclude pairs of the same cells.
      if ((coord_B_x >= 0 && coord_B_x < Nx && coord_B_y >= 0 && coord_B_y < Ny)
          && (coord_A_x != coord_B_x || coord_A_y != coord_B_y))
      {
        Graph graph1 = local_graphs[coord_A_x][coord_A_y];
        Graph graph2 = local_graphs[coord_B_x][coord_B_y];

        // Calculate full Jaccard indexes.
        CalculateJaccardIndex(graph1, graph2, indexes);

        // Store calculated Jaccard indexes.
        topo_indexes[0] = indexes[0];
        topo_indexes[1] = indexes[1];

        // Remove uncommon nodes (i.e., corresponding them graph edges) from graphs.
        // Resulting graphs should have the same set of nodes (units) but different topology.
        Graph graph1_reduced = graph1;
        Graph graph2_reduced = graph2;

        RemoveUncommonUnits(graph1_reduced, graph2);
        RemoveUncommonUnits(graph2_reduced, graph1);

        // Calculate reduced Jaccard indexes.
        CalculateJaccardIndex(graph1_reduced, graph2_reduced, indexes);

        // Store calculated Jaccard indexes.
        topo_indexes[2] = indexes[0];
        topo_indexes[3] = indexes[1];

        //------------------------------------------------------------------
        // Cell weight (according to the inverse distance between cells).
        double weight = 1. / sqrt(pow(coord_A_x - coord_B_x, 2.) + pow(coord_A_y - coord_B_y, 2.));

        // Compute the average local index.
        for (int i = 0; i < 4; i++)
        {
          topo_indexes_avg[coord_A_x][coord_A_y].ind[i] += weight * topo_indexes[i];
        }
        weight_sum = weight_sum + weight;
      }
    }

    // Normalize.
    for (int i = 0; i < 4; i++)
    {
      topo_indexes_avg[coord_A_x][coord_A_y].ind[i] /= weight_sum;
    }
  }

  //----------------------------------------------------------------------------------------
  // Calculate gradients.
  // Container that stores gradients of the averaged Jaccard indexes.
  std::vector< std::vector< TopologyIndexesGradient > > topo_indexes_gradient(Nx, std::vector<TopologyIndexesGradient>(Ny));

  if (Nx > 1 && Ny > 1)
  {
    for (int ix = 0; ix < Nx; ix++)
      for (int iy = 0; iy < Ny; iy++)
    {
      // Calculate gradients using forward difference for all indexes.
      for (int i = 0; i < 4; i++)
      {
        // Gradient in x-direction.
        if (ix != Nx - 1)
        {
          topo_indexes_gradient[ix][iy].grad[i].x
            = topo_indexes_avg[ix + 1][iy].ind[i] - topo_indexes_avg[ix][iy].ind[i];
        }
        else
        { // Gradient on the right boundary.
          topo_indexes_gradient[ix][iy].grad[i].x = topo_indexes_gradient[ix - 1][iy].grad[i].x;
        }

        // Gradient in y-direction.
        if (iy != Ny - 1)
        {
          topo_indexes_gradient[ix][iy].grad[i].x
            = topo_indexes_avg[ix][iy + 1].ind[i] - topo_indexes_avg[ix][iy].ind[i];
        }
        else
        { // Gradient on the right boundary.
          topo_indexes_gradient[ix][iy].grad[i].x = topo_indexes_gradient[ix][iy - 1].grad[i].x;
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // Write data to a file.
  std::ofstream file;
  std::string file_name = file_path + "/jaccard_indexes.txt";
  file.open(file_name.c_str());

  // Header (needed for data import in QGIS).
  file << "i j x y avg avg_weighted reduced_avg reduced_avg_weighted reduced_avg_grad reduced_avg_weighted_grad number_contacts density" << std::endl;

  for (int ix = 0; ix < Nx; ix++)
    for (int iy = 0; iy < Ny; iy++)
  {
    file << std::setprecision(14) << ix << " " << iy << " "
         << local_graph_coords[ix][iy].first << " "
         << local_graph_coords[ix][iy].second << " "
         << topo_indexes_avg[ix][iy].ind[0] << " "
         << topo_indexes_avg[ix][iy].ind[1] << " "
         << topo_indexes_avg[ix][iy].ind[2] << " "
         << topo_indexes_avg[ix][iy].ind[3] << " "
         << topo_indexes_gradient[ix][iy].grad[2].GetLength() << " "
         << topo_indexes_gradient[ix][iy].grad[3].GetLength() << " "
         << local_graphs[ix][iy].size() << " "
         << CalculateGraphDensity(local_graphs[ix][iy]) << std::endl;
  }
  file.close();
}

//=======================================================================================
void LocalTopologyAnalyzer::BuildLocalGraphs(const Parameters &par, Converter &converter)
{
  // The local map size.
  double dx = (double)(ConverterUtils::CoordToInteger(par.subregion_size_x));
  double dy = (double)(ConverterUtils::CoordToInteger(par.subregion_size_y));

  GraphWriter graph;

  std::vector<AABB> windows_local;

  for (int coord_x = 0; coord_x < Nx; coord_x++)
    for (int coord_y = 0; coord_y < Ny; coord_y++)
  {
    std::cout << "Processing a subregion: " << coord_x << " " << coord_y << std::endl;

    // Building a local clipping window.
    AABB window_local(0);
    window_local.TopLeft.X  = (cInt)(window.TopLeft.X + (coord_x) * dx);
    window_local.TopLeft.Y  = (cInt)(window.TopLeft.Y + (coord_y) * dy);
    window_local.BtmRight.X = (cInt)(window.TopLeft.X + (coord_x + 1) * dx);
    window_local.BtmRight.Y = (cInt)(window.TopLeft.Y + (coord_y + 1) * dy);

    // A postfix for filenames.
    std::string postfix = SSTR(coord_x) + "_" + SSTR(coord_y);

    // Clipping data to a local window.
    converter.ClipData(window_local);

    windows_local.push_back(window_local);

    // An image showing a local clipping window.
    //converter.WriteImage("polygons_local_" + postfix, window, 0, NULL);

    converter.BuildUnitsAndGroupsLists();

    Contacts contacts;
    converter.FindContacts(contacts, par.intersectPolygonsDistanceBuffer);

    converter.IdentifyPolygonContactTypes(contacts, par.angleEpsilon, par.distanceEpsilon);
    converter.SplitMixedPolygonContacts(contacts);

    // Building unit contacts.
    UnitContacts unit_contacts;
    converter.BuildUnitContactsList(contacts, unit_contacts);

    // Store the local graph (for further analysis).
    local_graphs[coord_x][coord_y] = unit_contacts;

    // Store the local graph world coordinates (of the center of the map subregion).
    double world_x = ConverterUtils::CoordToDouble(window_local.TopLeft.X + (cInt)(dx / 2.));
    double world_y = ConverterUtils::CoordToDouble(window_local.TopLeft.Y + (cInt)(dy / 2.));
    local_graph_coords[coord_x][coord_y] = std::pair<double, double>(world_x, world_y);

    std::string depositName = "";

    // Generating a complete graph.
    std::string file_graph = par.path_output + "/graph_" + postfix + ".gml";
    graph.WriteGraph(file_graph, unit_contacts, "UNIT_NAME", false, 1,
                     par.graph_edge_width_categories, par.graph_edge_direction_type,
                     depositName);
  }

  // An image showing polygons with topology grid.
  converter.ClipData(window);
  converter.WriteImage("polygons_with_grid", windows_local, 0, NULL);
}

}

