/******************************************************************************
* A class for analyzing the map topology (graphs).
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#include <algorithm>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cassert>

#include "graph.h"
#include "local_topology_analyzer.h"
#include "parameters_utils.h"

#include "topology_analyzer.h"

namespace ConverterLib {

TopologyAnalyzer::TopologyAnalyzer(AABB window):
    window(window)
{
}
//=======================================================================================

void TopologyAnalyzer::Initialize(const Parameters &par)
{
  // Objects' IDs to read. If empty, then read all objects.
  std::vector<int> faultsIdsToRead, polygonsIdsToRead, pointsIdsToRead;
  std::vector<std::string> pointsIdsToReadString;

  bool test = false;

  if (test)
  {
    // Some test cases with subset of real data.
    //faultsIdsToRead.push_back(69);
    //faultsIdsToRead.push_back(685);

    //polygonsIdsToRead.push_back(2);
    //polygonsIdsToRead.push_back(920);

    pointsIdsToReadString.push_back("S0225657");
    pointsIdsToReadString.push_back("S0021264");
  }

  // Set output data folder.
  converter.SetOutputFolder(par.path_output);

  // Setting global variables.
  converter.SILL_DESCRIPTION = getConstName(par.constNames, "SILL_STRING");
  converter.IGNEOUS_STRING = getConstName(par.constNames, "IGNEOUS_STRING");
  converter.VOLCANIC_STRING = getConstName(par.constNames, "VOLCANIC_STRING");

  // Reading the map data.
  converter.ReadData(par.path_geology, "POLYGON", par.constNames, polygonsIdsToRead);
  converter.ReadData(par.path_faults, "LINESTRING", par.constNames, faultsIdsToRead);
  if (par.path_points != "") {
      converter.ReadData(par.path_points, "POINT", par.constNames, pointsIdsToRead, pointsIdsToReadString);
  }

  // Initialize the clipping window.
  if (par.clipping_window[0] == par.clipping_window[2]
      || par.clipping_window[1] == par.clipping_window[3])
  {// Zero clipping window.
    // Set clipping window to the total map size.
    window = converter.GetAllPolygonsAABB();

    std::cout << std::setprecision(14)
              << "Zero input clipping window. Set the window to the total map size: "
              << window.TopLeft << " " << window.BtmRight << std::endl;
  }
  else
  {
    window = AABB(par.clipping_window[0], par.clipping_window[1], par.clipping_window[2], par.clipping_window[3]);
  }

  std::vector<AABB> windows;
  windows.push_back(window);

  // Write some images with data read.
  converter.WriteImage("polygons_read", windows, 0, NULL);
  converter.WriteImage("polygons_read_with_points", windows, 1, NULL);
  converter.WriteImage("polygons_read_with_faults", windows, 2, NULL);

  // Clip data according to input clipping window.
  converter.ClipData(window);

  converter.WriteImage("polygons_clipped", windows, 0, NULL);
  converter.WriteImage("polygons_clipped_with_points", windows, 1, NULL);
  converter.WriteImage("polygons_clipped_with_faults", windows, 2, NULL);
}
//=======================================================================================

void TopologyAnalyzer::AnalyzeLocalTopology(const Parameters &par)
{
  // Topology grid size (for local graphs).
  int Nx = (int)((window.BtmRight.X - window.TopLeft.X) / ConverterUtils::CoordToInteger(par.subregion_size_x));
  int Ny = (int)((window.BtmRight.Y - window.TopLeft.Y) / ConverterUtils::CoordToInteger(par.subregion_size_y));

  std::cout << "Topology grid size Nx, Ny = " << Nx << " " << Ny << std::endl;

  LocalTopologyAnalyzer locTopoAnalyzer(window, Nx, Ny);

  // Build local graphs (for each map subregion).
  locTopoAnalyzer.BuildLocalGraphs(par, converter);

  // Calculate (modified) Jaccard indexes.
  locTopoAnalyzer.BuildTopologyGrid(par.path_output);
}
//=======================================================================================================

void TopologyAnalyzer::WriteUnitFaultIntersectionList(const std::string &fileName,
              const UnitFaultIntersectionList &unitFaultIntersectionList)
{
  std::ofstream file;
  file.open(fileName.c_str());

  size_t counter = 0;
  for (UnitFaultIntersectionList::const_iterator it = unitFaultIntersectionList.begin();
      it != unitFaultIntersectionList.end(); ++it) {

    std::string unitName = it->first;
    file << counter << ", " << unitName << ", {";

    std::list<int> faultIdList = it->second;
    bool first = true;

    for (std::list<int>::const_iterator faultIdIt = faultIdList.begin();
         faultIdIt != faultIdList.end(); ++faultIdIt) {

      if (!first) {
        file << ", ";
      } else {
        first = false;
      }
      file << *faultIdIt;
    }
    file << "}" << std::endl;
    counter++;
  }
  file.close();
  std::cout << "Wrote " << counter << " items to unit-fault intersection list: " << fileName << std::endl;
}
//=======================================================================================

void TopologyAnalyzer::WriteFaultIntersectionList(const std::string &fileName,
              const FaultIntersectionList &faultIntersectionList)
{
  std::ofstream file;
  file.open(fileName.c_str());

  file << std::fixed << std::setprecision(1);

  size_t counter = 0;
  for (FaultIntersectionList::const_iterator it = faultIntersectionList.begin();
      it != faultIntersectionList.end(); ++it) {

    int faultId1 = it->first;
    file << counter << ", " << faultId1 << ", {";

    std::list<FaultIntersectionData> faultDataList = it->second;
    bool first = true;

    for (std::list<FaultIntersectionData>::const_iterator faultDataIt = faultDataList.begin();
        faultDataIt != faultDataList.end(); ++faultDataIt) {

      if (!first) {
        file << ", ";
      } else {
        first = false;
      }
      file << "(" << faultDataIt->first << ", " << faultDataIt->second.first << ", " << faultDataIt->second.second << ")";
    }
    file << "}" << std::endl;
    counter++;
  }
  file.close();
  std::cout << "Wrote " << counter << " items to fault-fault intersection list: " << fileName << std::endl;
}
//=======================================================================================

void TopologyAnalyzer::WritePointPolygonIntersectionList(const std::string &fileName,
        const PointPolygonIntersectionList &pointPolygonIntersectionList)
{
    std::ofstream file;
    file.open(fileName.c_str());

    // Adding a title.
    file << std::setprecision(16) <<
            "#ID,SiteCode,SiteCommo," <<
            "PositionX,PositionY,DistanceToContact," <<
            "ContainingPolygonID,NeighbourPolygonID," <<
            "ContainingPolygonCode,NeighbourPolygonCode," <<
            "ContactType,FaultID" << std::endl;

    size_t counter = 0;
    for (PointPolygonIntersectionList::const_iterator it = pointPolygonIntersectionList.begin();
        it != pointPolygonIntersectionList.end(); ++it) {

        std::string contactType;
        int faultID = -1;

        if (it->contactType == StratigraphicContact) {
            contactType = "Strat";
        }
        else if (it->contactType == FaultContact) {
            contactType = "Fault";
            faultID = it->faultObjectID;
        }
        else if (it->contactType == IntrusiveIgneousContact) {
            contactType = "IntrIgn";
        }
        else if (it->contactType == IgneousContact) {
            contactType = "Igneous";
        }
        else {
            assert(false);
        }

        double posX = ConverterUtils::CoordToDouble(it->point->paths[0][0].X);
        double posY = ConverterUtils::CoordToDouble(it->point->paths[0][0].Y);

        // Process site commodity.
        std::string commodity = it->point->site_commo;
        // replace all ',' to '_'
        std::replace(commodity.begin(), commodity.end(), ',', '_');
        // remove spaces
        std::string::iterator end_pos = std::remove(commodity.begin(), commodity.end(), ' ');
        commodity.erase(end_pos, commodity.end());

        // Write a row.
        file << std::setprecision(16) << std::defaultfloat <<
                it->point->id << "," <<
                it->point->site_code << "," <<
                commodity << "," <<
                posX << "," <<
                posY << "," <<
                std::setprecision(2) << std::fixed <<
                sqrt(it->distanceToContact2) << "," <<
                it->containingPolygon->id << "," <<
                it->neighbourPolygon->id << "," <<
                it->containingPolygon->code << "," <<
                it->neighbourPolygon->code << "," <<
                contactType << "," <<
                faultID << std::endl;

        counter++;
    }

    file.close();
    std::cout << "Wrote " << counter << " items to point polygon intersection list: " << fileName << std::endl;
}

void TopologyAnalyzer::AnalyzeGlobalTopology(const Parameters &par, const std::string &graph_comments)
{
  converter.BuildUnitsAndGroupsLists();

  //------------------------------------------------------------------
  {
    // Build and write to a file a unit-fault intersection list.
    UnitFaultIntersectionList unitFaultIntersectionList;
    IntersectFaultsWithPolygons(converter.GetFaults(), converter.GetPolygons(), unitFaultIntersectionList);

    std::string fileName = par.path_output + "/unit-fault-intersection.txt";
    WriteUnitFaultIntersectionList(fileName, unitFaultIntersectionList);
  }

  //------------------------------------------------------------------
  {
    // Build and write to a file a fault-fault intersection list.
    FaultIntersectionList faultIntersectionList;
    IntersectFaultsWithFaults(converter.GetFaults(), faultIntersectionList, par.faultFaultDistanceBuffer);

    std::string fileName = par.path_output + "/fault-fault-intersection.txt";
    WriteFaultIntersectionList(fileName, faultIntersectionList);
  }

  //-----------------------------------------------------------------------------------
  std::vector<AABB> windows;
  windows.push_back(window);

  Contacts contacts, fault_contacts, stratigraphic_contacts;

  double contact_length = converter.FindContacts(contacts);
  converter.WriteImage("polygons_with_all_contacts", windows, 0, &contacts);

  converter.IdentifyPolygonContactTypes(contacts, par.angleEpsilon, par.distanceEpsilon);
  converter.SplitMixedPolygonContacts(contacts);

  for (Contacts::const_iterator it = contacts.begin(); it != contacts.end(); ++it)
  {
    if (it->type == FaultContact) fault_contacts.push_back(*it);
    if (it->type == StratigraphicContact) stratigraphic_contacts.push_back(*it);
  }
  converter.WriteImage("polygons_with_fault_contacts", windows, 0, &fault_contacts);
  converter.WriteImage("polygons_with_stratigraphic_contacts", windows, 0, &stratigraphic_contacts);

  // For testing.
  converter.PrintLengths(contact_length);

  //-----------------------------------------------------------------------------------
  // Build and write point-polygon intersection list.
  PointPolygonIntersectionList pointPolygonIntersectionList;
  IntersectPointsWithPolygons(pointPolygonIntersectionList, converter.GetPoints(), converter.GetPolygons());

  // Finding the closest point contacts.
  FindClosestContactForPoints(pointPolygonIntersectionList, contacts, converter.GetFaults(),
                              converter, par.pointToContactDistanceBuffer);

  std::string fileNamePointPolygon = par.path_output + "/point-polygon-intersection.txt";
  WritePointPolygonIntersectionList(fileNamePointPolygon, pointPolygonIntersectionList);

  //-----------------------------------------------------------------------------------
  // Building unit contacts.
  UnitContacts unit_contacts;
  converter.BuildUnitContactsList(contacts, unit_contacts);
  converter.IdentifyIgneousUnitContacts(unit_contacts);
  converter.AddPointsToUnitContacts(unit_contacts, pointPolygonIntersectionList);
  converter.AddPointsToUnits(pointPolygonIntersectionList);

  //------------------------------------------------------------------------
  // Writing full graphs with all, strat, fault, and igneous contacts.
  //------------------------------------------------------------------------
  GraphWriter graph;

  for (size_t i = 0; i < par.depositListForGraphInfo.size(); i++) {
      std::string depositName = par.depositListForGraphInfo[i];

      // ALL contacts. -----------------------------------------------
      std::string file_graph_all = par.path_output + "/graph_all_" + depositName + ".gml";
      graph.WriteGraph(file_graph_all, unit_contacts, "UNIT_NAME", false, 1,
                       par.graph_edge_width_categories, par.graph_edge_direction_type,
                       depositName, graph_comments);

      // Stratigraphic unit contacts. --------------------------------
      UnitContacts unit_contacts_strat =
              Converter::FilterUnitStratFaultContacts(unit_contacts, StratigraphicContact, par.minFractionInMixedContact);

      std::string file_graph_strat = par.path_output + "/graph_strat_" + depositName + ".gml";
      graph.WriteGraph(file_graph_strat, unit_contacts_strat, "UNIT_NAME", false, 1,
                       par.graph_edge_width_categories, par.graph_edge_direction_type,
                       depositName, graph_comments);

      // Fault unit contacts. ----------------------------------------
      UnitContacts unit_contacts_fault =
              Converter::FilterUnitStratFaultContacts(unit_contacts, FaultContact, par.minFractionInMixedContact);

      std::string file_graph_fault = par.path_output + "/graph_fault_" + depositName + ".gml";
      graph.WriteGraph(file_graph_fault, unit_contacts_fault, "UNIT_NAME", false, 1,
                       par.graph_edge_width_categories, par.graph_edge_direction_type,
                       depositName, graph_comments);

       // Igneous unit contacts. -------------------------------------
      UnitContacts unit_contacts_ingeous = Converter::FilterIgneousUnitContacts(unit_contacts);

      std::string file_graph_ign = par.path_output + "/graph_igneous_" + depositName + ".gml";
      graph.WriteGraph(file_graph_ign, unit_contacts_ingeous, "UNIT_NAME", false, 1,
                       par.graph_edge_width_categories, par.graph_edge_direction_type,
                       depositName, graph_comments);
  }

  //----------------------------------------------------------------------------------------
  // Generating a partial graph (around one polygon).
  //----------------------------------------------------------------------------------------
  {
      Contacts contacts_filtered;
      converter.FilterContactsByPolygon(contacts, contacts_filtered, par.partial_graph_polygon_id, par.partial_graph_depth);

      unit_contacts.clear();
      converter.BuildUnitContactsList(contacts_filtered, unit_contacts);

      std::string file_graph = par.path_output + "/graph_partial.gml";
      graph.WriteGraph(file_graph, unit_contacts, "UNIT_NAME", false, 1,
                       par.graph_edge_width_categories, par.graph_edge_direction_type,
                       graph_comments);
  }
}

}

