/******************************************************************************
* converter.cpp
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#ifndef converter_cpp
#define converter_cpp

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <string>

#include <cassert>

#include "converter.h"
#include "svgbuilder.h"
#include "data_reader.h"
#include "intersector.h"

using namespace ClipperLib;

namespace ConverterLib {

void Converter::SetOutputFolder(const std::string &path)
{
  path_output = path;
}

//========================================================================================

void Converter::ReadData(const std::string &filename, const std::string& keyword,
                         const std::map<std::string, std::string>& constNames,
                         const std::vector<int>& idsToRead,
                         const std::vector<std::string>& idsToReadString)
{
  if (keyword == "POLYGON")
  {
    ReadDataObj(filename, keyword, constNames, polygons_read, idsToRead);
    CalculateBoundingBoxesObj(polygons_read);
    polygons = polygons_read;
  }
  else if (keyword == "LINESTRING")
  {
    ReadDataObj(filename, keyword, constNames, faults_read, idsToRead);
    CalculateBoundingBoxesObj(faults_read);
    faults = faults_read;
  }
  else if (keyword == "POINT")
  {
    ReadDataObj(filename, keyword, constNames, points_read, idsToRead, idsToReadString);
    points = points_read;
  }
  else
  {
    std::cout << "Unknown keyword in ReadData!" << std::endl;
  }
}
//========================================================================================

AABB Converter::GetAllPolygonsAABB() const
{
  Paths paths;
  for (Objects::const_iterator it = polygons_read.begin(); it != polygons_read.end(); ++it)
  {
    paths.push_back(AABB::AABB2Path(it->aabb));
  }
  AABB box = AABB::CalculateBoundingBox(paths);
  return box;
}

//========================================================================================
// Returns perimeter length of a multi-polygon.
static double get_polygon_length(const Object& polygon)
{
  double poly_length = 0.0;

  // Loop over all polygons within a multi-polygon.
  for (std::size_t i = 0; i < polygon.paths.size(); i++) {
    // Sanity check on the format of polygons: last vertex should be the same as first one.
    if (polygon.paths[i].back() != polygon.paths[i].front()) {
      std::cout << "Wrong polygons format: last vertex /= first vertex!" << std::endl;
      exit(0);
    }

    // Loop over all edges of a polygon.
    for (std::size_t j = 0; j < polygon.paths[i].size() - 1; j++) {
      IntPoint v1 = polygon.paths[i][j];
      IntPoint v2 = polygon.paths[i][j + 1];

      double length = ConverterUtils::GetDistance(v1, v2);
      poly_length += length;
    }
  }
  return poly_length;
}

void Converter::BuildUnitsAndGroupsLists()
{
  std::cout << "Building lists of units and groups..." << std::endl;

  units.clear();
  groups.clear();

  // Build units list.
  for (Objects::const_iterator poly_it = polygons.begin(); poly_it != polygons.end(); ++poly_it) {
    // Perimeter length of a multi-polygon.
    double poly_length = get_polygon_length(*poly_it);

    Units::iterator unit_it = units.find(poly_it->name);
    if (unit_it != units.end()) {
      // Adding perimeter length of another multi-polygon of the same unit type.
      unit_it->second.length += poly_length;
      unit_it->second.codes.push_back(poly_it->code);

    } else {
      // Adding a new unit.
      Unit unit;
      unit.id           = (int)(units.size());
      unit.name         = poly_it->name;
      unit.min_age      = poly_it->min_age;
      unit.max_age      = poly_it->max_age;

      std::string group = poly_it->group;
      if (group == "") {
          // Use unit name when group is empty.
          group = unit.name;
      }
      unit.group        = group;

      unit.rocktype1    = poly_it->rocktype1;
      unit.rocktype2    = poly_it->rocktype2;

      unit.length   = poly_length;
      unit.codes.push_back(poly_it->code);

      // Marking sills.
      std::string description = poly_it->description;
      std::size_t pos = description.find(SILL_DESCRIPTION);
      if (pos == std::string::npos) {
        unit.is_sill = false;
      } else {
        unit.is_sill = true;
      }

      units.insert(Units::value_type(unit.name, unit));
    }
  }

  std::cout << "Number of units: " << units.size() << std::endl;

  //---------------------------------------------------------------------
  // Calculate total unit perimeter length.
  double total_length = 0.0;
  for (Units::const_iterator unit_it = units.begin(); unit_it != units.end(); ++unit_it)
  {
    total_length += unit_it->second.length;
  }

  m_units_length = total_length;

  std::cout << std::setprecision(14) << "Total unit length = " << total_length << std::endl;

  //---------------------------------------------------------------------
  // Build groups list.
  for (Units::const_iterator unit_it = units.begin(); unit_it != units.end(); ++unit_it)
  {
    std::string group_name = unit_it->second.group;
    Groups::iterator group_it = groups.find(group_name);
    if (group_it == groups.end())
    {
      // Adding a new group.
      groups.insert(Groups::value_type(group_name, (int)(groups.size())));

      std::cout << groups.size() << " " << group_name << std::endl;
    }
  }

  //---------------------------------------------------------------------
  // Sort units by their length.

  std::vector<Unit> units_sorted;
  for (Units::const_iterator it = units.begin(); it != units.end(); ++it)
  {
    units_sorted.push_back(it->second);
  }

  // Perform sorting.
  std::sort(units_sorted.begin(), units_sorted.end());

  // Write sorted units into a file.
  std::ofstream file;
  file.open((this->path_output + "/units_sorted.txt").c_str());

  for (std::vector<Unit>::const_iterator unit_it = units_sorted.begin(); unit_it != units_sorted.end(); ++unit_it)
  {
    file << unit_it - units_sorted.begin() << " " << unit_it->length << " "
         << unit_it->min_age << " " << unit_it->max_age << " "
         << unit_it->is_sill << " " << unit_it->name << std::endl;
  }
  file.close();
}
//========================================================================================

void Converter::ClipData(const AABB &window)
{
  // Clipping polygons.
  polygons = ClipDataObj(polygons_read, window, true);
  std::cout << "Polygons left (after clipping): " << polygons.size() << std::endl;

  // Setting last vertex = first vertex, for polygons. As this gets broken after clipping.
  for (Objects::iterator poly_it = polygons.begin(); poly_it != polygons.end(); poly_it++)
  {
    for (Paths::iterator path_it = poly_it->paths.begin(); path_it != poly_it->paths.end(); path_it++)
    {
      if (path_it->front() != path_it->back())
      {
        path_it->push_back(path_it->front());
      }
    }
  }

  // Clipping faults.
  faults = ClipDataObj(faults_read, window, false);
  std::cout << "Faults left (after clipping): " << faults.size() << std::endl;

  // Clipping points.
  points.clear();
  for (Objects::const_iterator it = points_read.begin();
       it != points_read.end(); ++it) {
      const IntPoint& pt = it->paths[0][0];
      if (AABB::PointInsideBox(pt, window)) {
          points.push_back(*it);
      }
  }
  std::cout << "Points left (after clipping): " << points.size() << std::endl;

  //----------------------------------------------------------------------------------
  // Calculate window perimeter.
  Path window_path = AABB::AABB2Path(window);

  double perimeter = 0.0;
  for (Path::const_iterator it = window_path.begin(); it != window_path.end() - 1; it++) {
      perimeter += ConverterUtils::GetDistance(*it, *(it + 1));
  }
  m_window_length = perimeter;

  std::cout << std::setprecision(14) << "Window perimeter = " << perimeter << std::endl;
}
//========================================================================================

Objects Converter::ClipDataObj(const Objects &objects, const AABB &window, bool closed_path)
{
  // Skip clipping if zero window is set.
  if (window.BtmRight == window.TopLeft) return objects;

  // Building the window polygon.
  Paths window_paths;
  window_paths.push_back(AABB::AABB2Path(window));

  PolyTree polytree;
  Paths solution;

  Objects objects_clipped;

  for (Objects::const_iterator it = objects.begin(); it != objects.end(); ++it)
  {
    Clipper clpr;
    clpr.Clear();

    // Adding data to the clipper.
    clpr.AddPaths(it->paths, ptSubject, closed_path);
    clpr.AddPaths(window_paths, ptClip, true);

    // Do the intersection.
    if (!clpr.Execute(ctIntersection, polytree, pftEvenOdd, pftEvenOdd))
    {
      std::cout << "Error in Execute() in ClipDataObj()!" << std::endl;
      exit(0);
    }

    // Retrieve the solution.
    if (closed_path)
    {
      ClosedPathsFromPolyTree(polytree, solution);
    } else
    {
      OpenPathsFromPolyTree(polytree, solution);
    }

    if (solution.size() != 0)
    {// An object intersects the clipping window.
      // Add the clipped object.
      Object object = *it;
      object.paths = solution;
      objects_clipped.push_back(object);
    }
  }
  return objects_clipped;
}
//========================================================================================

double Converter::FindContacts(Contacts &contacts) const
{
  std::cout << "Finding contacts..." << std::endl;

  int numCommonEdges = 0;

  // N^2 loop for every pair of multipolygons.
  for (Objects::const_iterator it_i = polygons.begin(); it_i != polygons.end(); ++it_i)
  {
    for (Objects::const_iterator it_j = it_i + 1; it_j != polygons.end(); ++it_j)
    {
      // Bounding box test.
      if (AABB::BoundingBoxesOverlap(it_i->aabb, it_j->aabb))
      {
        numCommonEdges += IntersectPolygons(*it_i, *it_j, contacts);
      }
    }
  }
  std::cout << "numCommonEdges = " << numCommonEdges << std::endl;

  // Assign unique id to contacts.
  AssignContactIds(contacts);

  // Calculate contact lengths.
  double total_length = 0.0;
  for (Contacts::iterator it_contact = contacts.begin(); it_contact != contacts.end(); ++it_contact)
  {
    double length = 0.0;

    for (Path::const_iterator it_vertex = it_contact->path.begin();
         it_vertex != it_contact->path.end() - 1; ++it_vertex)
    {
      // Add up a contact segment length.
      length += ConverterUtils::GetDistance(*it_vertex, *(it_vertex + 1));
    }

    it_contact->length = length;
    total_length += length;
  }

  std::cout << "Total contact length = " << total_length << std::endl;
  return total_length;
}
//========================================================================================

const Units& Converter::GetUnits() const
{
  return units;
}

const Groups& Converter::GetGroups() const
{
  return groups;
}

const Objects& Converter::GetPolygons() const
{
  return polygons;
}

const Objects& Converter::GetFaults() const
{
  return faults;
}

const Objects& Converter::GetPoints() const
{
  return points;
}

void Converter::PrintLengths(double contacts_length) const
{
  // Test: two lengths should be the same, if there is a proper clipping window.
  double length1 = m_window_length + 2.0 * contacts_length;
  double length2 = m_units_length;
  std::cout << "PERIMETER LENGTHS: " << length1 << " " << length2 << " " << length1 - length2 << std::endl;
}
//========================================================================================

void Converter::BuildUnitContactsList(const Contacts &poly_contacts, UnitContacts &unit_contacts) const
{
  std::cout << "Finding unit-contacts..." << std::endl;

  // Loop over all polygon contacts.
  for (Contacts::const_iterator poly_contact_it = poly_contacts.begin();
       poly_contact_it != poly_contacts.end(); ++poly_contact_it)
  {
    if (poly_contact_it->obj1 == 0 || poly_contact_it->obj2 == 0)
    {
      std::cout << "Error: null object at contact!" << std::endl;
      exit(0);
    }

    std::string name1 = poly_contact_it->obj1->name;
    std::string name2 = poly_contact_it->obj2->name;

    // Skipping contacts between units with the same name (happens e.g. due to faults).
    if (name1 == name2) continue;

    int exist = false;

    for (UnitContacts::iterator unit_contact_it = unit_contacts.begin();
         unit_contact_it != unit_contacts.end(); ++unit_contact_it)
    {
      if ((unit_contact_it->unit1->name == name1 && unit_contact_it->unit2->name == name2)
          || (unit_contact_it->unit1->name == name2 && unit_contact_it->unit2->name == name1))
      {
      // This pair of unit names is already added.
        exist = true;

        // Add up the total contact length.
        unit_contact_it->total_length += poly_contact_it->length;
        //! TODO: move this block with the same one below to a function.
        // Add up partial lengths.
        if (poly_contact_it->type == StratigraphicContact)
        {
          unit_contact_it->stratigraphic_length += poly_contact_it->length;
        }
        else if (poly_contact_it->type == FaultContact)
        {
          unit_contact_it->faults_length += poly_contact_it->length;
        }
        else if (poly_contact_it->type == MixedContact)
        {
          unit_contact_it->mixed_length += poly_contact_it->length;
        }
        else
        {
          std::cout << "Wrong contact type! type = " << poly_contact_it->type << std::endl;
          exit(0);
        }

        // Store the corresponding polygon contact.
        unit_contact_it->contacts.push_back(*poly_contact_it);

        break;
      }
    }

    if (!exist)
    {
      // Adding a new unit-contact.
      UnitContact unit_contact;
      unit_contact.total_length = poly_contact_it->length;

      unit_contact.faults_length = 0.0;
      unit_contact.stratigraphic_length = 0.0;
      unit_contact.mixed_length = 0.0;

      // Set partial lengths.
      if (poly_contact_it->type == StratigraphicContact)
      {
        unit_contact.stratigraphic_length = poly_contact_it->length;
      }
      else if (poly_contact_it->type == FaultContact)
      {
        unit_contact.faults_length = poly_contact_it->length;
      }
      else if (poly_contact_it->type == MixedContact)
      {
        unit_contact.mixed_length = poly_contact_it->length;
      }
      else
      {
        std::cout << "Wrong contact type! type = " << poly_contact_it->type << std::endl;
        exit(0);
      }

      unit_contact.contacts.push_back(*poly_contact_it);

      // Setting pointers to the contacting units.
      Units::const_iterator unit_it;
      unit_it = units.find(name1);
      if (unit_it != units.end())
      {
        unit_contact.unit1 = &(unit_it->second);
      } else
      {
        std::cout << "Error: could not find a unit!" << std::endl;
        exit(0);
      }

      unit_it = units.find(name2);
      if (unit_it != units.end())
      {
        unit_contact.unit2 = &(unit_it->second);
      } else
      {
        std::cout << "Error: could not find a unit!" << std::endl;
        exit(0);
      }

      unit_contacts.push_back(unit_contact);
    }
  }

  //------------------------------------------------------------------------------
  // Determine unit contact type.

  int num_strat_contacts = 0;
  int num_falut_contacts = 0;
  int num_mixed_contacts = 0;

  // Loop over all units.
  for (UnitContacts::iterator unit_contact_it = unit_contacts.begin();
       unit_contact_it != unit_contacts.end();
       unit_contact_it++)
  {
    bool hasStratigraphicContacts = false;
    bool hasFaultContacts = false;
    bool hasMixedContacts = false;

    // Loop over all contacts within a unit.
    for (Contacts::const_iterator poly_contact = unit_contact_it->contacts.begin();
         poly_contact != unit_contact_it->contacts.end();
         poly_contact++)
    {
      if (poly_contact->type == StratigraphicContact) hasStratigraphicContacts = true;
      if (poly_contact->type == FaultContact) hasFaultContacts = true;
      if (poly_contact->type == MixedContact) hasMixedContacts = true;
    }

    // Should not have mixed contacts after polygon contacts splitting.
    assert(!hasMixedContacts);

    if (hasMixedContacts) unit_contact_it->type = MixedContact;
    if (hasStratigraphicContacts && hasFaultContacts) unit_contact_it->type = MixedContact;
    if (hasStratigraphicContacts && !hasFaultContacts && !hasMixedContacts) unit_contact_it->type = StratigraphicContact;
    if (!hasStratigraphicContacts && hasFaultContacts && !hasMixedContacts) unit_contact_it->type = FaultContact;

    if (unit_contact_it->type == StratigraphicContact) num_strat_contacts++;
    if (unit_contact_it->type == FaultContact) num_falut_contacts++;
    if (unit_contact_it->type == MixedContact) num_mixed_contacts++;
  }

  std::cout << "Number of unit-contacts = " << unit_contacts.size() << std::endl;
  std::cout << "Number of stratigraphic unit-contacts = " << num_strat_contacts << std::endl;
  std::cout << "Number of fault unit-contacts = " << num_falut_contacts << std::endl;
  std::cout << "Number of mixed unit-contacts = " << num_mixed_contacts << std::endl;

  //------------------------------------------------------------------------------
  // Write distribution of (relative and absolute) contact lengths.

  std::vector<double> relative_lengths;
  std::vector<double> absolute_lengths;

  // Loop over all unit contacts.
  for (UnitContacts::iterator unit_contact_it = unit_contacts.begin();
       unit_contact_it != unit_contacts.end();
       unit_contact_it++)
  {
    double length1 = unit_contact_it->unit1->length;
    double length2 = unit_contact_it->unit2->length;

    unit_contact_it->relative_length = unit_contact_it->total_length / (length1 + length2) / 2.0;

    relative_lengths.push_back(unit_contact_it->relative_length);
    absolute_lengths.push_back(unit_contact_it->total_length);
  }

  // Sort lengths.
  std::sort(relative_lengths.begin(), relative_lengths.end());
  std::sort(absolute_lengths.begin(), absolute_lengths.end());

  // Write sorted relative lengths to a file.
  std::ofstream file;
  file.open((this->path_output + "/units_contacts_relative_lengths.txt").c_str());
  for (std::vector<double>::const_iterator it = relative_lengths.begin(); it != relative_lengths.end(); ++it)
  {
    file << *it << std::endl;
  }
  file.close();

  // Write sorted absolute lengths to a file.
  file.open((this->path_output + "/units_contacts_absolute_lengths.txt").c_str());
  for (std::vector<double>::const_iterator it = absolute_lengths.begin(); it != absolute_lengths.end(); ++it)
  {
    file << *it << std::endl;
  }
  file.close();
}
//========================================================================================
int Converter::TestForIgneousContact(const std::string& rocktype1,
                                     const std::string& rocktype2,
                                     double age1, double age2) const
{
    int res = 0;

    bool contains_i_1 = (rocktype1.find(IGNEOUS_STRING) != std::string::npos);
    bool contains_i_2 = (rocktype2.find(IGNEOUS_STRING) != std::string::npos);

    bool contains_v_1 = (rocktype1.find(VOLCANIC_STRING) != std::string::npos);
    bool contains_v_2 = (rocktype2.find(VOLCANIC_STRING) != std::string::npos);

    bool igneous1 = contains_i_1 && !contains_v_1;
    bool igneous2 = contains_i_2 && !contains_v_2;

    if (igneous1 || igneous2) {
        // Intrusive igneous.
        res = 1;
    }

    if ((igneous1 && igneous2)
        || (igneous1 && (age1 < age2))
        || (igneous2 && (age2 < age1))) {
        // Igneous.
        res = 2;
    }

    return res;
}

void Converter::IdentifyIgneousUnitContacts(UnitContacts &unitContacts) const
{
    std::cout << "Identifying igneous contacts..." << std::endl;
    size_t numIntrusiveIgneous = 0;
    size_t numIgneous = 0;

    for (UnitContacts::iterator it = unitContacts.begin();
         it != unitContacts.end(); ++it) {
        // Build the rocktype form rocktype1 and rocktype2 fields (extracted from shapefile).
        std::string rocktype1 = it->unit1->rocktype1 +  " " + it->unit1->rocktype2;
        std::string rocktype2 = it->unit2->rocktype1 +  " " + it->unit2->rocktype2;

        int res = TestForIgneousContact(rocktype1, rocktype2,
                                        it->unit1->max_age, it->unit2->max_age);

        if (res > 0) {
            it->isIntrusiveIgneous = true;
            numIntrusiveIgneous++;
        }

        if (res == 2) {
            // A potentially (partially) igneous contact (unless its full length = fault length).
            it->isIgneous = true;
            numIgneous++;
        }
    }

    std::cout << "Number of intrusive igneous unit contacts: " << numIntrusiveIgneous << std::endl;
    std::cout << "Number of igneous unit contacts: " << numIgneous << std::endl;
}
//========================================================================================

void Converter::AddPointsToUnitContacts(UnitContacts &unitContacts,
                                        const PointPolygonIntersectionList& pointPolygonIntersectionList) const
{
    // Loop over unit contacts.
    for (UnitContacts::iterator it = unitContacts.begin();
         it != unitContacts.end(); ++it) {

        // Loop over points (deposits).
        for (std::size_t i = 0; i < pointPolygonIntersectionList.size(); i++) {
            const PointInPolygon& pip = pointPolygonIntersectionList[i];

            if (pip.onContact) {
                if ((it->unit1->name == pip.containingPolygon->name
                     && it->unit2->name == pip.neighbourPolygon->name)
                    || (it->unit2->name == pip.containingPolygon->name
                        && it->unit1->name == pip.neighbourPolygon->name)) {
                    // Adding deposit on this contact.
                    it->deposits.push_back(&pip);
                }
            }
        }
    }
}
//========================================================================================

void Converter::AddPointsToUnits(const PointPolygonIntersectionList& pointPolygonIntersectionList)
{
    // Loop over points (deposits).
    for (std::size_t i = 0; i < pointPolygonIntersectionList.size(); i++) {
        const PointInPolygon& pip = pointPolygonIntersectionList[i];
        if (!pip.onContact) {
            const std::string unitName = pip.containingPolygon->name;

            Units::iterator it = units.find(unitName);
            if (it != units.end()) {
                it->second.deposits.push_back(&pip);
            } else {
                assert(false);
            }
        }
    }
}
//========================================================================================

void Converter::IdentifyPolygonContactTypes(Contacts &contacts, double angleEpsilon, double distanceEpsilon) const
{
  std::cout << "Identifying contact types..." << std::endl;

  int num_fc = 0;
  int num_mc = 0;

  std::size_t total_num_intersections = 0;

  Contacts::iterator it = contacts.begin();

  // Loop over all contacts.
  while (it != contacts.end())
  {
    Contact contact = *it;

    Paths contact_paths;
    contact_paths.push_back(contact.path);
    AABB contact_aabb = AABB::CalculateBoundingBox(contact_paths);

    // Mark initially all contact vertexes (segments) as non-intersecting.
    for (Path::iterator vertex_it = contact.path.begin(); vertex_it != contact.path.end(); ++vertex_it)
    {
      vertex_it->type = 0;
    }

    // Loop over all faults.
    for (Objects::const_iterator fault_it = faults.begin(); fault_it != faults.end(); ++fault_it)
    {
      // Bounding boxes overlap test (for optimization).
      if (AABB::BoundingBoxesOverlap(contact_aabb, fault_it->aabb))
      {
        // Mark coinciding contact-with-fault segments.
        size_t num_intersections0 = IntersectContactWithFault(contact, *fault_it, angleEpsilon, distanceEpsilon);

        if (num_intersections0 == contact.path.size() - 1)
        {
          //std::cout << "Full contact of polygons: " << contact.obj1->id << ", " << contact.obj2->id <<
          //             " with a fault: " << fault_it->id << std::endl;
          break;
        }
      }
    }

    std::size_t num_intersections = 0;
    for (Path::const_iterator vertex_it = contact.path.begin(); vertex_it != contact.path.end() - 1; ++vertex_it)
    {
      if (vertex_it->type == 1) num_intersections++;
    }

    total_num_intersections += num_intersections;

    if (num_intersections == contact.path.size() - 1)
    {
    // All contact segments coincide with fault(s). Mark as fault-contact.
      num_fc++;
      contact.type = FaultContact;
    }
    else if (num_intersections > 0)
    {
    // Some contact segments coincide with fault(s). Mark as mixed-contact.
      num_mc++;
      contact.type = MixedContact;
    }
    else
    {
      contact.type = StratigraphicContact;
    }

    // Replace the contact with a modified contact.
    *it = contact;

    it++;
  }

  std::cout << "Total number of intersections: " << total_num_intersections << std::endl;
  std::cout << "Number of all contacts: " << contacts.size() << std::endl;
  std::cout << "Number of fault-contacts: " << num_fc << std::endl;
  std::cout << "Number of mixed-contacts: " << num_mc << std::endl;
}
//========================================================================================

void Converter::SplitMixedPolygonContacts(Contacts &contacts)
{
  std::cout << "Splitting mixed contacts..." << std::endl;

  Contacts split_contacts;

  Contacts::iterator contact_it = contacts.begin();

  // Loop over all contacts.
  while (contact_it != contacts.end())
  {
    if (contact_it->type != MixedContact)
    {
      ++contact_it;
      continue;
    }

    Contact new_contact;
    new_contact.obj1 = contact_it->obj1;
    new_contact.obj2 = contact_it->obj2;
    new_contact.length =  0.0;

    double total_length = 0.0;

    // Loop over all vertexes of a contact.
    for (Path::const_iterator vertex_it = contact_it->path.begin();
         vertex_it != contact_it->path.end() - 1; ++vertex_it)
    {
      new_contact.path.push_back(*vertex_it);
      new_contact.length += ConverterUtils::GetDistance(*vertex_it, *(vertex_it + 1));

      if (vertex_it->type != (vertex_it + 1)->type
          || vertex_it == contact_it->path.end() - 2)
      {
        new_contact.path.push_back(*(vertex_it + 1));

        if (vertex_it->type == 0)
        {
          new_contact.type = StratigraphicContact;
        } else
        {
          new_contact.type = FaultContact;
        }

        split_contacts.push_back(new_contact);
        total_length += new_contact.length;

        new_contact.length = 0.0;
        new_contact.path.clear();
      }
    }
    // Remove a mixed contact.
    contact_it = contacts.erase(contact_it);
  }

  // Adding the split contacts.
  contacts.insert(contacts.end(), split_contacts.begin(), split_contacts.end());

  // Assign unique id to contacts.
  AssignContactIds(contacts);

  std::cout << "Number of contacts after splitting: " << contacts.size() << std::endl;
}
//========================================================================================

void Converter::AssignContactIds(Contacts &contacts)
{
  for (Contacts::iterator contact_it = contacts.begin();
       contact_it != contacts.end(); ++contact_it)
  {
    contact_it->id = (int)(contact_it - contacts.begin());
  }
}
//========================================================================================

// Use poly_last_id to avoid moving backward in the graph, during recursion, for optimization.
void Converter::FilterContactsByPolygon(const Contacts &contacts, Contacts &contacts_filtered,
                                        int poly_id, int depth, int poly_last_id, bool last_call)
{
  // Building a list of polygon ids, that are in contact with a given polygon.
  // Also adding the corresponding contacts to a container of filtered contacts.
  std::vector<int> poly_next_ids;
  for (Contacts::const_iterator contact_it = contacts.begin();
       contact_it != contacts.end(); ++contact_it)
  {
    // Sanity check.
    if (contact_it->id == - 1)
    {
      std::cout << "Wrong contact id in FilterContacts!" << std::endl;
      exit(0);
    }

    if (contact_it->obj1 != NULL && contact_it->obj2 != NULL)
    {
      if (contact_it->obj1->id == poly_id && contact_it->obj2->id != poly_last_id)
      {
        contacts_filtered.push_back(*contact_it);
        poly_next_ids.push_back(contact_it->obj2->id);
      }
      else if (contact_it->obj2->id == poly_id && contact_it->obj1->id != poly_last_id)
      {
        contacts_filtered.push_back(*contact_it);
        poly_next_ids.push_back(contact_it->obj1->id);
      }
    }
    else
    {
      std::cout << "Error: null object at contact!" << std::endl;
    }
  }

  if (depth > 1)
  {
    // Search for the higher degree contacts, i.e., for a graph A-B-C, nodes A and C are in 2nd-degree contact.
    for (std::vector<int>::const_iterator it = poly_next_ids.begin(); it != poly_next_ids.end(); ++it)
    {
      bool last_call_next = (last_call && (it == poly_next_ids.end() - 1));
      FilterContactsByPolygon(contacts, contacts_filtered, *it, depth - 1, poly_id, last_call_next);
    }
  }
  else
  {
    if (last_call)
    { // This is the last call of a recursive function, so the container of filtered contacts is complete.
      std::cout << "Number of filtered contacts: " << contacts_filtered.size() << std::endl;

      // Remove duplicate contacts (via assigning to a std::set).
      std::set<Contact> s(contacts_filtered.begin(), contacts_filtered.end());
      contacts_filtered.assign(s.begin(), s.end());

      std::cout << "Number of unique filtered contacts: " << contacts_filtered.size() << std::endl;
    }
  }
}
//========================================================================================

UnitContacts Converter::FilterUnitStratFaultContacts(const UnitContacts &contacts, ContactType contactType,
                                                     double minFractionInMixedContact)
{
    assert(contactType == StratigraphicContact || contactType == FaultContact);

    UnitContacts filteredContacts;

    // Loop over all contacts.
    for (UnitContacts::const_iterator it = contacts.cbegin();
         it != contacts.cend(); ++it) {
        bool keepContact = false;

        if (it->type == contactType) {
            keepContact = true;
        }

        if (it->type == MixedContact) {
            if (contactType == StratigraphicContact) {
                if (it->stratigraphic_length / it->total_length >= minFractionInMixedContact) {
                    keepContact = true;
                }
            }
            if (contactType == FaultContact) {
                if (it->faults_length / it->total_length >= minFractionInMixedContact) {
                    keepContact = true;
                }
            }
        }

        if (keepContact) {
            filteredContacts.push_back(*it);
        }
    }
    return filteredContacts;
}
//========================================================================================

UnitContacts Converter::FilterIgneousUnitContacts(const UnitContacts &contacts)
{
    UnitContacts filteredContacts;

    // Loop over all contacts.
    for (UnitContacts::const_iterator it = contacts.cbegin();
         it != contacts.cend(); ++it) {

        if (it->isIntrusiveIgneous && it->type != FaultContact) {
            filteredContacts.push_back(*it);
        }
    }
    return filteredContacts;
}
//========================================================================================

void Converter::WriteImage(const std::string &filename, const std::vector<AABB> &windows, int add_objects,
                           const Contacts *contacts) const
{
  Paths subject, clip, solution;
  PolyTree polytree;

  clip.clear();
  subject.clear();
  polytree.Clear();

  // Adding all polygons to clip object for visualization.
  for (Objects::const_iterator it = polygons.begin(); it != polygons.end(); ++it)
  {
    clip.insert(clip.end(), it->paths.begin(), it->paths.end());
  }

  // Add window rectangles.
  for (std::vector<AABB>::const_iterator it = windows.begin(); it != windows.end(); ++it)
  {
    AABB window = *it;
    if (window.TopLeft != window.BtmRight)
      subject.push_back(AABB::AABB2Path(window));
  }

  if (add_objects != 0)
  { // Adding points/faults to the image.
    const Objects *objects = 0;
    if (add_objects == 1)
    {
      objects = &points;
    }
    else
    {
      objects = &faults;
    }
    for (Objects::const_iterator it = objects->begin(); it != objects->end(); it++)
    {
      subject.insert(subject.end(), it->paths.begin(), it->paths.end());
    }
  }

  if (contacts != 0)
  { // Inserting all contacts.
    for (Contacts::const_iterator it_contact = contacts->begin(); it_contact != contacts->end(); ++it_contact)
    {
      subject.push_back(it_contact->path);
    }
  }

  // Write SVG graphics file.
  SvgBuilder::ClipperData *clt = new SvgBuilder::ClipperData();

  clt->title = this->path_output + "/" + filename;
  clt->subj2 = subject;
  clt->clip = clip;
  clt->solution = solution;
  clt->polytree = polytree;
  clt->pft = pftEvenOdd;

  //! TODO: pass as input parameters.
  bool showPolyVertices = true;
  int polyVerticesType = 0;

  SvgBuilder::DoSVG(clt, showPolyVertices, polyVerticesType);
}
//========================================================================================

void Converter::CalculateBoundingBoxesObj(Objects &objects)
{
  // Loop over all objects
  for (std::size_t i = 0; i < objects.size(); i++)
  {
    if (objects[i].paths.size() == 0) continue;
    if (objects[i].paths[0].size() == 0) continue;

    objects[i].aabb = AABB::CalculateBoundingBox(objects[i].paths);
  }
}
//========================================================================================

void Converter::MakePolygonFromInts(const cInt *ints, int size, Path &p, double scale)
{
  p.clear();
  p.reserve(size / 2);
  for (int i = 0; i < size; i += 2)
    p.push_back(IntPoint((cInt)(ints[i] * scale), (cInt)(ints[i + 1] * scale)));
}

}

#endif

