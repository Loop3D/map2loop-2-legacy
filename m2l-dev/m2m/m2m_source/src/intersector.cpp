/******************************************************************************
* A class for intersecting faults & polygons.
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#include <iostream>
#include <exception>
#include <cassert>

#include "clipper.hpp"

#include "intersector.h"

namespace ConverterLib {

int PointInMultipolygon(const IntPoint &pt, const Object &multipolygon)
{
    assert(multipolygon.paths.size() > 0);
    std::vector<int> inPolygonResults(multipolygon.paths.size());

    // Assume that the first polygon is the "main" one, and the followings are its holes.
    // So a point must be inside the "main" polygon, but not inside any of its holes.
    for (std::size_t i = 0; i < multipolygon.paths.size(); i++) {
        const Path &polygon = multipolygon.paths[i];
        inPolygonResults[i] = ClipperLib::PointInPolygon(pt, polygon);
    }

    //----------------------------------------------------------------------
    // Check if a point lies on boundary.
    for (std::size_t i = 0; i < multipolygon.paths.size(); i++) {
        if (inPolygonResults[i] == -1) {
            return -1;
        }
    }

    //---------------------------------------------------------------------
    // Check if a point is inside a multipolygon, and not inside its hole.
    bool inside = false;
    for (std::size_t i = 0; i < multipolygon.paths.size(); i++) {
        if (inPolygonResults[i] == 1) {
        // Point is inside.
            inside = true;

            const Path &polygon = multipolygon.paths[i];

            // This way it is tested for a hole in DoSVG().
            bool hole = !ClipperLib::Orientation(polygon);

            if (hole) {
                // Point is inside a hole.
                return 0;
            }
        }
    }
    if (inside) {
        return 1;
    } else {
        return 0;
    }
}

void IntersectPointsWithPolygons(PointPolygonIntersectionList &pointPolygonIntersectionList,
                                 const Objects &points, const Objects &polygons)
{
    std::cout << "Intersecting points with polygons..." << std::endl;

    for (std::size_t i = 0; i < points.size(); i++) {
        assert(points[i].paths.size() == 1);
        assert(points[i].paths[0].size() == 1);

        bool foundIntersectingPolygon = false;

        const IntPoint &pt = (points[i].paths[0])[0];

        size_t numInside = 0;
        for (std::size_t j = 0; j < polygons.size(); j++) {
            int res = PointInMultipolygon(pt, polygons[j]);

            if (res > 0) {
                foundIntersectingPolygon = true;
                numInside++;

                // Sanity check.
                if (numInside > 1) {
                    throw std::runtime_error("A point " + points[i].site_code + " is inside more than one multipolygon! This should not happen.");
                }
            }
            else if (res == -1) {
                foundIntersectingPolygon = true;
                std::cout << "Point is on boundary: " << points[i].site_code << ", " << polygons[j].id << std::endl;
            }

            if (res != 0) {
                PointInPolygon pip;
                pip.point = &points[i];
                pip.containingPolygon = &polygons[j];
                pip.onContact = (res == -1);

                pointPolygonIntersectionList.push_back(pip);
            }
        }

        if (!foundIntersectingPolygon) {
            std::cout << "WARNING: Did not find an intersecting polygon for a point with SITE CODE = " << points[i].site_code << std::endl;
        }
    }
}

void FindClosestContactForPoints(PointPolygonIntersectionList &pointPolygonIntersectionList,
                                 const Contacts& contacts, const Objects& faults,
                                 const Converter& converter, double pointToContactDistanceBuffer)
{
    std::cout << "Finding closest contacts (for points)..." << std::endl;

    //==================================================================================
    // 1. First check points against the contacts.
    //==================================================================================
    for (std::size_t i = 0; i < pointPolygonIntersectionList.size(); i++) {
        PointInPolygon& pip = pointPolygonIntersectionList[i];

        const IntPoint& pt = pip.point->paths[0][0];

        double minDist2;
        int closestContactIndex = -1;
        int faultID1 = -1;
        //int faultID2 = -1;

        // Finding a closest contact to a point.
        for (std::size_t j = 0; j < contacts.size(); j++) {
            const Contact& contact = contacts[j];

            if (contact.obj1 == pip.containingPolygon || contact.obj2 == pip.containingPolygon) {
            // This is a contact of a polygon where the point is located.

                // Need to find the (minimum) distance from a point to the contact.
                // Loop over all segments of a contact.
                for (Path::const_iterator it = contact.path.begin();
                     it != contact.path.end() - 1; ++it) {
                    double dist2 = ConverterUtils::DistanceFromSegmentSqrd(pt, *it, *(it + 1));

                    if (closestContactIndex < 0 || dist2 < minDist2) {
                        minDist2 = dist2;
                        closestContactIndex = j;

                        faultID1 = it->faultID;
                        //faultID2 = (it + 1)->faultID;
                    }
                }
            }
        }
        assert(closestContactIndex >= 0);

        const Contact& closestContact = contacts[closestContactIndex];

        if (closestContact.obj1 == pip.containingPolygon) {
            pip.neighbourPolygon = closestContact.obj2;
        }
        else if (closestContact.obj2 == pip.containingPolygon) {
            pip.neighbourPolygon = closestContact.obj1;
        }
        else {
            std::cout << "Wrong closest contact!" << std::endl;
            assert(false);
        }

        // Finding contact type -----------------------------
        std::string rocktype1 = closestContact.obj1->rocktype1 +  " " + closestContact.obj1->rocktype2;
        std::string rocktype2 = closestContact.obj2->rocktype1 +  " " + closestContact.obj2->rocktype2;
        double age1 = closestContact.obj1->max_age;
        double age2 = closestContact.obj2->max_age;

        int res = converter.TestForIgneousContact(rocktype1, rocktype2, age1, age2);

        if (res > 0) {
        // Igneouse/intrusive contacts get priority over fault contacts.
            if (res == 2) {
                pip.contactType = IgneousContact;
            } else {
                pip.contactType = IntrusiveIgneousContact;
            }
        } else {
            pip.contactType = closestContact.type;
        }
        //---------------------------------------------------

        if (pip.contactType == FaultContact) {
            assert(faultID1 >= 0);
            pip.faultObjectID = faultID1;
        }

        pip.distanceToContact2 = minDist2;
    }

    //==================================================================================
    // 2. Treat faults separately.
    // Note: faults are also considered as contacts, even if it is a fault inside a unit (i.e, between A and A).
    // Our conventional contacts don't contain A-A fault contacts, so we test against faults separately.
    //==================================================================================
    for (std::size_t i = 0; i < pointPolygonIntersectionList.size(); i++) {
        PointInPolygon& pip = pointPolygonIntersectionList[i];

        const IntPoint& pt = pip.point->paths[0][0];

        for (std::size_t j = 0; j < faults.size(); j++) {
            const Object &fault = faults[j];
            assert(fault.paths.size() == 1);

            for (Path::const_iterator it = fault.paths[0].begin();
                 it != fault.paths[0].end() - 1; ++it) {

                double dist2 = ConverterUtils::DistanceFromSegmentSqrd(pt, *it, *(it + 1));

                if (dist2 <= pip.distanceToContact2) {
                // Found a fault that is closer to the point than any of the contact.
                // So this should be a fault splitting the same unit, i.e., A-A contact.
                    pip.contactType = FaultContact;
                    pip.faultObjectID = fault.id;
                    pip.neighbourPolygon = pip.containingPolygon; // fault contact between A and A.
                    pip.distanceToContact2 = dist2;
                }
            }
        }
        // Setting the "onContact" flag.
        if (pip.distanceToContact2 <= pointToContactDistanceBuffer * pointToContactDistanceBuffer) {
            pip.onContact = true;
        } else {
            assert(pip.onContact == false);
        }
    }
}
//=====================================================================================================

void IntersectFaultsWithPolygons(const Objects &faults, const Objects &polygons,
                                 UnitFaultIntersectionList &unitFaultIntersectionList)
{
  std::cout << "Intersecting faults with polygons..." << std::endl;

  // Loop over all faults.
  for (std::size_t j = 0; j < faults.size(); j++) {
    // Loop over all multipolygons.
    for (std::size_t i = 0; i < polygons.size(); i++) {
      double distanceEpsilon = 5.;
      if (FaultAndPolygonIntersecting(faults[j], polygons[i], distanceEpsilon)) {

        std::string unitName = polygons[i].name;
        int faultId = faults[j].id;

        unitFaultIntersectionList[unitName].push_back(faultId);
      }
    }
  }
  // Remove duplicate fault IDs.
  for (UnitFaultIntersectionList::iterator it = unitFaultIntersectionList.begin();
       it != unitFaultIntersectionList.end(); ++it) {
    it->second.sort();
    it->second.unique();
  }
}
//=======================================================================================

bool FaultAndPolygonIntersecting(const Object &fault, const Object &polygon, double distanceEpsilon)
{
  // Bounding box check.
  if (!AABB::BoundingBoxesOverlap(fault.aabb, polygon.aabb)) {
    return 0;
  }

  //------------------------------------------------------------------------------------
  // 1. Clip a fault with a multipolygon and check if solution is non empty.
  // (Note: fault aligning with polygon boundary is not included in the solution).

  // Adding data to the clipper.
  ClipperLib::Clipper clpr;
  clpr.Clear();
  clpr.AddPaths(fault.paths, ClipperLib::ptSubject, false);
  clpr.AddPaths(polygon.paths, ClipperLib::ptClip, true);

  ClipperLib::PolyTree polytree;

  // Do the intersection.
  if (!clpr.Execute(ClipperLib::ctIntersection, polytree, ClipperLib::pftEvenOdd, ClipperLib::pftEvenOdd)) {
    throw std::runtime_error("Error in FaultAndPolygonIntersecting: Execute!");
  }

  Paths solution;
  OpenPathsFromPolyTree(polytree, solution);

  if (solution.size() > 0) {
    //-----------------------------------------------------------------------
    // Patch to not count fault-unit intersections where the fault only slightly intersects the unit.
    if (solution[0].size() == 2) {
      double dist2 = ConverterUtils::GetDistanceSquared(solution[0][0], solution[0][1]);
      if (dist2 <= distanceEpsilon * distanceEpsilon) {
          return false;
      }
    }
    //-----------------------------------------------------------------------
    return true;

  } else {
    //------------------------------------------------------------------------------------
    // 2. Check if fault points lie on the multipolygon boundary.
    size_t points_on_boundary = 0;

    assert(fault.paths.size() == 1);

    // Loop over all points of a fault.
    for (Path::const_iterator it = fault.paths[0].begin();
         it != fault.paths[0].end(); ++it) {

      // Loop over all polygons inside a multipolygon.
      for (Paths::const_iterator polyIt = polygon.paths.begin();
           polyIt != polygon.paths.end(); ++polyIt) {

        // returns 0 if false, +1 if true, -1 if point ON polygon boundary
        int res = ClipperLib::PointInPolygon(*it, *polyIt);

        if (res < 0) points_on_boundary++;
      }
    }
    // Consider intersection when more than one fault points lie on the polygon boundary.
    return (points_on_boundary > 1);
  }

  return 0;
}
//=====================================================================================================

void IntersectFaultsWithFaults(const Objects &faults, FaultIntersectionList &faultIntersectionList,
                               double distanceBuffer)
{
  std::cout << "Intersecting faults with faults..." << std::endl;

  for (std::size_t i = 0; i < faults.size(); i++) {
    for (std::size_t j = i + 1; j < faults.size(); j++) {

      double intersectionAngle;
      int intersectionType;

      if (FaultsAreIntersecting(faults[i], faults[j], intersectionType, intersectionAngle, distanceBuffer)) {

        int faultId1, faultId2;

        if (intersectionType == FAULT1_STOPS_ON_FAULT2 || intersectionType == FAULT2_STOPS_ON_FAULT1) {
          if (intersectionType == FAULT2_STOPS_ON_FAULT1) {
            faultId1 = faults[i].id;
            faultId2 = faults[j].id;
          } else {
            faultId1 = faults[j].id;
            faultId2 = faults[i].id;
          }
          // Adding one intersection: for faultId1 adding the list of faults that stop on it.
          std::pair<std::string, double> intersectionData = std::make_pair("T", intersectionAngle);
          faultIntersectionList[faultId1].push_back(std::make_pair(faultId2, intersectionData));

        } else if (intersectionType == FAULTS_CROSSING) {
          faultId1 = faults[i].id;
          faultId2 = faults[j].id;

          // Adding both intersections (A vs B, and B vs A).
          std::pair<std::string, double> intersectionData = std::make_pair("X", intersectionAngle);
          faultIntersectionList[faultId1].push_back(std::make_pair(faultId2, intersectionData));
          faultIntersectionList[faultId2].push_back(std::make_pair(faultId1, intersectionData));

        } else {
          throw std::runtime_error("Unknown intersection type!");
        }
      }
    }
  }
}
//=======================================================================================

bool FaultsAreIntersecting(const Object &fault1, const Object &fault2,
                           int &intersectionType, double &intersectionAngle,
                           double distanceBuffer)
{
  // Bounding box check.
  if (!AABB::BoundingBoxesOverlap(fault1.aabb, fault2.aabb)) {
    intersectionAngle = 0.;
    intersectionType = FAULTS_DONT_INTERSECT;
    return false;
  }

  {
    //------------------------------------------------------------------------------------
    // 1. Checking for T-type intersections (touching in one point / faultA stops on faultB).
    Path faultA, faultB;

    for (int f = 1; f <= 2; f++) {

      if (f == 1) {
        faultA = fault1.paths[0];
        faultB = fault2.paths[0];
      } else {
        faultA = fault2.paths[0];
        faultB = fault1.paths[0];
      }

      assert(faultA.size() > 1);
      assert(faultB.size() > 1);

      // Distance buffer from a point to segment (to find if a fault stops on another fault).
      const double distanceBufferSquared = distanceBuffer * distanceBuffer;

      // Loop over all segments of faultB.
      for (Path::const_iterator it = faultB.begin();
           it != faultB.end() - 1; ++it) {

        IntPoint p1 = *it;
        IntPoint p2 = *(it + 1);

        bool intersect = false;
        TEdge2 faultASegment;

        if (ConverterUtils::DistanceFromSegmentSqrd(faultA.front(), p1, p2) <= distanceBufferSquared) {
          intersect = true;
          // Form an edge from the first two points.
          faultASegment = TEdge2(*(faultA.begin() + 1), *(faultA.begin()));

        } else if (ConverterUtils::DistanceFromSegmentSqrd(faultA.back(), p1, p2) <= distanceBufferSquared) {
          intersect = true;
          // Form an edge from the last two points.
          faultASegment = TEdge2(*(faultA.end() - 2), *(faultA.end() - 1));
        }

        if (intersect) {
        // faultA stops on faultB.

          TEdge2 faultBSegment = TEdge2(*it, *(it + 1));

          intersectionAngle = ConverterUtils::CalculateEdgesAngle(faultASegment, faultBSegment);
          intersectionAngle = ConverterUtils::NormalizeAngle90(intersectionAngle);

          if (f == 1) {
            intersectionType = FAULT1_STOPS_ON_FAULT2;
          } else {
            intersectionType = FAULT2_STOPS_ON_FAULT1;
          }
          return true;
        }
      }
    }
  }

  {
    //------------------------------------------------------------------------------------
    // 2. Use clipper to intersect two faults.
    // Note: this intersection normally does not include T-type intersections (fault A stops on fault B).
    // Note: This type of intersection should not be present in a perfect map.

    // Adding data to the clipper.
    ClipperLib::Clipper clpr;
    clpr.Clear();
    clpr.AddPaths(fault1.paths, ClipperLib::ptSubject, false);
    clpr.AddPaths(fault2.paths, ClipperLib::ptSubject, false);

    ClipperLib::PolyTree polytree;

    // Do the intersection.
    if (!clpr.Execute(ClipperLib::ctIntersection, polytree, ClipperLib::pftEvenOdd, ClipperLib::pftEvenOdd)) {
      throw std::runtime_error("Error in FaultAndFaultIntersecting: Execute!");
    }

    if (clpr.intersections.size() > 0) {
      // Consider the first intersection (if there are more than one).
      ClipperLib::Intersection intersection = clpr.intersections[0];
      intersectionAngle = ConverterUtils::CalculateEdgesAngle(intersection.FaultEdge, intersection.PolygonEdge);
      intersectionAngle = ConverterUtils::NormalizeAngle90(intersectionAngle);

      intersectionType = FAULTS_CROSSING;
      return true;
    }
  }

  intersectionAngle = 0.;
  intersectionType = FAULTS_DONT_INTERSECT;
  return false;

}
//=======================================================================================

std::size_t IntersectContactWithFault(Contact &contact, const Object &fault,
                                      double angleEpsilon, double distanceEpsilon)
{
  // Ignore contacts & faults of less than two points.
  if (contact.path.size() < 2) { return false; }
  if (fault.paths[0].size() < 2) { return false; }

  const double distanceEpsilon_squared = distanceEpsilon * distanceEpsilon;

  std::size_t num_intersections = 0;

  // Loop over all segments of a contact.
  for (Path::iterator contact_vertex_it = contact.path.begin();
       contact_vertex_it != contact.path.end() - 1; ++contact_vertex_it)
  {
    // Form a contact segment.
    TEdge2 contactSegment = TEdge2(*contact_vertex_it, *(contact_vertex_it + 1));

    // Optimization: skip segments whose nodes are not lying inside of a bounding box of a fault.
    if (!AABB::PointInsideBox(contactSegment.Top, fault.aabb)
        && !AABB::PointInsideBox(contactSegment.Bot, fault.aabb))
    {
      continue;
    }

    bool coincide = false;

    // Loop over all segments of a fault.
    for (Path::const_iterator fault_vertex_it = fault.paths[0].begin();
         fault_vertex_it != fault.paths[0].end() - 1; ++fault_vertex_it)
    {
      // Form a fault segment.
      TEdge2 faultSegment = TEdge2(*fault_vertex_it, *(fault_vertex_it + 1));

      // Contact segment: C1-C2.
      // Fault segment:   F1-F2.
      TEdge2 segment_C1_C2, segment_F1_F2;
      bool intersect = true;

      // Build two segments C1-C2 and F1-F2 so that C1 coincides with F1 (within epsilon).
      if (ConverterUtils::GetDistanceSquared(contactSegment.Top, faultSegment.Top) <= distanceEpsilon_squared)
      {
        segment_C1_C2 = contactSegment;
        segment_F1_F2 = faultSegment;
      }
      else if (ConverterUtils::GetDistanceSquared(contactSegment.Top, faultSegment.Bot) <= distanceEpsilon_squared)
      {
        segment_C1_C2 = contactSegment;
        segment_F1_F2 = faultSegment.Inverse();
      }
      else if (ConverterUtils::GetDistanceSquared(contactSegment.Bot, faultSegment.Top) <= distanceEpsilon_squared)
      {
        segment_C1_C2 = contactSegment.Inverse();
        segment_F1_F2 = faultSegment;
      }
      else if (ConverterUtils::GetDistanceSquared(contactSegment.Bot, faultSegment.Bot) <= distanceEpsilon_squared)
      {
        segment_C1_C2 = contactSegment.Inverse();
        segment_F1_F2 = faultSegment.Inverse();
      }
      else
      {
        intersect = false;
      }

      if (intersect)
      {
        // Calculate the angle between segments.
        double angle = ConverterUtils::CalculateEdgesAngle(segment_C1_C2, segment_F1_F2);

        if ((std::fabs(angle) < angleEpsilon) || (std::fabs(angle) > 360.0 - angleEpsilon))
        {
          // Two segments intersect and are almost parallel.
          // Consider that a fault and a contact segment coincide.
          coincide = true;
          break;
        }
      }
    }

    if (coincide)
    {
      // Mark a corresponding segment vertex as intersecting.
      contact_vertex_it->type = 1;
      contact_vertex_it->faultID = fault.id;
      num_intersections++;
    }
  }

  return num_intersections;
}
//=======================================================================================

int IntersectPolygons(const Object &poly1, const Object &poly2, Contacts &contacts)
{
  // 1. Identifying common edges (those that are present in both multi-polygons).

    Object poly = poly1;

    size_t numCommonEdges = 0;

    // Loop over all polygons within the first multi-polygon.
    for (Paths::iterator path1 = poly.paths.begin(); path1 != poly.paths.end(); path1++) {
        // Loop over all edges of a polygon.
        for (Path::iterator vertex1 = path1->begin(); vertex1 != path1->end() - 1; vertex1++) {
            const IntPoint& v11 = *vertex1;
            const IntPoint& v12 = *(vertex1 + 1);

            bool found = false;

            // Loop over all polygons within the second multi-polygon.
            for (Paths::const_iterator path2 = poly2.paths.begin(); path2 != poly2.paths.end(); path2++) {
                // Loop over all edges of a polygon.
                for (Path::const_iterator vertex2 = path2->begin(); vertex2 != path2->end() - 1; vertex2++) {
                    const IntPoint& v21 = *vertex2;
                    const IntPoint& v22 = *(vertex2 + 1);

                    // Compare two edges (exact map coordinates case).
//                    if ((v11 == v21 && v12 == v22)
//                        || (v11 == v22 && v12 == v21)) {
//                        found = true;
//                        break;
//                    }

                    //---------------------------------------------------------------------------------
                    // Patch for non-exact maps.
                    // Consider vertexes between polygons to coincide if they located within epsilon distance.
                    double eps2 = 0.1; // epsilon distance (squared), in meters.

                    // Compare two edges: (v11, v12) vs (v21, v22).
                    double dist2_11_21 = ConverterUtils::GetDistanceSquared(v11, v21);
                    if (dist2_11_21 <= eps2) {
                        double dist2_12_22 = ConverterUtils::GetDistanceSquared(v12, v22);
                        if (dist2_12_22 <= eps2) {
                          found = true;
                          break;
                        }
                    }
                    // Compare two edges: (v11, v12) vs (v22, v21).
                    double dist2_11_22 = ConverterUtils::GetDistanceSquared(v11, v22);
                    if (dist2_11_22 <= eps2) {
                        double dist2_12_21 = ConverterUtils::GetDistanceSquared(v12, v21);
                        if (dist2_12_21 <= eps2) {
                          found = true;
                          break;
                        }
                    }
                    //---------------------------------------------------------------------------------
                }
                if (found) break;
            }

            if (!found) {
                // Mark the edge as not common.
                vertex1->type = - 1;
            } else {
                // Mark the edge as common.
                vertex1->type = 0;
                numCommonEdges++;
            }
        }
    }

  // 2. Building contacts (polylines).
  // Note that two polygons can have more than one contact.

  Contact contact;

  contact.obj1 = &poly1;
  contact.obj2 = &poly2;

  size_t numContactEdges = 0;
  std::size_t numContacts0 = contacts.size();

  // Loop over all polygons within the first multi-polygon.
  for (Paths::const_iterator path = poly.paths.begin(); path != poly.paths.end(); path++)
  {
    // Loop over all edges of a polygon.
    for (Path::const_iterator vertex = path->begin(); vertex != path->end() - 1; vertex++)
    {
      if (vertex->type != - 1)
      {
        contact.path.push_back(*vertex);
      }

      if (vertex->type == - 1 && contact.path.size() > 0)
      {
        // Close the edge.
        contact.path.push_back(*vertex);
        // Add a new contact.
        contacts.push_back(contact);

        numContactEdges += contact.path.size() - 1;
        contact.path.clear();
      }
    }

    if (contact.path.size() > 0)
    {// Last contact is not added, so adding it now.
      bool joined = false;
      if (contacts.size() > numContacts0)
      {
        // Note: path->back() != contact.path.back(), as the above loop iterates until path->end() - 1.
        if (path->back() == contacts[numContacts0].path.front())
        {
        // Joining last and first contacts as they are the same contact split between first and last polygon vertex.
          joined = true;

          contacts[numContacts0].path.insert(contacts[numContacts0].path.begin(),
                                             contact.path.begin(), contact.path.end());

          numContactEdges += contact.path.size();
        }
      }

      if (!joined)
      {
        // Close the last edge.
        contact.path.push_back(path->back());
        // Add a new contact.
        contacts.push_back(contact);

        numContactEdges += contact.path.size() - 1;
      }

      contact.path.clear();
    }
  }

  // Sanity check.
  if (numCommonEdges != numContactEdges)
  {
    std::cout << "Sanity check failed in IntersectPolygons!" << " numCommonEdges = "
              << numCommonEdges << " numContactEdges = " << numContactEdges
              << " for polygons " << poly1.id << " and " << poly2.id << std::endl;
    exit(0);
  }

  return (int)(numCommonEdges);
}
//========================================================================================
}

