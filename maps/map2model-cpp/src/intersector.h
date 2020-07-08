/******************************************************************************
* A class for intersecting faults & polygons.
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#ifndef intersector_h
#define intersector_h

#include <list>

#include "converter_types.h"
#include "converter.h"

namespace ConverterLib {

const int FAULTS_DONT_INTERSECT     = 0;
const int FAULTS_CROSSING           = 1;
const int FAULT1_STOPS_ON_FAULT2    = 2;
const int FAULT2_STOPS_ON_FAULT1    = 3;

// Maps Unit name to list of fault IDs.
typedef std::map<std::string, std::list<int> > UnitFaultIntersectionList;

// Intersecting fault ID, and intersection type & angle.
typedef std::pair<int, std::pair<std::string, double> > FaultIntersectionData;

// Maps Fault ID to list of intersecting fault IDs, and intersection types & angles.
typedef std::map<int, std::list<FaultIntersectionData> > FaultIntersectionList;

// Identifies the polygon ID where a point is located.
void IntersectPointsWithPolygons(PointPolygonIntersectionList &pointPolygonIntersectionList,
                                 const Objects &points, const Objects &polygons);

// Find the closest contacts for points.
// Note: faults are also considered as contacts, even if it is a fault inside a unit (i.e, between A and A).
// Our conventional contacts don't contain A-A fault contacts, so we test against faults separately.
void FindClosestContactForPoints(PointPolygonIntersectionList &pointPolygonIntersectionList,
                                 const Contacts& contacts, const Objects& faults,
                                 const Converter& converter, double pointToContactDistanceBuffer);

// Intersects every multipolygon (i.e., polygons with the same id) with every fault.
void IntersectFaultsWithPolygons(const Objects &faults, const Objects &polygons,
                                 UnitFaultIntersectionList &unitFaultIntersectionList);

// Intersects faults.
void IntersectFaultsWithFaults(const Objects &faults, FaultIntersectionList &faultIntersectionList,
                               double distanceBuffer);

// Intersects a fault with a polygon.
// Returns whether they are intersecting.
bool FaultAndPolygonIntersecting(const Object &fault, const Object &polygon, double distanceEpsilon = 0.);

// Intersects a fault with a fault.
// Returns intersection type, and intersection angle.
bool FaultsAreIntersecting(const Object &fault1, const Object &fault2,
                           int &intersectionType, double &intersectionAngle,
                           double distanceBuffer = 1.e-17);

// Check if a contact (partially) coincides with a fault, and mark intersecting contact segments.
// Returns the number of coinciding segments with a given fault.
std::size_t IntersectContactWithFault(Contact &contact, const Object &fault,
                                      double angleEpsilon, double distanceEpsilon);

// Finds common parts of boundaries between two multi-polygons.
// Returns the number of common edges.
int IntersectPolygons(const Object &poly1, const Object &poly2, Contacts &contacts, double distanceBuffer);

};

#endif
