/******************************************************************************
* tester.h
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#ifndef tester_h
#define tester_h

#include <cassert>

#include "intersector.h"

using namespace ClipperLib;

namespace ConverterLib {

#define _unused(x) ((void)(x))

// Testing intersection FAULTS vs FAULT
void Test_FaultsAreIntersecting()
{
  //--------------------------------------------------------------
  // Test setup
  //--------------------------------------------------------------
  double intersectionAngle;
  int intersectionType;

  _unused(intersectionAngle);
  _unused(intersectionType);

  Path path_fault1;
  path_fault1.push_back(IntPoint(0, 0));
  path_fault1.push_back(IntPoint(10, 0));
  path_fault1.push_back(IntPoint(10, 10));

  Object fault1, fault2;
  fault1.id = 1;
  fault2.id = 2;

  fault1.paths.push_back(path_fault1);

  Path path_fault2;

  //-----------------------------------------------
  // Crossing vertical edge.
  path_fault2.push_back(IntPoint(0, 5));
  path_fault2.push_back(IntPoint(20, 5));
  fault2.paths.push_back(path_fault2);

  assert(FaultsAreIntersecting(fault1, fault2, intersectionType, intersectionAngle));
  assert(intersectionType == FAULTS_CROSSING);
  assert(intersectionAngle == 90.);

  path_fault2.clear();
  fault2.paths.clear();

  //-----------------------------------------------
  // Crossing horizontal edge.
  path_fault2.push_back(IntPoint(5, -5));
  path_fault2.push_back(IntPoint(5, 5));
  fault2.paths.push_back(path_fault2);

  assert(FaultsAreIntersecting(fault1, fault2, intersectionType, intersectionAngle));
  assert(intersectionType == FAULTS_CROSSING);
  assert(intersectionAngle == 90.);

  path_fault2.clear();
  fault2.paths.clear();

  //-----------------------------------------------
  // Crossing two edges (straight).
  path_fault2.push_back(IntPoint(3, -3));
  path_fault2.push_back(IntPoint(13, 7));
  fault2.paths.push_back(path_fault2);

  assert(FaultsAreIntersecting(fault1, fault2, intersectionType, intersectionAngle));
  assert(intersectionType == FAULTS_CROSSING);
  assert(intersectionAngle == 45.);

  path_fault2.clear();
  fault2.paths.clear();

  //-----------------------------------------------
  // Crossing two edges (non straight).
  path_fault2.push_back(IntPoint(5, -5));
  path_fault2.push_back(IntPoint(5, 5));
  path_fault2.push_back(IntPoint(15, 5));
  fault2.paths.push_back(path_fault2);

  assert(FaultsAreIntersecting(fault1, fault2, intersectionType, intersectionAngle));
  assert(intersectionType == FAULTS_CROSSING);
  assert(intersectionAngle == 90.);

  path_fault2.clear();
  fault2.paths.clear();

  //-----------------------------------------------
  // Touching at one point (middle of the edge).
  path_fault2.push_back(IntPoint(5, -5));
  path_fault2.push_back(IntPoint(5, 0));
  fault2.paths.push_back(path_fault2);

  assert(FaultsAreIntersecting(fault1, fault2, intersectionType, intersectionAngle));
  assert(intersectionType == FAULT2_STOPS_ON_FAULT1);
  assert(intersectionAngle == 90.);

  path_fault2.clear();
  fault2.paths.clear();

  //-----------------------------------------------
  // Touching at one point (at the edge vertex).
  path_fault2.push_back(IntPoint(5, -5));
  path_fault2.push_back(IntPoint(0, 0));
  fault2.paths.push_back(path_fault2);

  assert(FaultsAreIntersecting(fault1, fault2, intersectionType, intersectionAngle));
  assert(intersectionType == FAULT1_STOPS_ON_FAULT2);
  assert(intersectionAngle == 45.);

  path_fault2.clear();
  fault2.paths.clear();

  //-----------------------------------------------
  // Touching at one point.
  path_fault2.push_back(IntPoint(0, -5));
  path_fault2.push_back(IntPoint(0, 5));
  fault2.paths.push_back(path_fault2);

  assert(FaultsAreIntersecting(fault1, fault2, intersectionType, intersectionAngle));
  assert(intersectionType == FAULT1_STOPS_ON_FAULT2);
  assert(intersectionAngle == 90.);

  path_fault2.clear();
  fault2.paths.clear();

  //-----------------------------------------------
  // Touching at one point.
  path_fault2.push_back(IntPoint(-2, -1));
  path_fault2.push_back(IntPoint(2, 1));
  fault2.paths.push_back(path_fault2);

  assert(FaultsAreIntersecting(fault1, fault2, intersectionType, intersectionAngle));
  assert(intersectionType == FAULT1_STOPS_ON_FAULT2);
  assert(intersectionAngle > 26.56 && intersectionAngle < 26.57); // arctan 1/2

  path_fault2.clear();
  fault2.paths.clear();

  //-----------------------------------------------
  // Touching at one point.
  path_fault2.push_back(IntPoint(2, 0));
  path_fault2.push_back(IntPoint(4, 1));
  fault2.paths.push_back(path_fault2);

  assert(FaultsAreIntersecting(fault1, fault2, intersectionType, intersectionAngle));
  assert(intersectionType == FAULT2_STOPS_ON_FAULT1);
  assert(intersectionAngle > 26.56 && intersectionAngle < 26.57); // arctan 1/2

  path_fault2.clear();
  fault2.paths.clear();

  //-----------------------------------------------
  // Crossing zero edges.
  path_fault2.push_back(IntPoint(0, 5));
  path_fault2.push_back(IntPoint(9, 5));
  fault2.paths.push_back(path_fault2);

  assert(!FaultsAreIntersecting(fault1, fault2, intersectionType, intersectionAngle));
  assert(intersectionType == FAULTS_DONT_INTERSECT);

  path_fault2.clear();
  fault2.paths.clear();
}

// Testing intersection FAULTS vs POLYGONS
void Test_FaultAndPolygonIntersecting()
{
  //--------------------------------------------------------------
  // Test setup
  //--------------------------------------------------------------
  Path path_poly;
  path_poly.push_back(IntPoint(0, 0));
  path_poly.push_back(IntPoint(10, 0));
  path_poly.push_back(IntPoint(10, 10));
  path_poly.push_back(IntPoint(0, 10));
  path_poly.push_back(IntPoint(0, 0));

  Object fault, polygon;
  polygon.paths.push_back(path_poly);

  Path path_fault;

  //-----------------------------------------------
  // Crossing 2 edges.
  path_fault.push_back(IntPoint(5, -10));
  path_fault.push_back(IntPoint(5, +20));
  fault.paths.push_back(path_fault);

  assert(FaultAndPolygonIntersecting(fault, polygon));

  path_fault.clear();
  fault.paths.clear();

  //-----------------------------------------------
  // Crossing 1 edge (top).
  path_fault.push_back(IntPoint(5, 2));
  path_fault.push_back(IntPoint(5, +20));
  fault.paths.push_back(path_fault);

  assert(FaultAndPolygonIntersecting(fault, polygon));

  path_fault.clear();
  fault.paths.clear();

  //-----------------------------------------------
  // Crossing 0 edges (inside a polygon).
  path_fault.push_back(IntPoint(5, 2));
  path_fault.push_back(IntPoint(5, 8));
  fault.paths.push_back(path_fault);

  assert(FaultAndPolygonIntersecting(fault, polygon));

  path_fault.clear();
  fault.paths.clear();

  //-----------------------------------------------
  // Crossing 0 edges (outside a polygon).
  path_fault.push_back(IntPoint(15, 2));
  path_fault.push_back(IntPoint(15, 8));
  fault.paths.push_back(path_fault);

  assert(!FaultAndPolygonIntersecting(fault, polygon));

  path_fault.clear();
  fault.paths.clear();

  //-----------------------------------------------
  // Crossing 1 edge (left).
  path_fault.push_back(IntPoint(-1, 1));
  path_fault.push_back(IntPoint(1, 9));
  fault.paths.push_back(path_fault);

  assert(FaultAndPolygonIntersecting(fault, polygon));

  path_fault.clear();
  fault.paths.clear();

  //-----------------------------------------------
  // Coincides with polygon edge (left).
  path_fault.push_back(IntPoint(0, 10));
  path_fault.push_back(IntPoint(0, 0));
  fault.paths.push_back(path_fault);

  assert(FaultAndPolygonIntersecting(fault, polygon));

  path_fault.clear();
  fault.paths.clear();

  //-----------------------------------------------
  // Partially coincides with polygon edge (one point).
  path_fault.push_back(IntPoint(0, 5));
  path_fault.push_back(IntPoint(0, 15));
  fault.paths.push_back(path_fault);

  assert(!FaultAndPolygonIntersecting(fault, polygon));

  path_fault.clear();
  fault.paths.clear();

  //-----------------------------------------------
  // Partially coincides with polygon edge (two points).
  path_fault.push_back(IntPoint(0, 5));
  path_fault.push_back(IntPoint(0, 7));
  path_fault.push_back(IntPoint(0, 15));
  fault.paths.push_back(path_fault);

  assert(FaultAndPolygonIntersecting(fault, polygon));

  path_fault.clear();
  fault.paths.clear();

  //-----------------------------------------------
  // Touches polygon edge (bottom).
  path_fault.push_back(IntPoint(5, -10));
  path_fault.push_back(IntPoint(5, 0));
  fault.paths.push_back(path_fault);

  assert(!FaultAndPolygonIntersecting(fault, polygon));

  path_fault.clear();
  fault.paths.clear();

}

}

#endif
