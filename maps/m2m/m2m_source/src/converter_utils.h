/******************************************************************************
* A class with utils for ConverterLib.
*
* Author: Vitaliy Ogarko (2016).
*******************************************************************************/

#ifndef converter_utils_h
#define converter_utils_h

#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "clipper.hpp"

using ClipperLib::TEdge2;
using ClipperLib::cInt;
using ClipperLib::IntPoint;

namespace ConverterLib {

//! Converts integer to string.
std::string SSTR(int x);

class ConverterUtils
{
public:
  //! Converts latitude/longitude coordinate to an integer.
  static cInt CoordToInteger(double coord);

  //! Converts latitude/longitude coordinate to a double.
  static double CoordToDouble(cInt coord);

  //! Converts degrees to radians.
  static double ToRadians(double degrees);

  //! Converts radians to degrees.
  static double ToDegrees(double radians);

  //! Calculates distance in meters between two points given in latitude & longitude, using Haversine formula.
  static double LatLongDistance(double lat1, double lng1, double lat2, double lng2);

  //! Calculates the distance in meters between two IntPoint points in meters.
  static double GetDistance(const IntPoint &p1, const IntPoint &p2);

  //! Calculates the distance squared in meters between two IntPoint points in meters.
  static double GetDistanceSquared(const IntPoint &p1, const IntPoint &p2);

  // Calculates the distance from point to line defined by two points.
  static double DistanceFromLineSqrd(const IntPoint &pt, const IntPoint &ln1, const IntPoint &ln2);

  // Calculates the distance from point to segment defined by two points.
  static double DistanceFromSegmentSqrd(const IntPoint &pt, const IntPoint &ln1, const IntPoint &ln2);

  //! Calculates the distance in meters between two IntPoint points in lat/long.
  static double GetDistanceLatLong(const IntPoint &p1, const IntPoint &p2);

  //! Returns angle between two edges.
  static double CalculateEdgesAngle(TEdge2 &e1, TEdge2 &e2);

  //! Returns normalized angle to [0...90] range.
  static double NormalizeAngle90(double angle);

};

}

#endif
