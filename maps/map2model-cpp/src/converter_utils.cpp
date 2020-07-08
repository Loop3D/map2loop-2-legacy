/******************************************************************************
* A class with utils for ConverterLib.
*
* Author: Vitaliy Ogarko (2016).
*******************************************************************************/

#ifndef converter_utils_cpp
#define converter_utils_cpp

#include <cmath>
#include <sstream>

#include "converter_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ConverterLib {

//! Converts integer to string.
std::string SSTR(int x) {
    std::stringstream ss;
    ss << x;
    std::string str = ss.str();
    return str;
}

const double CONVERSION_PRECISION = 1e10;

cInt ConverterUtils::CoordToInteger(double coord)
{
  return (cInt)(coord * CONVERSION_PRECISION);
}

double ConverterUtils::CoordToDouble(cInt coord)
{
  return ((double)coord / CONVERSION_PRECISION);
}

double ConverterUtils::ToRadians(double degrees)
{
  return degrees * M_PI / 180.0;
}

double ConverterUtils::ToDegrees(double radians)
{
  return 180.0 * radians / M_PI;
}
//========================================================================================

// Tested on LatLongDistance(42., -82., 42., -83.), correct result is 82.63 km.
double ConverterUtils::LatLongDistance(double lat1, double lng1, double lat2, double lng2)
{
  double earthRadius = 6371000.0; // meters
  double dLat = ToRadians(lat2 - lat1);
  double dLng = ToRadians(lng2 - lng1);
  double a = sin(dLat / 2) * sin(dLat / 2) +
             cos(ToRadians(lat1)) * cos(ToRadians(lat2)) *
             sin(dLng / 2) * sin(dLng / 2);
  double c = 2 * atan2(sqrt(a), sqrt(1. - a));
  double dist = earthRadius * c;
  return dist;
}
//========================================================================================

double ConverterUtils::GetDistanceLatLong(const IntPoint &p1, const IntPoint &p2)
{
  double distance = LatLongDistance(CoordToDouble(p1.Y),
                                    CoordToDouble(p1.X),
                                    CoordToDouble(p2.Y),
                                    CoordToDouble(p2.X));

  return distance;
}
//========================================================================================

double ConverterUtils::GetDistanceSquared(const IntPoint &p1, const IntPoint &p2)
{
  double distance_squared = pow(CoordToDouble(p1.X) - CoordToDouble(p2.X), 2.0)
                            + pow(CoordToDouble(p1.Y) - CoordToDouble(p2.Y), 2.0);
  return distance_squared;
}
//========================================================================================

double ConverterUtils::DistanceFromLineSqrd(const IntPoint &pt, const IntPoint &ln1, const IntPoint &ln2) {
    double y1 = CoordToDouble(ln1.Y);
    double y2 = CoordToDouble(ln2.Y);
    double x1 = CoordToDouble(ln1.X);
    double x2 = CoordToDouble(ln2.X);

    double x0 = CoordToDouble(pt.X);
    double y0 = CoordToDouble(pt.Y);

    // See https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    double nom = (y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1;
    double denom2 = (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1);
    return nom * nom / denom2;
}
//========================================================================================

double ConverterUtils::DistanceFromSegmentSqrd(const IntPoint &pt, const IntPoint &ln1, const IntPoint &ln2) {
    double y1 = CoordToDouble(ln1.Y);
    double y2 = CoordToDouble(ln2.Y);
    double x1 = CoordToDouble(ln1.X);
    double x2 = CoordToDouble(ln2.X);

    double x0 = CoordToDouble(pt.X);
    double y0 = CoordToDouble(pt.Y);

    //------------------------------------------------------
    // See https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment

    double A = x0 - x1;
    double B = y0 - y1;
    double C = x2 - x1;
    double D = y2 - y1;

    double dot = A * C + B * D;
    double len_sq = C * C + D * D;
    double param = dot / len_sq;

    double xx, yy;

    if (param < 0 || (x1 == x2 && y1 == y2)) {
        xx = x1;
        yy = y1;
    }
    else if (param > 1) {
        xx = x2;
        yy = y2;
    }
    else {
        xx = x1 + param * C;
        yy = y1 + param * D;
    }

    double dx = x0 - xx;
    double dy = y0 - yy;

    return (dx * dx + dy * dy);
}
//========================================================================================

double ConverterUtils::GetDistance(const IntPoint &p1, const IntPoint &p2)
{
  double distance = sqrt(GetDistanceSquared(p1, p2));
  return distance;
}
//========================================================================================

double ConverterUtils::CalculateEdgesAngle(TEdge2 &e1, TEdge2 &e2)
{
  IntPoint vector1 = IntPoint(e1.Top.X - e1.Bot.X, e1.Top.Y - e1.Bot.Y);
  IntPoint vector2 = IntPoint(e2.Top.X - e2.Bot.X, e2.Top.Y - e2.Bot.Y);

  double angle = atan2(vector2.Y, vector2.X) - atan2(vector1.Y, vector1.X);

  // Normalize angle to the range 0 .. 2 * Pi.
  if (angle < 0) angle += 2 * M_PI;

  return ToDegrees(angle);
}
//========================================================================================

double ConverterUtils::NormalizeAngle90(double angle)
{
  double newAngle = angle;

  if (newAngle >= 180.) {
    newAngle -= 180.;
  }
  if (newAngle > 90.) {
    newAngle = 180. - newAngle;
  }
  return newAngle;
}

}

#endif

