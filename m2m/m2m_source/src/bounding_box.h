/******************************************************************************
* A class of bounding box used in the ConverterLib.
*
* Author: Vitaliy Ogarko (2016).
*******************************************************************************/

#ifndef bounding_box_h
#define bounding_box_h

#include "clipper.hpp"
#include "converter_utils.h"

using ClipperLib::IntPoint;
using ClipperLib::Path;
using ClipperLib::Paths;

namespace ConverterLib {

// Bounding box (axis-aligned).
class AABB
{
public:
  //! Coordinates of the corners.
  IntPoint TopLeft;
  IntPoint BtmRight;

  AABB(IntPoint _TopLeft = IntPoint(0, 0),
       IntPoint _BtmRight = IntPoint(0, 0)):
       TopLeft(_TopLeft), BtmRight(_BtmRight) {};

  AABB(int a = 0): TopLeft(IntPoint(a, a)), BtmRight(IntPoint(a, a)) {};

  AABB(double tx = 0., double ty = 0.,
       double bx = 0., double by = 0.):
       TopLeft(IntPoint(ConverterUtils::CoordToInteger(tx), ConverterUtils::CoordToInteger(ty))),
       BtmRight(IntPoint(ConverterUtils::CoordToInteger(bx), ConverterUtils::CoordToInteger(by))) {};

  bool operator==(const AABB &other) const {
     if (TopLeft == other.TopLeft && BtmRight == other.BtmRight) {
       return true;
     } else {
       return false;
     }
  }

  bool operator!=(const AABB &other) const {
    return !(*this == other);
  }

  //! Checks if bounding boxes do overlap.
  // Note that the +y axis points down,
  // and the +x axis is directed to the right (i.e. typical screen/pixel coordinates).
  static bool BoundingBoxesOverlap(const AABB &box1, const AABB &box2)
  {
    return !(box2.TopLeft.X > box1.BtmRight.X
             || box2.BtmRight.X < box1.TopLeft.X
             || box2.TopLeft.Y > box1.BtmRight.Y
             || box2.BtmRight.Y < box1.TopLeft.Y);
  }

  //! Converts the AABB to Path type.
  static Path AABB2Path(const AABB &box)
  {
    Path path;
    path.push_back(box.TopLeft);
    path.push_back(IntPoint(box.BtmRight.X, box.TopLeft.Y));
    path.push_back(box.BtmRight);
    path.push_back(IntPoint(box.TopLeft.X, box.BtmRight.Y));
    path.push_back(path.front());

    return path;
  }

  //! Calculates a bounding box for given paths.
  // Note that the +y axis points down,
  // and the +x axis is directed to the right (i.e. typical screen/pixel coordinates).
  static AABB CalculateBoundingBox(const Paths &paths)
  {
    IntPoint TopLeft  = paths[0][0];
    IntPoint BtmRight = paths[0][0];

    // Loop over all paths (e.g. polygons) within an object (e.g. multipolygon).
    for (std::size_t j = 0; j < paths.size(); j++)
    {
      // Loop over all vertices.
      for (std::size_t k = 0; k < paths[j].size(); k++)
      {
        IntPoint v = paths[j][k];

        if (v.X < TopLeft.X) TopLeft.X = v.X;
        if (v.Y < TopLeft.Y) TopLeft.Y = v.Y;

        if (v.X > BtmRight.X) BtmRight.X = v.X;
        if (v.Y > BtmRight.Y) BtmRight.Y = v.Y;
      }
    }

    return AABB(TopLeft, BtmRight);
  }

  //! Check if a point is inside a box.
  static bool PointInsideBox(const IntPoint &point, const AABB &box)
  {
    if (point.X >= box.TopLeft.X
        && point.X <= box.BtmRight.X
        && point.Y >= box.TopLeft.Y
        && point.Y <= box.BtmRight.Y)
    {
      return true;
    }
    return false;
  }
};

}

#endif
