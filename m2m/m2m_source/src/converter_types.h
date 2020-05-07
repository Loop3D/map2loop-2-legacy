/******************************************************************************
* Different types used in ConverterLib.
*
* Author: Vitaliy Ogarko (2016).
*******************************************************************************/

#ifndef converter_types_h
#define converter_types_h

#include <vector>
#include <map>
#include <string>

#include "converter_utils.h"
#include "bounding_box.h"
#include "clipper.hpp"

using ClipperLib::Paths;

namespace ConverterLib {

// Fault or multi-polygon or point.
class Object
{
public:
  // Unique id.
  int id;
  // Min/max age ma (for polygons).
  double min_age, max_age;
  // Unit name (for polygons).
  std::string name;
  // Group name (for polygons).
  std::string group;
  // ENGUINT code (for polygons).
  std::string code;
  // Description (for polygons).
  std::string description;
  // Rocktypes (for polygons).
  std::string rocktype1, rocktype2;

  // Site code (for deposits).
  std::string site_code;
  // Site commodity (for deposits).
  std::string site_commo;

  // Polygons withing a multi-polygon, or a fault, or point coordinates.
  Paths paths;
  // Bounding box (AABB), Y-axis points down.
  AABB aabb;

  Object(AABB _aabb = AABB(0)):
         aabb(_aabb) {};
};

enum ContactType { StratigraphicContact, FaultContact, MixedContact, NotSpecifiedContact,
                   IntrusiveIgneousContact, IgneousContact};

// Litho contact between two polygons.
class Contact
{
public:
  // Unique contact id.
  int id;
  // Pointers to polygons that share a contact.
  const Object *obj1, *obj2;
  // Coordinates of vertexes.
  Path path;
  // The contact length (in meters).
  double length;
  // Contact type.
  ContactType type;

  // For fault contact type.
  int faultID;

  Contact(int _id = - 1):
          id(_id) {};

  bool operator< (const Contact &other) const
  {
    return id < other.id;
  }
};

struct PointInPolygon {
    const Object* point;

    bool onContact;
    // Containing and neighbor polygons.
    const Object *containingPolygon, *neighbourPolygon;
    // Fault, Strat, etc.
    ContactType contactType;
    // Distance (squared) to the closest contact segment.
    double distanceToContact2;
    // For fault contacts
    int faultObjectID;
};
typedef std::vector<PointInPolygon> PointPolygonIntersectionList;

// Class representing a geological unit.
class Unit
{
public:
  // Unique id only within a local (clipped) map.
  int id;
  // Unit name (unique identifier).
  std::string name;
  // Group name.
  std::string group;
  // ENGUINT codes of polygons that share this unit.
  std::vector<std::string> codes;
  // Perimeter length (total length of all polygon segments that belong to this unit).
  double length;
  // Min/max age (ma).
  double min_age, max_age;
  // A flag for sill-units.
  bool is_sill;
  // Rocktypes.
  std::string rocktype1, rocktype2;
  // Deposits located on this unit.
  std::vector<const PointInPolygon*> deposits;

  bool operator< (const Unit &other) const
  {
    return length < other.length;
  }
};

typedef std::vector<Contact> Contacts;

// Litho contact between two units.
class UnitContact
{
public:
  // Units that are in contact.
  const Unit *unit1, *unit2;
  // Polygon contacts which the unit contact is made of.
  Contacts contacts;
  // Total contact length, length of fault, stratigraphic and mixed contacts.
  double total_length, faults_length, stratigraphic_length, mixed_length;
  // Relative contact length w.r.t. (average) unit length.
  double relative_length;
  // Contact type.
  ContactType type;
  // Flags for igneous contacts;
  bool isIgneous;
  bool isIntrusiveIgneous;
  // Deposits located on this contact.
  std::vector<const PointInPolygon*> deposits;
};

typedef std::vector<Object> Objects;
// Container of units in the clipped map (not all map units stored).
typedef std::map<std::string, Unit> Units;
typedef std::map<std::string, const Unit*> Units_p;
typedef std::vector<UnitContact> UnitContacts;
typedef std::map<std::string, int> Groups;

}

#endif
