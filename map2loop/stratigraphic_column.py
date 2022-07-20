import pandas
import numpy


class StratigraphicColumn(object):
    def __init__(self):
        self.column = []

        self.groups = []

        # Create empty dataframes for units
        stratigraphicUnitColumns = numpy.dtype(
            [
                ("layerId", int),
                ("name", str),
                ("minAge", float),
                ("maxAge", float),
                ("group", str),
                ("thickness", float),
                ("colour", str),
            ]
        )
        self.stratigraphicUnits = pandas.DataFrame(
            numpy.empty(0, dtype=stratigraphicUnitColumns)
        )
        self.stratigraphicUnits = self.stratigraphicUnits.set_index("name")

        lithologyUnitColumns = numpy.dtype(
            [
                ("layerId", int),
                ("name", str),
                ("minAge", float),
                ("maxAge", float),
                ("group", str),
                ("thickness", float),
                ("colour", str),
            ]
        )
        self.lithologyUnits = pandas.DataFrame(
            numpy.empty(0, dtype=lithologyUnitColumns)
        )
        self.lithologyUnits = self.lithologyUnits.set_index("name")

    def findStratigraphicUnit(self, id):
        if type(id) == int:
            return self.stratigraphicUnits[self.stratigraphicUnits["layerId"] == id]
        elif type(id) == str:
            return self.stratigraphicUnits[self.stratigraphicUnits.index == id]
        else:
            print("ERROR: Unknown identifier type used to find stratigraphic unit")

    def findLithologyUnit(self, id):
        if type(id) == int:
            return self.lithologyUnits[self.lithologyUnits["layerId"] == id]
        elif type(id) == str:
            return self.lithologyUnits[self.lithologyUnits.index == id]
        else:
            print("ERROR: Unknown identifier type used to find lithology unit")

    def addStratigraphicUnit(self, unit):
        if type(unit) == pandas.DataFrame or type(unit) == dict:
            if "name" in unit.keys():
                if unit["name"] in self.stratigraphicUnits.index:
                    print("Replacing stratigraphic unit", unit["name"])
                self.stratigraphicUnits.loc[unit["name"]] = unit
            else:
                print("No name field in stratigraphic unit", unit)
        else:
            print("Cannot add unit to dataframe with type", type(unit))

    def addLithologyUnit(self, unit):
        if type(unit) == pandas.DataFrame or type(unit) == dict:
            if "name" in unit.keys():
                if unit["name"] in self.lithologyUnits.index:
                    print("Replacing lithology unit", unit["name"])
                self.lithologyUnits.loc[unit["name"]] = unit
            else:
                print("No name field in lithology unit", unit)
        else:
            print("Cannot add unit to dataframe with type", type(unit))
