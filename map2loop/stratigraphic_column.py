import pandas as pd
import numpy as np
import geopandas as gpd


class StratigraphicColumn(object):
    def __init__(self):
        self.column = []

        self.groups = []

        # Create empty dataframes for units
        self.stratigraphicUnitColumns = np.dtype(
            [
                ("layerId", int),
                ("name", str),
                ("minAge", float),
                ("maxAge", float),
                ("group", str),
                ("thickness", float),
                ("colour", str),
                ("indexInGroup", int),
                ("groupNum", int),
                ("NumInGroup", int),
            ]
        )
        self.stratigraphicUnits = pd.DataFrame(
            np.empty(0, dtype=self.stratigraphicUnitColumns)
        )
        self.stratigraphicUnits = self.stratigraphicUnits.set_index("name")

        self.lithologyUnitColumns = np.dtype(
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
        self.lithologyUnits = pd.DataFrame(np.empty(0, dtype=self.lithologyUnitColumns))
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
        if type(unit) == pd.DataFrame or type(unit) == dict:
            if "name" in unit.keys():
                if unit["name"] in self.stratigraphicUnits.index:
                    print("Replacing stratigraphic unit", unit["name"])
                self.stratigraphicUnits.loc[unit["name"]] = unit
            else:
                print("No name field in stratigraphic unit", unit)
        else:
            print("Cannot add unit to dataframe with type", type(unit))

    def addLithologyUnit(self, unit):
        if type(unit) == pd.DataFrame or type(unit) == dict:
            if "name" in unit.keys():
                if unit["name"] in self.lithologyUnits.index:
                    print("Replacing lithology unit", unit["name"])
                self.lithologyUnits.loc[unit["name"]] = unit
            else:
                print("No name field in lithology unit", unit)
        else:
            print("Cannot add unit to dataframe with type", type(unit))

    def populate(self, geology_map_data: gpd.GeoDataFrame):
        if geology_map_data.shape[0] == 0:
            return
        geology_data = geology_map_data.copy()
        geology_data = geology_data.drop_duplicates(subset=["UNIT_NAME"])
        geology_data = geology_data.reset_index(drop=True)
        # geology_data = geology_data.dropna(subset=["UNIT_NAME"])

        self.stratigraphicUnits = pd.DataFrame(
            np.empty(geology_data.shape[0], dtype=self.stratigraphicUnitColumns)
        )
        self.stratigraphicUnits["layerId"] = np.arange(geology_data.shape[0])
        self.stratigraphicUnits["name"] = geology_data["UNIT_NAME"]
        self.stratigraphicUnits["minAge"] = geology_data["MIN_AGE"]
        self.stratigraphicUnits["maxAge"] = geology_data["MAX_AGE"]
        self.stratigraphicUnits["group"] = geology_data["GROUP"]
        self.stratigraphicUnits["thickness"] = -1
        self.stratigraphicUnits["colour"] = ""
        self.stratigraphicUnits["indexInGroup"] = -1

    def set_stratigraphic_unit_parameter_by_name(
        self, name: str, parameter: str, value
    ):
        self.stratigraphicUnits.iloc[self.stratigraphicUnits["name"] == name][
            parameter
        ] = value

    def sort_from_relationship_list(relationshipList: list):
        pass
