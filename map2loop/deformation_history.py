import pandas
import numpy


class DeformationHistory(object):
    def __init__(self):
        self.history = []

        # Create empty fault and fold dataframes
        faultColumns = numpy.dtype(
            [
                ("faultId", int),
                ("name", str),
                ("minAge", float),
                ("maxAge", float),
                ("group", str),
                ("avgDisplacement", float),
                ("avgDownthrowDir", float),
                ("influenceDirection", float),
                ("verticalRadius", float),
                ("horizontalRadius", float),
                ("colour", str),
                ("centreX", float),
                ("centreY", float),
                ("centreZ", float),
                ("avgSlipDirX", float),
                ("avgSlipDirY", float),
                ("avgSlipDirZ", float),
                ("avgNormalX", float),
                ("avgNormalY", float),
                ("avgNormalZ", float),
            ]
        )
        self.faults = pandas.DataFrame(numpy.empty(0, dtype=faultColumns))
        self.faults = self.faults.set_index("name")

        foldColumns = numpy.dtype(
            [
                ("faultId", int),
                ("name", str),
                ("minAge", float),
                ("maxAge", float),
                ("periodic", bool),
                ("wavelength", float),
                ("amplitude", float),
                ("asymmetry", bool),
                ("asymmetryShift", float),
                ("secondaryWavelength", float),
                ("secondaryAmplitude", float),
            ]
        )
        self.folds = pandas.DataFrame(numpy.empty(0, dtype=foldColumns))
        self.folds = self.folds.set_index("name")

    def findfault(self, id):
        if type(id) == int:
            return self.faults[self.faults["faultId"] == id]
        elif type(id) == str:
            return self.faults[self.faults.index == id]
        else:
            print("ERROR: Unknown identifier type used to find fault")

    def findfold(self, id):
        if type(id) == int:
            return self.folds[self.folds["foldId"] == id]
        elif type(id) == str:
            return self.folds[self.folds.index == id]
        else:
            print("ERROR: Unknown identifier type used to find fold")

    def addFault(self, fault):
        if type(fault) == pandas.DataFrame or type(fault) == dict:
            if "name" in fault.keys():
                if fault["name"] in self.faults.index:
                    print("Replacing fault", fault["name"])
                self.faults[fault["name"]] = fault
            else:
                print("No name field in fault", fault)
        else:
            print("Cannot add fault to dataframe with type", type(fault))

    def addFold(self, fold):
        if type(fold) == pandas.DataFrame or type(fold) == dict:
            if "name" in fold.keys():
                if fold["name"] in self.folds.index:
                    print("Replacing fold", fold["name"])
                self.folds[fold["name"]] = fold
            else:
                print("No name field in fold", fold)
        else:
            print("Cannot add fold to dataframe with type", type(fold))
