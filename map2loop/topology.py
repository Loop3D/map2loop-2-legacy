import os
import re
from map2loop.m2l_enums import Datatype, VerboseLevel
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import numpy as np
from math import degrees, acos
import warnings
import beartype
from map2loop.m2l_interpolation import interpolate_contacts

from map2loop.stratigraphic_column import StratigraphicColumn

from . import m2l_utils
from .m2l_utils import display
from .config import Config


class Topology(object):
    @beartype.beartype
    def __init__(self, config: Config):
        self.asud_strat_file = "https://gist.githubusercontent.com/yohanderose/3b257dc768fafe5aaf70e64ae55e4c42/raw/8598c7563c1eea5c0cd1080f2c418dc975cc5433/ASUD.csv"
        if config.run_flags["aus"]:
            # Attempt to open asud database and abort if it fails
            try:
                asud = pd.read_csv(self.asud_strat_file, sep=",")
                asud["over"] = asud["over"].map(lambda x: x.lower())
                asud["under"] = asud["under"].map(lambda x: x.lower())
                asud = asud.replace("[ -/?]", "_", regex=True)
                self.asud = asud
            except Exception as e:
                warnings.warn(str(e))
                self.asud = None
        else:
            self.asud = None
        self.graph = nx.read_gml(config.strat_graph_filename, label="id")
        self.graph_path = config.graph_path

    ####################################
    # parse stratigraphy GML file to get number of series and series names
    #
    # get_series(path_in,id_label)
    # Args:
    # path_in path to gml file
    # id_label code for eacah node that defines node type (group or formation)
    ####################################
    def get_group_labels(self):
        return {
            x: self.graph.nodes[x]["LabelGraphics"]["text"]
            for x, y in self.graph.nodes(data=True)
            if "isGroup" in y
        }

    ####################################
    # parse stratigraphy GML file to save units for each series
    #
    # Save out csv of maximum and minimum ages of formations within each group so that groups can be better sorted.
    # save_units(G,tmp_path,glabels)
    # Args:
    # G networkx format stratigraphy graph
    # tmp_path directory of temporary outputs glabels list of group names
    #
    # The choice of what constitutes basic unit and what a group of units is defined in the c_l codes. Not even sure we need two levels but it seemed like a good idea at the time. Text outputs list alternate topologies for series and surfaces, which if confirmed by comapring max-min ages will be a nice source of uncertainty.
    ####################################
    @beartype.beartype
    def save_units(self, config: Config, stratColumn: StratigraphicColumn):
        if config.verbose_level != VerboseLevel.NONE:
            print("Generating topology graph display and unit groups...")

        group_labels = self.get_group_labels()

        # For each group create a file with the permutation of possible orderings
        for p in group_labels.keys():
            # Take a copy of the master graph
            GD = self.graph.subgraph(
                [
                    x
                    for x, y in self.graph.nodes(data=True)
                    if "gid" in y and y["gid"] == p
                ]
            ).copy()

            # Store node labels
            labels = {}
            for node in GD.nodes():
                labels[node] = self.graph.nodes[node]["LabelGraphics"]["text"]

            # Look for cycles and remove one edge to break it
            # For Australian data see if ASUD helps
            cycles = nx.simple_cycles(GD)

            if config.run_flags["aus"]:
                for cy in cycles:
                    found = False

                    len_cy = len(cy)
                    # For each edge in cycle
                    for i in range(len_cy - 1):
                        # Skip if either node is a group node
                        if (
                            "isGroup" in GD.nodes[cy[i]]
                            or "isGroup" in GD.nodes[cy[i + 1]]
                        ):
                            continue

                        # Search ASUD for specific edge
                        # Note this relies on direct match of label to ASUD labels
                        # to avoid removing a false negative
                        glabel_0 = GD.nodes[cy[i]]["LabelGraphics"]["text"].lower()
                        glabel_1 = GD.nodes[cy[i + 1]]["LabelGraphics"]["text"].lower()
                        edge = self.asud.loc[
                            (
                                (self.asud["over"] == glabel_0)
                                & (self.asud["under"] == glabel_1)
                            )
                            | (
                                (self.asud["under"] == glabel_0)
                                & (self.asud["over"] == glabel_1)
                            )
                        ]
                        # If not found and edge exists remove it, possible false negative
                        if len(edge) == 0:
                            if GD.has_edge(cy[i], cy[i + 1]):
                                warning_msg = (
                                    "map2loop warning 1: Stratigraphic relationship: "
                                    + str(GD.nodes[cy[i]]["LabelGraphics"]["text"])
                                    + " overlies "
                                    + str(GD.nodes[cy[i + 1]]["LabelGraphics"]["text"])
                                    + " removed to prevent cycle"
                                )
                                warnings.warn(warning_msg)
                                GD.remove_edge(cy[i], cy[i + 1])
                                found = True

                    if not found:
                        # Special case for edge between end element and first element in list
                        if (
                            "isGroup" in GD.nodes[len_cy - 1]
                            or "isGroup" in GD.nodes[0]
                        ):
                            continue
                        glabel_0 = GD.nodes[cy[len_cy - 1]]["LabelGraphics"][
                            "text"
                        ].lower()
                        glabel_1 = GD.nodes[cy[0]]["LabelGraphics"]["text"].lower()
                        edge = self.asud.loc[
                            (
                                (self.asud["over"] == glabel_0)
                                & (self.asud["under"] == glabel_1)
                            )
                            | (
                                (self.asud["under"] == glabel_0)
                                & (self.asud["over"] == glabel_1)
                            )
                        ]
                        if len(edge) == 0:
                            if GD.has_edge(cy[len_cy - 1], cy[0]):
                                warning_msg = (
                                    "map2loop warning 1: Stratigraphic relationship: "
                                    + str(
                                        GD.nodes[cy[len_cy - 1]]["LabelGraphics"][
                                            "text"
                                        ]
                                    )
                                    + " overlies "
                                    + str(GD.nodes[cy[0]]["LabelGraphics"]["text"])
                                    + " removed to prevent cycle"
                                )
                                warnings.warn(warning_msg)
                                GD.remove_edge(cy[len_cy - 1], cy[0])
                                found = True

                    if not found:
                        # As all links in cycle have been found (unlikely I know)
                        # Remove the first edge in the cycle
                        warning_msg = (
                            "map2loop warning 2: Stratigraphic relationship: "
                            + str(GD.nodes[cy[0]]["LabelGraphics"]["text"])
                            + " overlies "
                            + str(GD.nodes[cy[1]]["LabelGraphics"]["text"])
                            + " removed to prevent cycle"
                        )
                        warnings.warn(warning_msg)
                        if GD.has_edge(cy[len_cy - 1], cy[0]):
                            if GD.has_edge(cy[0], cy[1]):
                                GD.remove_edge(cy[0], cy[1])

            else:
                # No help from ASUD so just remove the first edge in cycle
                for cy in cycles:
                    warning_msg = (
                        "map2loop warning 3: Stratigraphic relationship: "
                        + str(GD.nodes[cy[0]]["LabelGraphics"]["text"])
                        + " overlies "
                        + str(GD.nodes[cy[1]]["LabelGraphics"]["text"])
                        + " removed to prevent cycle"
                    )
                    warnings.warn(warning_msg)
                    if GD.has_edge(cy[0], cy[1]):
                        GD.remove_edge(cy[0], cy[1])

            if config.verbose_level == VerboseLevel.ALL:
                plt.figure(p + 1)  # display strat graph for one group
                plt.title(group_labels[p])
                plt.tight_layout()
                nx.draw_networkx(
                    GD, pos=nx.kamada_kawai_layout(GD), arrows=True, with_labels=False
                )
                nx.draw_networkx_labels(
                    GD,
                    pos=nx.kamada_kawai_layout(GD),
                    labels=labels,
                    font_size=12,
                    font_family="sans-serif",
                )
                plt.show()

            # TODO: Use a class numpy array to store these values not files
            # one possible sorted directional graphs
            nlist = list(nx.topological_sort(GD))

            ordered_unit_labels = [
                self.graph.nodes[x]["LabelGraphics"]["text"] for x in nlist
            ]

            # Set ordered index in stratigraphic column based on ordered units
            for i in np.arange(len(ordered_unit_labels)):
                stratColumn.stratigraphicUnits.loc[
                    (stratColumn.stratigraphicUnits["name"] == ordered_unit_labels[i]),
                    "indexInGroup",
                ] = i

            # Output to file
            # "Choice 0", 1stUnit, 2ndUnit, 3rdUnit, etc
            ordered_unit_labels.insert(0, "Choice 0")
            data = np.array(ordered_unit_labels)
            group_filename = os.path.join(config.tmp_path, group_labels[p] + ".csv")
            np.savetxt(group_filename, data[None], fmt="%s", delimiter=",")

    ####################################
    # save out a list of max/min/ave ages of all formations in a group
    #
    # abs_age_groups(geol,tmp_path)
    # Args:
    # geol geopandas GeoDataFrame of geology polygons
    # tmp_path
    #
    # Save out csv of maximum and minimum ages of formations within each group so that groups can be better sorted.
    ####################################
    @beartype.beartype
    def abs_age_groups(self, geol: gpd.GeoDataFrame, tmp_path: str):
        # Extract unique group names from geology map
        groupList = geol["GROUP"].unique()
        groupList = groupList[groupList != "None"]

        # Extract min and max ages for each group based on all formations in that group
        groups = pd.DataFrame([], columns=["group_", "min", "max"])
        groups["group_"] = groupList
        groups["min"] = [min(geol[geol["GROUP"] == x]["MIN_AGE"]) for x in groupList]
        groups["max"] = [max(geol[geol["GROUP"] == x]["MAX_AGE"]) for x in groupList]

        # Extract formations without a group as a separate group in themselves
        if len(geol[geol["GROUP"] == "None"]):
            ages = geol[geol["GROUP"] == "None"]["UNIT_NAME", "MIN_AGE", "MAX_AGE"]
            ages.rename(
                columns={"UNIT_NAME": "group_", "MIN_AGE": "min", "MAX_AGE": "max"}
            )
            groups = pd.concat([ages, groups])

        # calculate average ages and sort by the average
        groups["ave"] = (groups["min"] + groups["max"]) / 2
        groups = groups.sort_values(by="ave")

        # re-index, output to file and return for internal use
        groups.reset_index(inplace=True, drop=True)
        groups.to_csv(os.path.join(tmp_path, "age_sorted_groups.csv"))

        return groups

    ####################################
    # save out tables of groups and sorted formation data
    #
    # save_group(G,tmp_path,glabels,geol_clip,c_l)
    # G networkx format stratigraphy graph
    # tmp_path directory of temporary outputs glabels list of group names
    # geol_clip path to clipped geology layer c_l dictionary of codes and labels specific to input geo information layers
    #
    # Takes stratigraphy graph created by map2model c++ code to generate list of groups found in the region of interest
    # Uses first of each possible set of toplogies per unit and per group, which is arbitrary. On the other hand we
    # are not checking relative ages again to see if this helps reduce ambiguity, which I think it would.
    ####################################
    @beartype.beartype
    def save_group(self, config: Config, map_data, stratColumn: StratigraphicColumn):
        glabels = self.get_group_labels()

        geol = map_data.get_map_data(Datatype.GEOLOGY)
        ordered_groups = self.abs_age_groups(geol, config.tmp_path)
        ordered_groups.set_index("group_", inplace=True)
        local_geol = geol.copy()
        local_geol.drop_duplicates(subset="UNIT_NAME", inplace=True)
        local_geol.set_index("UNIT_NAME", inplace=True)

        if config.verbose_level == VerboseLevel.ALL:
            display(ordered_groups)

        # Make a sub graph copy of only the group nodes
        Gp = self.graph.subgraph(
            [x for x, y in self.graph.nodes(data=True) if "isGroup" in y]
        ).copy()

        # Add group edges if any formation within a group has an edge with a corresponding
        # edge in another group the direction is based on the average age of each group
        for e in self.graph.edges:
            if self.graph.nodes[e[0]]["gid"] != self.graph.nodes[e[1]]["gid"]:
                glabel_0 = self.graph.nodes[e[0]]["LabelGraphics"]["text"]
                glabel_1 = self.graph.nodes[e[1]]["LabelGraphics"]["text"]

                if str(local_geol.loc[glabel_0]["GROUP"]) == "None":
                    grp0 = glabel_0
                else:
                    grp0 = local_geol.loc[glabel_0]["GROUP"]
                if str(local_geol.loc[glabel_1]["GROUP"]) == "None":
                    grp1 = glabel_1
                else:
                    grp1 = local_geol.loc[glabel_1]["GROUP"]

                if grp0 in glabels.values() and grp1 in glabels.values():
                    if (
                        ordered_groups.loc[grp0]["ave"]
                        < ordered_groups.loc[grp1]["ave"]
                    ):
                        Gp.add_edge(
                            self.graph.nodes[e[0]]["gid"], self.graph.nodes[e[1]]["gid"]
                        )

        # remove duplicate edges with opposite directions
        GpBuffer = Gp.copy()
        for e in GpBuffer.edges:
            for f in GpBuffer.edges:
                # arbitrary choice to ensure edge is not completely removed
                if e[0] == f[1] and e[1] == f[0] and e[0] < f[0]:
                    Gp.remove_edge(e[0], e[1])

        if config.verbose_level == VerboseLevel.ALL:
            plt.figure(1)  # display strat graph for one group
            plt.title("Groups")
            if len(glabels) > 1:
                nx.draw_networkx(
                    Gp, pos=nx.kamada_kawai_layout(Gp), arrows=True, with_labels=False
                )
                labels = {x: Gp.nodes[x]["LabelGraphics"]["text"] for x in Gp.nodes()}
                nx.draw_networkx_labels(
                    Gp,
                    pos=nx.kamada_kawai_layout(Gp),
                    labels=labels,
                    font_size=12,
                    font_family="sans-serif",
                )
            plt.show()

        # As permutations of groups can be very large if groups (with an unknown
        # order) > 10 we take the first possibility
        if len(ordered_groups) > 10:
            # first possible sorted directional graphs
            glist = [list(nx.topological_sort(Gp))]
            if config.verbose_level != VerboseLevel.NONE:
                print("group choices: 1 (more than 10 groups)")
        else:
            # all possible sorted directional graphs
            glist = list(nx.all_topological_sorts(Gp))
            if config.verbose_level != VerboseLevel.NONE:
                print("group choices:", len(glist))

        # glist is a list of lists
        f = open(os.path.join(config.tmp_path, "groups.csv"), "w")
        glen = len(glist)
        if glen > 100:
            glen = 100
        for n in range(0, glen):
            f.write("Choice " + str(n))
            for m in range(0, len(glist[0])):
                f.write(
                    "," + str(self.graph.nodes[glist[n][m]]["LabelGraphics"]["text"])
                )
            f.write("\n")
        f.close()

        nx.write_gml(Gp, os.path.join(config.tmp_path, "groups.gml"))

        contents = np.genfromtxt(
            os.path.join(config.tmp_path, "groups.csv"), delimiter=",", dtype="U100"
        )
        self.contents = contents
        self.ordered_groups = ordered_groups

        # Set Group order as an index and add the number of formations in each group
        orderedGroupNames = ordered_groups.index[:]
        groupIndexes = dict(
            zip(orderedGroupNames, range(1, len(orderedGroupNames) + 1))
        )
        groupNumElements = dict(
            zip(
                orderedGroupNames,
                [
                    len(
                        stratColumn.stratigraphicUnits[
                            stratColumn.stratigraphicUnits.group == x
                        ]
                    )
                    for x in orderedGroupNames
                ],
            )
        )
        stratColumn.stratigraphicUnits["groupNum"] = stratColumn.stratigraphicUnits[
            "group"
        ].map(groupIndexes)
        stratColumn.stratigraphicUnits["numInGroup"] = stratColumn.stratigraphicUnits[
            "group"
        ].map(groupNumElements)

        # Reorder strat column by group and formation ordering
        stratColumn.stratigraphicUnits.sort_values(
            by=["groupNum", "indexInGroup"], inplace=True
        )

        # Output to all_sorts.csv file
        df = stratColumn.stratigraphicUnits[
            ["groupNum", "indexInGroup", "numInGroup", "name", "group"]
        ].copy()
        df.rename(
            columns={
                "groupNum": "group number",
                "indexInGroup": "index in group",
                "numInGroup": "number in group",
                "name": "code",
            },
            inplace=True,
        )
        df.reset_index(drop=True, inplace=True)
        df.index.name = "index"
        df.to_csv(os.path.join(config.tmp_path, "all_sorts.csv"))

    ####################################
    # save out fault fault relationship information as array
    #
    # parse_fault_relationships(graph_path,tmp_path,output_path)
    # graph_path path to graphs
    # tmp_path path to tmp directory
    # output_path path to output directory
    #
    # Saves fault vs unit, group and fault relationship tables using outputs from map2model c++ code
    ####################################
    @beartype.beartype
    def parse_fault_relationships(
        self, config: Config, mapData, stratColumn: StratigraphicColumn
    ):
        # Get list of fault ids in fault dimensions as not all faults needed
        df = pd.read_csv(
            os.path.join(config.output_path, "fault_dimensions.csv"), delimiter=","
        )
        faultInfo = df[["Fault"]].copy()
        faultInfo["FaultId"] = faultInfo["Fault"].replace("Fault_", "")
        faultInfo = faultInfo.reset_index()

        # Create int dataframe of unit and fault intersections
        # Parse unit fault intersections from map2model output
        df = pd.read_csv(
            os.path.join(config.graph_path, "unit-fault-intersection.txt"),
            delimiter="{",
            header=None,
        )
        df[0] = list(df[0].str.replace("^[0-9]*, ", "", regex=True))
        df[0] = list(df[0].str.replace(", ", ""))
        df[1] = list(df[1].str.replace("}", "", regex=False))
        df.loc[:, 1] = df.loc[:, 1].str.split(",")
        df.rename(columns={0: "code", 1: "FaultList"}, inplace=True)

        # Populate Fault columns with intersections
        for faultName in faultInfo["Fault"]:
            df[faultName] = 0
        for i in range(len(df)):
            for j in range(len(faultInfo)):
                if faultInfo["FaultId"][j] in df.loc[i, "FaultList"]:
                    df.loc[i, faultInfo["Fault"][j]] = 1
        df.drop(columns="FaultList", inplace=True)
        unitFaultIntersections = df.set_index("code")
        unitFaultIntersections.to_csv(
            os.path.join(config.output_path, "unit-fault-relationships.csv")
        )

        # Get group and supergroup summary
        # TODO: Get these from strat column
        summary = pd.read_csv(os.path.join(config.tmp_path, "all_sorts_clean.csv"))

        # Output group fault relationships table
        df = unitFaultIntersections.join(
            summary[["code", "group"]].set_index("code"), on="code"
        )
        df.set_index("group", inplace=True)
        groupFaultIntersections = df.groupby("group", sort=False).max()
        groupFaultIntersections.index.name = "group"
        groupFaultIntersections.to_csv(
            os.path.join(config.output_path, "group-fault-relationships.csv")
        )

        # Output supergroup fault relationships table
        df = unitFaultIntersections.join(
            summary[["code", "supergroup"]].set_index("code"), on="code"
        )
        df.set_index("supergroup", inplace=True)
        supergroupFaultIntersections = df.groupby("supergroup", sort=False).max()
        supergroupFaultIntersections.index.name = "supergroup"
        supergroupFaultIntersections.to_csv(
            os.path.join(config.output_path, "supergroup-fault-relationships.csv")
        )

        # Parse fault fault intersections from map2model output
        # Note: The faults in this file do not always match the unit fault intersection faults!!!
        df = pd.read_csv(
            os.path.join(config.graph_path, "fault-fault-intersection.txt"),
            delimiter="{",
            header=None,
        )
        df[0] = list(df[0].str.replace("^[0-9]*, ", "", regex=True))
        df[0] = list(df[0].str.replace(",", "", regex=False))
        df[0] = list(df[0].str.replace(" ", "", regex=False))
        df[0] = "Fault_" + df[0]
        df[1] = list(df[1].str.replace("}", "", regex=False))
        df[1] = [re.findall("\(.*?\)", i) for i in df[1]]
        for j in range(len(df)):
            df[1][j] = [i.strip("()").replace(" ", "").split(",") for i in df[1][j]]
        df.set_index(0, inplace=True)

        # Export fault network edges (include to fault nodes that have been culled)
        faultIntersectionList = []
        for fault1 in df.index:
            for j in df[1][fault1]:
                faultIntersectionList.append(
                    [str(fault1), "Fault_" + str(j[0]), j[2], j[1]]
                )
        faultIntersectionEdges = pd.DataFrame(
            data=faultIntersectionList, columns=["fault_1", "fault_2", "angle", "topol"]
        )
        faultIntersectionEdges = faultIntersectionEdges.set_index("fault_1")
        faultIntersectionEdges.to_csv(
            os.path.join(config.tmp_path, "fault_network_edges.csv")
        )

        # Create fault to fault array of only restricted fault list
        df2 = faultInfo[["Fault"]].copy()
        df2.rename(columns={"Fault": "fault_id"}, inplace=True)
        df2.set_index("fault_id", inplace=True)
        for faultName in faultInfo["Fault"]:
            df2[faultName] = 0
        for idx, row in faultIntersectionEdges.iterrows():
            if idx in df2.index and row["fault_2"] in df2.columns:
                df2[row["fault_2"]][idx] = 1
        faultFaultIntersections = df2.copy()
        faultFaultIntersections.to_csv(
            os.path.join(config.output_path, "fault-fault-relationships.csv")
        )

        # Create directional graph network of fault network
        G = nx.DiGraph()
        for fault in faultInfo["Fault"]:
            G.add_node(fault)
        for fault1 in faultFaultIntersections.index:
            for fault2 in faultInfo["Fault"]:
                if faultFaultIntersections[fault2][fault1]:
                    G.add_edge(fault1, fault2)

        # Display graph of fault network
        if (
            config.verbose_level == VerboseLevel.ALL
            and len(faultFaultIntersections) > 0
        ):
            fig, ax = plt.subplots()
            nx.draw(G, ax=ax, with_labels=True, font_weight="bold")
            plt.title("Fault Network")
            plt.show()

        # check for and remove first edge from fault cycle
        cycles = list(nx.simple_cycles(G))
        for i in range(0, len(cycles)):
            if G.has_edge(cycles[i][0], cycles[i][1]):
                G.remove_edge(cycles[i][0], cycles[i][1])
                print("fault cycle removed:", cycles[i][0], cycles[i][1])

        # add angle and type of fault intersection to graph
        df = faultIntersectionEdges
        for i in range(len(df)):
            if (
                G.has_node(df.index[i])
                and G.has_node(df["fault_2"][i])
                and G.has_edge(df.index[i], df["fault_2"][i])
            ):
                G[df.index[i]][df["fault_2"][i]]["angle"] = df["angle"][i]
                G[df.index[i]][df["fault_2"][i]]["topol"] = df["topol"][i]

        # export graph of fault network to tmp path
        nx.write_gml(G, os.path.join(config.tmp_path, "fault_network.gml"))

    @beartype.beartype
    def super_groups_and_groups(group_girdle, config: Config, map_data, workflow: dict):
        group_girdle = pd.DataFrame.from_dict(group_girdle, orient="index")
        group_girdle.columns = ["plunge", "bearing", "num orientations"]
        group_girdle.sort_values(by="num orientations", ascending=False, inplace=True)

        l, m, n = m2l_utils.ddd2dircos(
            group_girdle.iloc[0]["plunge"], group_girdle.iloc[0]["bearing"]
        )
        super_group = pd.DataFrame(
            [[group_girdle[0:1].index[0], "Super_Group_0", l, m, n]],
            columns=["Group", "Super_Group", "l", "m", "n"],
        )
        super_group.set_index("Group", inplace=True)

        # geol = gpd.read_file(os.path.join(config.tmp_path, 'geol_clip.shp'))
        local_geol = map_data.get_map_data(Datatype.GEOLOGY).copy()
        local_geol.drop_duplicates(subset="UNIT_NAME", keep="first", inplace=True)
        local_geol = local_geol.set_index("GROUP")

        sg_index = 0
        for i in range(0, len(group_girdle)):
            # if(config.c_l['intrusive'] in geol.loc[group_girdle.iloc[i].name.replace("_"," ")]["ROCKTYPE1"]
            # and config.c_l['sill'] not in geol.loc[group_girdle.iloc[i].name.replace("_"," ")]["DESCRIPTION"]):
            if group_girdle.iloc[i].name in local_geol["UNIT_NAME"]:
                if (
                    config.c_l["intrusive"]
                    in local_geol.loc[group_girdle.iloc[i].name]["ROCKTYPE1"]
                    and config.c_l["sill"]
                    not in local_geol.loc[group_girdle.iloc[i].name]["DESCRIPTION"]
                ):
                    sg_index = sg_index + 1
                    sgname = "Super_Group_" + str(sg_index)
                    super_group_new = pd.DataFrame(
                        [[group_girdle[i : i + 1].index[0], sgname, l, m, n]],
                        columns=["Group", "Super_Group", "l", "m", "n"],
                    )
                    super_group_new.set_index("Group", inplace=True)
                    super_group = pd.concat([super_group, super_group_new])

                elif group_girdle.iloc[i]["num orientations"] > 5:
                    l, m, n = m2l_utils.ddd2dircos(
                        group_girdle.iloc[i]["plunge"], group_girdle.iloc[i]["bearing"]
                    )

                    found = False
                    sg_i = 0
                    for ind, sg in super_group.iterrows():

                        c = sg["l"] * l + sg["m"] * m + sg["n"] * n
                        if c > 1:
                            c = 1
                        c = degrees(acos(c))
                        if c < config.run_flags["misorientation"] and not found:
                            found = True
                            sgname = "Super_Group_" + str(sg_i)
                            super_group_old = pd.DataFrame(
                                [[group_girdle[i : i + 1].index[0], sgname, l, m, n]],
                                columns=["Group", "Super_Group", "l", "m", "n"],
                            )
                            super_group_old.set_index("Group", inplace=True)
                            super_group = pd.concat([super_group, super_group_old])
                        sg_i = sg_i + 1
                    if not found:

                        sg_index = sg_index + 1
                        sgname = "Super_Group_" + str(sg_index)
                        super_group_new = pd.DataFrame(
                            [[group_girdle[i : i + 1].index[0], sgname, l, m, n]],
                            columns=["Group", "Super_Group", "l", "m", "n"],
                        )
                        super_group_new.set_index("Group", inplace=True)
                        super_group = pd.concat([super_group, super_group_new])
                else:
                    # not enough orientations to test, so lumped with group with most orientations
                    sgname = "Super_Group_" + str(0)
                    super_group_old = pd.DataFrame(
                        [[group_girdle[i : i + 1].index[0], sgname, l, m, n]],
                        columns=["Group", "Super_Group", "l", "m", "n"],
                    )
                    super_group_old.set_index("Group", inplace=True)
                    super_group = pd.concat([super_group, super_group_old])

        use_gcode3 = []
        for ind, sg in super_group.iterrows():
            clean = ind.replace(" ", "_").replace("-", "_")
            use_gcode3.append(clean)
        if workflow["cover_map"]:
            use_gcode3.append("cover")

        sg2 = set(super_group["Super_Group"])
        super_groups = []
        if workflow["cover_map"]:
            super_groups.append(["cover"])
        for s in sg2:
            temp = []
            for ind, sg in super_group.iterrows():
                if s == sg["Super_Group"]:
                    temp.append(ind)
            super_groups.append(temp)

        f = open(os.path.join(config.tmp_path, "super_groups.csv"), "w")
        for sg in super_groups:
            for s in sg:
                f.write(str(s) + ",")
            f.write("\n")
        f.close()
        return (super_groups, use_gcode3)

    @beartype.beartype
    def use_asud(self, config: Config):
        if config.verbose_level != VerboseLevel.NONE:
            print("Resolving ambiguities using ASUD...", end="\toutput_dir:")
        if type(self.asud) == None:
            print("Failed to load ASUD data, continuing without ASUD help")
            return

        G = self.graph
        Gp = G.copy().to_directed()
        for e in G.edges:
            glabel_0 = G.nodes[e[0]]["LabelGraphics"]["text"]
            glabel_1 = G.nodes[e[1]]["LabelGraphics"]["text"]

            edge = self.asud.loc[
                (self.asud["over"] == glabel_0) & (self.asud["under"] == glabel_1)
            ]
            if len(edge) > 0 and (
                not ("isGroup" in Gp.nodes[e[0]]) and not ("isGroup" in Gp.nodes[e[1]])
            ):
                if Gp.has_edge(e[1], e[0]):
                    Gp.remove_edge(e[1], e[0])
                if Gp.has_edge(e[0], e[1]):
                    Gp.remove_edge(e[0], e[1])
                Gp.add_edge(e[0], e[1])

            edge = self.asud.loc[
                (self.asud["under"] == glabel_0) & (self.asud["over"] == glabel_1)
            ]
            if len(edge) > 0 and (
                not ("isGroup" in Gp.nodes[e[0]]) and not ("isGroup" in Gp.nodes[e[1]])
            ):
                if Gp.has_edge(e[0], e[1]):
                    Gp.remove_edge(e[0], e[1])
                if Gp.has_edge(e[1], e[0]):
                    Gp.remove_edge(e[1], e[0])
                Gp.add_edge(e[1], e[0])

        # recode = {}
        # i = 0
        # for n in Gp.nodes:
        # if "isGroup" in Gp.nodes[n]:
        # recode[n] = i
        # i = i + 1

        # for n in Gp.nodes:
        # if "isGroup" not in Gp.nodes[n]:
        # Gp.nodes[n]["gid"] = recode[Gp.nodes[n]["gid"]]

        # check for cycle and arbitrarily remove first of the edges
        try:
            cycles = list(nx.simple_cycles(G))
            for c in cycles:
                G.remove_edge(c[0], c[1])
                warning_msg = (
                    'map2loop warning: The stratigraphic relationship: "'
                    + c[0]
                    + " overlies "
                    + c[1]
                    + '" was removed as it conflicts with another relationship'
                )
                warnings.warn(warning_msg)
        except Exception:
            print("no cycles")

        nx.write_gml(Gp, os.path.join(self.graph_path, "ASUD_strat.gml"))
        self.graph = Gp
        if config.verbose_level != VerboseLevel.NONE:
            print("Done.")

    ####################################
    # combine multiple outputs into single graph that contains all infor needed by LoopStructural
    #
    # make_Loop_graph(tmp_path,output_path)
    # tmp_path path to tmp directory
    # output_path path to output directory
    #
    # Returns networkx graph
    ####################################

    @beartype.beartype
    def make_Loop_graph(config: Config, map_data, point_data):
        Gloop = nx.DiGraph()
        strats = map_data.get_map_data(Datatype.GEOLOGY).copy()
        strats.drop_duplicates(subset=["UNIT_NAME"], inplace=True)
        strats["UNIT_NAME"] = strats["UNIT_NAME"].replace(" ", "_").replace("-", "_")
        strats.set_index(["UNIT_NAME"], inplace=True)

        # Load faults and stratigraphy
        Gf = nx.read_gml(os.path.join(config.tmp_path, "fault_network.gml"))
        Astrat = pd.read_csv(
            os.path.join(config.tmp_path, "all_sorts_clean.csv"), sep=","
        )

        # add formation stratigraphy to graph as nodes
        Astrat = Astrat.set_index("code")
        for ind, s in Astrat.iterrows():
            if s.name != "cover" and s.name != "cover_up":
                Gloop.add_node(
                    s.name,
                    s_colour=s["colour"],
                    ntype="formation",
                    group=s["group"],
                    StratType=s["strat_type"],
                    uctype=s["uctype"],
                    GroupNumber=s["group number"],
                    IndexInGroup=s["index in group"],
                    NumberInGroup=s["number in group"],
                    MinAge=strats.loc[s.name]["MIN_AGE"],
                    MaxAge=strats.loc[s.name]["MAX_AGE"],
                )
            else:
                Gloop.add_node(
                    s.name,
                    s_colour=s["colour"],
                    ntype="formation",
                    group=s["group"],
                    StratType=s["strat_type"],
                    uctype=s["uctype"],
                    GroupNumber=s["group number"],
                    IndexInGroup=s["index in group"],
                    NumberInGroup=s["number in group"],
                    MinAge=0,
                    MaxAge=1,
                )

        # add formation-formation stratigraphy to graph as edges
        i = 0
        for ind, s in Astrat.iterrows():
            if ind != Astrat.index[-1]:
                Gloop.add_edge(Astrat.iloc[i].name, Astrat.iloc[i + 1].name)
                Gloop[Astrat.iloc[i].name][Astrat.iloc[i + 1].name][
                    "etype"
                ] = "formation_formation"
            i = i + 1

        # add faults to graph as nodes
        Af_d = pd.read_csv(
            os.path.join(config.output_path, "fault_dimensions.csv"), sep=","
        )
        for ind, f in Af_d.iterrows():
            Gloop.add_node(f["Fault"], ntype="fault")

        # add fault centroid to node
        Afgeom = pd.read_csv(os.path.join(config.output_path, "faults.csv"), sep=",")

        pd.to_numeric(Afgeom["X"], downcast="float")
        pd.to_numeric(Afgeom["Y"], downcast="float")
        pd.to_numeric(Afgeom["Z"], downcast="float")

        for n in Gloop.nodes:
            if "Fault" in n:
                subset = Afgeom[Afgeom["formation"] == n]
                Gloop.nodes[n]["Xmean"] = subset["X"].mean()
                Gloop.nodes[n]["Ymean"] = subset["Y"].mean()
                Gloop.nodes[n]["Zmean"] = subset["Z"].mean()

        for e in Gf.edges:
            if Gloop.has_node(e[0]) and Gloop.has_node(e[1]):
                Gloop.add_edge(e[0], e[1])
                Gloop[e[0]][e[1]]["Angle"] = Gf.edges[e]["angle"]
                Gloop[e[0]][e[1]]["Topol"] = Gf.edges[e]["topol"]
                Gloop[e[0]][e[1]]["etype"] = "fault_fault"

        # add fault dimension info to fault nodes
        Af_d = pd.read_csv(
            os.path.join(config.output_path, "fault_dimensions.csv"), sep=","
        )
        fault_l = pd.DataFrame(Af_d["Fault"])
        Af_d = Af_d.set_index("Fault")
        fault_l = fault_l.set_index("Fault")
        fault_l["incLength"] = Af_d["incLength"]

        for n in Gloop.nodes:
            if "Fault" in n and n in Af_d.index:
                if Af_d.loc[n].name in Gloop.nodes:
                    Gloop.nodes[n]["HorizontalRadius"] = Af_d.loc[n]["HorizontalRadius"]
                    Gloop.nodes[n]["VerticalRadius"] = Af_d.loc[n]["VerticalRadius"]
                    Gloop.nodes[n]["InfluenceDistance"] = Af_d.loc[n][
                        "InfluenceDistance"
                    ]
                    Gloop.nodes[n]["IncLength"] = Af_d.loc[n]["incLength"]
                    Gloop.nodes[n]["f_colour"] = Af_d.loc[n]["colour"]

        # add fault orientation info and clustering on orientation and length to fault nodes
        Af_o = pd.read_csv(
            os.path.join(config.output_path, "fault_orientations.csv"), sep=","
        )
        Af_o = Af_o.drop_duplicates(subset="formation")
        Af_o = Af_o.set_index("formation")

        pd.to_numeric(Af_o["dip"], downcast="float")
        pd.to_numeric(Af_o["DipDirection"], downcast="float")
        from sklearn import cluster
        from sklearn.exceptions import ConvergenceWarning

        def convert_dd_lmn(dip, dipdir):
            r_dd = np.radians(dipdir)
            r_90_dip = np.radians(90 - dip)
            sin_dd = np.sin(r_dd)
            cos_dd = np.cos(r_dd)

            cos_90_dip = np.cos(r_90_dip)
            sin_90_dip = np.sin(r_90_dip)

            l = sin_dd * cos_90_dip
            m = cos_dd * cos_90_dip
            n = sin_90_dip
            return (l, m, n)

        l, m, n = convert_dd_lmn(Af_o["dip"], Af_o["DipDirection"])
        faults = pd.DataFrame(l, columns=["l"])
        faults["m"] = m
        faults["n"] = n

        # with warnings.catch_warnings():
        #     warnings.simplefilter(action="ignore", category=ConvergenceWarning)
        #     # warnings.simplefilter("ignore")
        #     if len(Af_d) >= config.run_flags["fault_orientation_clusters"]:
        #         k_means_o = cluster.KMeans(
        #             n_clusters=config.run_flags["fault_orientation_clusters"],
        #             max_iter=50,
        #             random_state=1,
        #         )
        #         k_means_o.fit(faults)
        #         labels_o = k_means_o.labels_
        #         clusters_o = pd.DataFrame(labels_o, columns=["cluster_o"])

        #     if len(Af_d) >= config.run_flags["fault_length_clusters"]:
        #         k_means_l = cluster.KMeans(
        #             n_clusters=config.run_flags["fault_length_clusters"],
        #             max_iter=50,
        #             random_state=1,
        #         )
        #         k_means_l.fit(fault_l)
        #         labels_l = k_means_l.labels_
        #         clusters_l = pd.DataFrame(labels_l, columns=["cluster_l"])

        Af_o = Af_o.reset_index()
        # if len(Af_d) >= config.run_flags["fault_orientation_clusters"]:
        #     Af_o["cluster_o"] = clusters_o["cluster_o"]
        # if len(Af_d) >= config.run_flags["fault_length_clusters"]:
        #     Af_o["cluster_l"] = clusters_l["cluster_l"]
        Af_o = Af_o.set_index(keys="formation")
        Af_o.to_csv(os.path.join(config.output_path, "fault_clusters.csv"))
        for n in Gloop.nodes:
            if "Fault" in n and n in Af_o:
                if Af_o.loc[n].name in Gloop.nodes:
                    Gloop.nodes[n]["Dip"] = Af_o.loc[n]["dip"]
                    Gloop.nodes[n]["DipDirection"] = Af_o.loc[n]["DipDirection"]
                    Gloop.nodes[n]["DipPolarity"] = Af_o.loc[n]["DipPolarity"]
                    # if len(Af_d) >= config.run_flags["fault_orientation_clusters"]:
                    #     Gloop.nodes[n]["OrientationCluster"] = Af_o.loc[n]["cluster_o"]
                    # else:
                    #     Gloop.nodes[n]["OrientationCluster"] = -1
                    # if len(Af_d) >= config.run_flags["fault_length_clusters"]:
                    #     Gloop.nodes[n]["LengthCluster"] = Af_o.loc[n]["cluster_l"]
                    # else:
                    #     Gloop.nodes[n]["LengthCluster"] = -1

        # add centrality measures to fault nodes

        Gfcc = nx.closeness_centrality(Gf)
        Gfbc = nx.betweenness_centrality(Gf)

        for n in Gloop.nodes:
            if "Fault" in n:
                if n in Gfcc.keys():
                    Gloop.nodes[n]["ClosenessCentrality"] = Gfcc[n]
                    Gloop.nodes[n]["BetweennessCentrality"] = Gfbc[n]
                else:
                    Gloop.nodes[n]["ClosenessCentrality"] = -1
                    Gloop.nodes[n]["BetweennessCentrality"] = -1

        # add formation thickness info to formation nodes
        As_t = pd.read_csv(
            os.path.join(config.output_path, "formation_summary_thicknesses.csv"),
            sep=",",
        )
        As_t = As_t.set_index("formation")

        for n in Gloop.nodes:
            if "Fault" not in n and "Point_data" not in n:
                if As_t.loc[n].name in Gloop.nodes:
                    if np.isnan(As_t.loc[n]["thickness std"]):
                        Gloop.nodes[n]["ThicknessStd"] = -1
                    else:
                        Gloop.nodes[n]["ThicknessStd"] = As_t.loc[n]["thickness std"]
                    if np.isnan(As_t.loc[n]["thickness median"]):
                        Gloop.nodes[n]["ThicknessMedian"] = -1
                    else:
                        Gloop.nodes[n]["ThicknessMedian"] = As_t.loc[n][
                            "thickness median"
                        ]

                    Gloop.nodes[n]["ThicknessMethod"] = As_t.loc[n]["method"]

        # add group-formation stratigraphy to graph as edges
        for ind, s in Astrat.iterrows():
            if not s["group"] + "_gp" in Gloop.nodes():
                Gloop.add_node(s["group"] + "_gp", ntype="group")
            Gloop.add_edge(s["group"] + "_gp", s.name)
            Gloop[s["group"] + "_gp"][s.name]["etype"] = "group_formation"

        # add group-fault edges to graph
        Af_s = pd.read_csv(
            os.path.join(config.output_path, "group-fault-relationships.csv"), sep=","
        )
        for ind, s in Af_s.iterrows():
            if s["group"] + "_gp" in Gloop.nodes:
                for col in Af_s.columns[1:]:
                    if col in Gloop.nodes:
                        if Af_s.loc[ind, col] == 1:
                            Gloop.add_edge(col, s["group"] + "_gp")
                            Gloop[col][s["group"] + "_gp"]["etype"] = "fault_group"

        # add group-group edges to graph
        Ag_g = pd.read_csv(
            os.path.join(config.tmp_path, "groups_clean.csv"),
            header=None,
            index_col=None,
        )
        i = 0
        for ind, g in Ag_g.iterrows():
            if ind != Ag_g.index[-1]:
                if (
                    Ag_g.iloc[i][0] + "_gp" in Gloop.nodes
                    and Ag_g.iloc[i + 1][0] + "_gp" in Gloop.nodes
                ):
                    Gloop.add_edge(Ag_g.iloc[i][0] + "_gp", Ag_g.iloc[i + 1][0] + "_gp")
                    Gloop[Ag_g.iloc[i][0] + "_gp"][Ag_g.iloc[i + 1][0] + "_gp"][
                        "etype"
                    ] = "group_group"
            i = i + 1

        # add supergroups as nodes and supergroup-group relationships as edges
        sgi = 0
        supergroups = {}
        with open(os.path.join(config.tmp_path, "super_groups.csv")) as sgf:
            lines = sgf.readlines()
            for l in lines:
                Gloop.add_node("supergroup_{}".format(sgi), ntype="supergroup")
                for g in l.split(","):
                    g = g.replace("-", "_").replace(" ", "_").rstrip()
                    if g:
                        supergroups[g] = "supergroup_{}".format(sgi)
                        if g + "_gp" in Gloop.nodes():
                            Gloop.add_edge(supergroups[g], g + "_gp")
                            Gloop[supergroups[g]][g + "_gp"][
                                "etype"
                            ] = "supergroup_group"
                sgi += 1

        # add all geolocated data to a single unconnected node
        Gloop.add_node("Point_data", ntype="points", data=point_data)

        with map_data.get_map_data(Datatype.DTM).open() as dtm:
            dtm_data = dtm.read(1)
            bounds = dtm.bounds
            minx = bounds.left
            miny = bounds.bottom
            maxx = bounds.right
            maxy = bounds.top
            xscale = (maxx - minx) / dtm_data.shape[1]
            yscale = (maxy - miny) / dtm_data.shape[0]
            Gloop.add_node(
                "DTM_data",
                ntype="dtm",
                data=str(dtm_data.tolist()),
                shape=str(dtm_data.shape),
                minx=minx,
                miny=miny,
                maxx=maxx,
                maxy=maxy,
                xscale=xscale,
                yscale=yscale,
            )

        Gloop.add_node("bbox", ntype="bbox", data=str(config.bbox))

        Gloop.add_node(
            "dst_crs", ntype="dst_crs", data=str(map_data.working_projection)
        )

        config_metadata = {
            "project_path": config.project_path,
            "geology_file": map_data.get_filename(Datatype.GEOLOGY),
            "fault_file": map_data.get_filename(Datatype.FAULT),
            "fold_file": map_data.get_filename(Datatype.FOLD),
            "structure_file": map_data.get_filename(Datatype.STRUCTURE),
            "mindep_file": map_data.get_filename(Datatype.MINERAL_DEPOSIT),
            "section_files": map_data.get_filename(Datatype.SECTION),
            "drillhole_file": map_data.get_filename(Datatype.DRILLHOLE),
            "dtb_grid_file": map_data.get_filename(Datatype.DTB_GRID),
            "cover_map_file": map_data.get_filename(Datatype.COVER_MAP),
            "bbox_3d": config.bbox_3d,
            "step_out": config.step_out,
            "dtm_crs": config.dtm_crs,
            "project_crs": map_data.working_projection,
        }
        Gloop.add_node(
            "metadata",
            ntype="metadata",
            run_flags=str(config.run_flags),
            c_l=str(config.c_l),
            config=str(config_metadata),
        )
        return Gloop

    def colour_Loop_graph(output_path, prefix):
        with open(os.path.join(output_path, prefix + ".gml")) as graph:
            new_graph = open(os.path.join(output_path, prefix + "_colour.gml"), "w")
            lines = graph.readlines()
            for l in lines:
                if "s_colour" in l:
                    new_graph.write(
                        "    "
                        + l.strip().replace("s_colour", "graphics [ fill ")
                        + " ]\n"
                    )
                elif "f_colour" in l:
                    new_graph.write(
                        "    "
                        + l.strip().replace(
                            "f_colour", 'graphics [ type "ellipse" fill '
                        )
                        + " ]\n"
                    )
                elif 'ntype "supergroup"' in l:
                    new_graph.write('    graphics [ type "triangle" fill "#FF0000" ]\n')
                    new_graph.write(l)
                elif 'ntype "group"' in l:
                    new_graph.write('    graphics [ type "triangle" fill "#FF9900" ]\n')
                    new_graph.write(l)
                elif (
                    'ntype "points"' in l
                    or 'ntype "metadata"' in l
                    or 'ntype "dtm"' in l
                    or 'ntype "bbox"' in l
                    or 'ntype "dst_crs"' in l
                ):
                    new_graph.write('    graphics [ type "octagon" fill "#00FF00" ]\n')
                    new_graph.write(l)
                elif 'etype "formation_formation"' in l:
                    new_graph.write(
                        '    graphics [ style "line" arrow "last" fill "#6666FF" ]\n'
                    )
                    new_graph.write(l)
                elif 'etype "fault_group"' in l:
                    new_graph.write(
                        '    graphics [ style "line" arrow "last" fill "#666600" ]\n'
                    )
                    new_graph.write(l)
                elif 'etype "fault_formation"' in l:
                    new_graph.write(
                        '    graphics [ style "line" arrow "last" fill "#0066FF" ]\n'
                    )
                    new_graph.write(l)
                elif 'etype "fault_fault"' in l:
                    new_graph.write(
                        '    graphics [ style "line" arrow "last" fill "#000000" ]\n'
                    )
                    new_graph.write(l)
                elif 'etype "group_formation"' in l:
                    new_graph.write(
                        '    graphics [ style "line" arrow "last" fill "#FF6600" ]\n'
                    )
                    new_graph.write(l)
                else:
                    new_graph.write(l)
            new_graph.close()
