# from map2loop import m2l_topology
import networkx as nx
import random
import numpy as np
import pandas as pd
import os
import sys
import geopandas as gpd
import rasterio
from rasterio import plot
from rasterio.plot import show
from rasterio.mask import mask
from rasterio.transform import from_origin
from rasterio.io import MemoryFile
import matplotlib
from . import m2l_utils
import LoopProjectFile
import beartype
from .config import Config
import re


##########################################################################
# Save out and compile taskfile needed to generate geomodeller model using the geomodellerbatch engine
#
# loop2geomodeller(test_data_path,tmp_path,output_path,save_faults,compute_etc)
# Args:
# test_data_path root directory of test data
# tmp_path directory of temporary outputs
# output_path directory of outputs
# ave_faults flag for saving faults or not
# compute_etc flag for actual calculations or just project output
#
# Creates geomodeller taskfile files from varous map2loop outputs
##########################################################################
def loop2geomodeller(
    model_name,
    test_data_path,
    tmp_path,
    output_path,
    dtm_file,
    bbox,
    model_top,
    model_base,
    save_faults,
    compute_etc,
    workflow,
):

    f = open(os.path.join(test_data_path, model_name, "m2l.taskfile"), "w")
    f.write("#---------------------------------------------------------------\n")
    f.write("#-----------------------Project Header-----------------------\n")
    f.write("#---------------------------------------------------------------\n")
    f.write('name: "UWA_Intrepid"\n')
    f.write('description: "Automate_batch_Model"\n')
    f.write("    GeomodellerTask {\n")
    f.write("    CreateProject {\n")
    f.write('        name: "Hamersley"\n')
    f.write('        author: "Mark"\n')
    f.write('        date: "23/10/2019  0: 0: 0"\n')
    f.write('        projection { map_projection: "GDA94 / MGA50"}\n')
    f.write('        version: "2.0"\n')
    f.write("        units: meters\n")
    f.write("        precision: 1.0\n")
    f.write("        Extents {\n")
    f.write("            xmin: " + str(bbox[0]) + "\n")
    f.write("            ymin: " + str(bbox[1]) + "\n")
    f.write("            zmin: " + str(model_base) + "\n")
    f.write("            xmax: " + str(bbox[2]) + "\n")
    f.write("            ymax: " + str(bbox[3]) + "\n")
    f.write("            zmax: " + str(model_top) + "\n")
    f.write("        }\n")
    f.write("        deflection2d: 0.001\n")
    f.write("        deflection3d: 0.001\n")
    f.write("        discretisation: 10.0\n")
    f.write("        referenceTop: false\n")
    f.write("        CustomDTM {\n")
    f.write("            Extents {\n")
    f.write("            xmin: " + str(bbox[0]) + "\n")
    f.write("            ymin: " + str(bbox[1]) + "\n")
    f.write("            xmax: " + str(bbox[2]) + "\n")
    f.write("            ymax: " + str(bbox[3]) + "\n")
    f.write("            }\n")
    f.write('            name: "Topography"\n')
    f.write("            filename {\n")
    f.write('                Grid_Name: "' + dtm_file + '"\n')
    f.write("            }\n")
    f.write("            nx: 10\n")
    f.write("            ny: 10\n")
    f.write("        }\n")
    f.write("    }\n")
    f.write("}\n")

    orientations = pd.read_csv(
        os.path.join(output_path, "orientations_clean.csv"), sep=","
    )
    contacts = pd.read_csv(os.path.join(output_path, "contacts_clean.csv"), sep=",")
    all_sorts = pd.read_csv(os.path.join(tmp_path, "all_sorts_clean.csv"), sep=",")
    print(os.getcwd())

    supergroups = {}
    sgi = 0
    with open(os.path.join(tmp_path, "super_groups.csv")) as sgf:
        lines = sgf.readlines()
        for l in lines:
            for g in l.split(","):
                g = g.replace("-", "_").replace(" ", "_").rstrip()
                if g:
                    supergroups[g] = "supergroup_{}".format(sgi)
            sgi += 1
    sg_set = set(supergroups.values())

    empty_fm = []

    for indx, afm in all_sorts.iterrows():
        foundcontact = False
        for indx2, acontact in contacts.iterrows():
            if acontact["formation"] in afm["code"]:
                foundcontact = True
                break
        foundorientation = False
        for indx3, ano in orientations.iterrows():
            if ano["formation"] in afm["code"]:
                foundorientation = True
                break
        if not foundcontact or not foundorientation:
            empty_fm.append(afm["code"])

    # print(empty_fm)
    asc = pd.read_csv(os.path.join(tmp_path, "all_sorts_clean.csv"), sep=",")

    all_sorts = np.genfromtxt(
        os.path.join(tmp_path, "all_sorts_clean.csv"), delimiter=",", dtype="U100"
    )
    nformations = len(all_sorts)

    f.write("#---------------------------------------------------------------\n")
    f.write("#-----------------------Create Formations-----------------------\n")
    f.write("#---------------------------------------------------------------\n")

    for ind, row in asc.iterrows():
        if not row["code"] in empty_fm:
            f.write("GeomodellerTask {\n")
            f.write("CreateFormation {\n")

            ostr = '    name: "' + row["code"].replace("\n", "") + '"\n'
            f.write(ostr)
            r, g, b = m2l_utils.hextoints(row["colour"])
            ostr = "    red: " + str(int(r)) + "\n"
            f.write(ostr)

            ostr = "    green: " + str(int(g)) + "\n"
            f.write(ostr)

            ostr = "    blue: " + str(int(b)) + "\n"
            f.write(ostr)

            f.write("    }\n")
            f.write("}\n")

    if True:  # Basal formation if only one supergroup
        f.write("GeomodellerTask {\n")
        f.write("CreateFormation {\n")

        ostr = '    name: "Basal_Formation"\n'
        f.write(ostr)
        r, g, b = m2l_utils.hextoints(row["colour"])
        ostr = "    red: " + str(int(128)) + "\n"
        f.write(ostr)

        ostr = "    green: " + str(int(128)) + "\n"
        f.write(ostr)

        ostr = "    blue: " + str(int(128)) + "\n"
        f.write(ostr)

        f.write("    }\n")
        f.write("}\n")

    f.write("#---------------------------------------------------------------\n")
    f.write("#-----------------------Set Stratigraphic Pile------------------\n")
    f.write("#---------------------------------------------------------------\n")

    for sg in sg_set:
        # for i in range (nformations-1,0,-1):
        f.write("GeomodellerTask {\n")
        f.write("SetSeries {\n")

        ostr = '    name: "' + sg + '"\n'
        f.write(ostr)

        ostr = "    position: 1\n"
        f.write(ostr)

        ostr = '    relation: "erode"\n'
        f.write(ostr)

        f.write("    }\n")
        f.write("}\n")

        for j in range(nformations - 1, 0, -1):
            #        for j in range(1,nformations):
            if supergroups[all_sorts[j, 5]] == sg:
                if not all_sorts[j][4] in empty_fm:
                    f.write("GeomodellerTask {\n")
                    f.write("AddFormationToSeries {\n")

                    ostr = '    series: "' + sg + '"\n'
                    f.write(ostr)

                    ostr = '    formation: "' + all_sorts[j][4] + '"\n'
                    f.write(ostr)

                    f.write("    }\n")
                    f.write("}\n")

    if True:  # Basal series/formation if only one supergroup
        f.write("GeomodellerTask {\n")
        f.write("SetSeries {\n")

        ostr = '    name: "Basal_Series"\n'
        f.write(ostr)

        ostr = "    position: 1\n"
        f.write(ostr)

        ostr = '    relation: "erode"\n'
        f.write(ostr)

        f.write("    }\n")
        f.write("}\n")

        f.write("GeomodellerTask {\n")
        f.write("AddFormationToSeries {\n")

        ostr = '    series: "Basal_Series"\n'
        f.write(ostr)

        ostr = '    formation: "Basal_Formation"\n'
        f.write(ostr)

        f.write("    }\n")
        f.write("}\n")

    if save_faults:
        output_path = os.path.join(test_data_path, "output")

        faults_len = pd.read_csv(os.path.join(output_path, "fault_dimensions.csv"))

        n_allfaults = len(faults_len)

        fcount = 0
        for i in range(0, n_allfaults):
            f.write("GeomodellerTask {\n")
            f.write("CreateFault {\n")
            r, g, b = m2l_utils.hextoints(str(faults_len.iloc[i]["colour"]))
            ostr = '    name: "' + faults_len.iloc[i]["Fault"] + '"\n'
            f.write(ostr)

            ostr = "    red: " + str(r) + "\n"
            f.write(ostr)

            ostr = "    green: " + str(g) + "\n"
            f.write(ostr)

            ostr = "    blue: " + str(b) + "\n"
            f.write(ostr)

            f.write("    }\n")
            f.write("}\n")
            fcount = fcount + 1

            f.write("GeomodellerTask {\n")
            f.write("    Set3dFaultLimits {\n")
            f.write('        Fault_name: "' + faults_len.iloc[i]["Fault"] + '"\n')
            f.write(
                "        Horizontal: "
                + str(faults_len.iloc[i]["HorizontalRadius"])
                + "\n"
            )
            f.write(
                "        Vertical: " + str(faults_len.iloc[i]["VerticalRadius"]) + "\n"
            )
            f.write(
                "        InfluenceDistance: "
                + str(faults_len.iloc[i]["InfluenceDistance"])
                + "\n"
            )
            f.write("    }\n")
            f.write("}\n")

    f.write("#---------------------------------------------------------------\n")
    f.write("#-----------------------Import 3D contact data ---Base Model----\n")
    f.write("#---------------------------------------------------------------\n")

    contacts = pd.read_csv(os.path.join(output_path, "contacts_clean.csv"), sep=",")
    all_sorts = pd.read_csv(os.path.join(tmp_path, "all_sorts_clean.csv"), sep=",")
    # all_sorts.set_index('code',  inplace = True)
    # display(all_sorts)

    for inx, afm in all_sorts.iterrows():
        # print(afm[0])
        if not afm["code"] in empty_fm:
            f.write("GeomodellerTask {\n")
            f.write("    Add3DInterfacesToFormation {\n")
            f.write('          formation: "' + str(afm["code"]) + '"\n')

            for indx2, acontact in contacts.iterrows():
                if acontact["formation"] in afm["code"]:
                    ostr = (
                        "              point {x:"
                        + str(acontact["X"])
                        + "; y:"
                        + str(acontact["Y"])
                        + "; z:"
                        + str(acontact["Z"])
                        + "}\n"
                    )
                    f.write(ostr)
            f.write("    }\n")
            f.write("}\n")
    f.write("#---------------------------------------------------------------\n")
    f.write("#------------------Import 3D orientation data ---Base Model-----\n")
    f.write("#---------------------------------------------------------------\n")

    orientations = pd.read_csv(
        os.path.join(output_path, "orientations_clean.csv"), sep=","
    )
    all_sorts = pd.read_csv(os.path.join(tmp_path, "all_sorts_clean.csv"), sep=",")
    # all_sorts.set_index('code',  inplace = True)
    # display(all_sorts)

    for inx, afm in all_sorts.iterrows():
        # print(groups[agp])
        if not afm["code"] in empty_fm:
            f.write("GeomodellerTask {\n")
            f.write("    Add3DFoliationToFormation {\n")
            f.write('          formation: "' + str(afm["code"]) + '"\n')
            for indx2, ano in orientations.iterrows():
                if ano["formation"] in afm["code"]:
                    f.write("           foliation {\n")
                    ostr = (
                        "                  Point3D {x:"
                        + str(ano["X"])
                        + "; y:"
                        + str(ano["Y"])
                        + "; z:"
                        + str(ano["Z"])
                        + "}\n"
                    )
                    f.write(ostr)
                    ostr = "                  direction: " + str(ano["azimuth"]) + "\n"
                    f.write(ostr)
                    ostr = "                  dip: " + str(ano["dip"]) + "\n"
                    f.write(ostr)
                    if ano["polarity"] == 1:
                        ostr = "                  polarity: Normal_Polarity\n"
                    else:
                        ostr = "                  polarity: Reverse_Polarity\n"
                    f.write(ostr)
                    ostr = "           }\n"
                    f.write(ostr)
            f.write("    }\n")
            f.write("}\n")

    f.write("#---------------------------------------------------------------\n")
    f.write("#-----------------------Import 3D fault data ---Base Model------\n")
    f.write("#---------------------------------------------------------------\n")

    contacts = pd.read_csv(os.path.join(output_path, "faults.csv"), sep=",")
    faults = pd.read_csv(os.path.join(output_path, "fault_dimensions.csv"), sep=",")

    for indx, afault in faults.iterrows():
        f.write("GeomodellerTask {\n")
        f.write("    Add3DInterfacesToFormation {\n")
        f.write('          formation: "' + str(afault["Fault"]) + '"\n')
        for indx2, acontact in contacts.iterrows():
            if acontact["formation"] == afault["Fault"]:
                ostr = (
                    "              point {x:"
                    + str(acontact["X"])
                    + "; y:"
                    + str(acontact["Y"])
                    + "; z:"
                    + str(acontact["Z"])
                    + "}\n"
                )
                f.write(ostr)
        f.write("    }\n")
        f.write("}\n")

    f.write("#---------------------------------------------------------------\n")
    f.write("#------------------Import 3D fault orientation data ------------\n")
    f.write("#---------------------------------------------------------------\n")

    orientations = pd.read_csv(
        os.path.join(output_path, "fault_orientations.csv"), sep=","
    )
    faults = pd.read_csv(os.path.join(output_path, "fault_dimensions.csv"), sep=",")

    for indx, afault in faults.iterrows():
        f.write("GeomodellerTask {\n")
        f.write("    Add3DFoliationToFormation {\n")
        f.write('          formation: "' + str(afault["Fault"]) + '"\n')
        for indx2, ano in orientations.iterrows():
            if ano["formation"] == afault["Fault"]:
                f.write("           foliation {\n")
                ostr = (
                    "                  Point3D {x:"
                    + str(ano["X"])
                    + "; y:"
                    + str(ano["Y"])
                    + "; z:"
                    + str(ano["Z"])
                    + "}\n"
                )
                f.write(ostr)
                ostr = "                  direction: " + str(ano["DipDirection"]) + "\n"
                f.write(ostr)
                if ano["dip"] == -999:
                    ostr = (
                        "                  dip: " + str(random.randint(60, 90)) + "\n"
                    )
                else:
                    ostr = "                  dip: " + str(ano["dip"]) + "\n"
                f.write(ostr)
                if ano["DipPolarity"] == 1:
                    ostr = "                  polarity: Normal_Polarity\n"
                else:
                    ostr = "                  polarity: Reverse_Polarity\n"
                f.write(ostr)
                ostr = "           }\n"
                f.write(ostr)
        f.write("    }\n")
        f.write("}\n")

    if save_faults:
        G = nx.read_gml(os.path.join(tmp_path, "fault_network.gml"), label="label")
        # nx.draw(G, with_labels=True, font_weight='bold')
        edges = list(G.edges)
        # for i in range(0,len(edges)):
        # print(edges[i][0],edges[i][1])
        cycles = list(nx.simple_cycles(G))
        # display(cycles)
        f.write("#---------------------------------------------------------------\n")
        f.write("#-----------------------Link faults with faults ----------------\n")
        f.write("#---------------------------------------------------------------\n")
        f.write("GeomodellerTask {\n")
        f.write("    LinkFaultsWithFaults {\n")

        for i in range(0, len(edges)):
            found = False
            for j in range(0, len(cycles)):
                if edges[i][0] == cycles[j][0] and edges[i][1] == cycles[j][1]:
                    found = True  # fault pair is first two elements in a cycle list so don't save to taskfile
            if not found:
                ostr = (
                    '        FaultStopsOnFaults{ fault: "'
                    + edges[i][1]
                    + '"; stopson: "'
                    + edges[i][0]
                    + '"}\n'
                )
                f.write(ostr)

        f.write("    }\n")
        f.write("}\n")

    if save_faults:
        all_fault_group = np.genfromtxt(
            os.path.join(output_path, "group-fault-relationships.csv"),
            delimiter=",",
            dtype="U100",
        )
        ngroups = len(all_fault_group)
        all_fault_group = np.transpose(all_fault_group)
        nfaults = len(all_fault_group)

        f.write("#---------------------------------------------------------------\n")
        f.write("#-----------------------Link series with faults ----------------\n")
        f.write("#---------------------------------------------------------------\n")
        f.write("GeomodellerTask {\n")
        f.write("    LinkFaultsWithSeries {\n")

        for i in range(1, nfaults):
            first = True
            for j in range(1, ngroups):
                if all_fault_group[i, j] == str(1):
                    if first:
                        ostr = (
                            '    FaultSeriesLinks{ fault: "'
                            + all_fault_group[i, 0]
                            + '"; series: ['
                        )
                        f.write(ostr)
                        ostr = '"' + supergroups[all_fault_group[0, j]] + '"'
                        f.write(ostr)
                        last_sg = supergroups[all_fault_group[0, j]]
                        first = False
                    else:
                        if not supergroups[all_fault_group[0, j]] == last_sg:
                            ostr = ', "' + supergroups[all_fault_group[0, j]] + '"'
                            last_sg = supergroups[all_fault_group[0, j]]
                            f.write(ostr)
            if not first:
                ostr = "]}\n"
                f.write(ostr)

        f.write("    }\n")
        f.write("}\n")

    f.write("GeomodellerTask {\n")
    f.write("    SaveProjectAs {\n")
    f.write('        filename: "./' + model_name + '.xml"\n')
    f.write("    }\n")
    f.write("}\n")
    f.close()

    if compute_etc:
        f = open(os.join.path(test_data_path, model_name, "m2l_compute.taskfile"), "w")
        f.write("#---------------------------------------------------------------\n")
        f.write("#----------------------------Load Model----------------------\n")
        f.write("#---------------------------------------------------------------\n")
        f.write("GeomodellerTask {\n")
        f.write("    OpenProjectNoGUI {\n")
        f.write('        filename: "./' + model_name + '.xml"\n')
        f.write("    }\n")
        f.write("}\n")

        f.write("#---------------------------------------------------------------\n")
        f.write("#----------------------------Compute Model----------------------\n")
        f.write("#---------------------------------------------------------------\n")
        f.write("\n")
        f.write("GeomodellerTask {\n")
        f.write("    ComputeModel {\n")
        f.write("        SeriesList {\n")
        f.write('            node: "All" \n')
        f.write("        }\n")
        f.write("        SectionList {\n")
        f.write('            node: "All"\n')
        f.write("        }\n")
        f.write("        FaultList {\n")
        f.write('            node: "All"\n')
        f.write("        }\n")
        f.write("        radius: 10.0\n")
        f.write("    }\n")
        f.write("}\n")

        f.write("#---------------------------------------------------------------\n")
        f.write("#-----------------------Add geophysical Properties--------------\n")
        f.write("#---------------------------------------------------------------\n")
        f.write("\n")
        f.write("\n")
        f.write("\n")
        f.write("#---------------------------------------------------------------\n")
        f.write("#--------------------------Export Lithology Voxet---------------\n")
        f.write("#---------------------------------------------------------------\n")
        f.write("GeomodellerTask {\n")
        f.write("    SaveLithologyVoxet {\n")
        f.write("        nx: 25\n")
        f.write("        ny: 25\n")
        f.write("        nz: 40\n")
        f.write('        LithologyVoxetFileStub: "./Litho_Voxet/LithoVoxet.vo"\n')
        f.write("    }\n")
        f.write("}\n")
        f.write("#---------------------------------------------------------------\n")
        f.write("#--------------------------Save As Model------------------------\n")
        f.write("#---------------------------------------------------------------\n")
        f.write("\n")

        f.write("GeomodellerTask {\n")
        f.write("    SaveProjectAs {\n")
        f.write('        filename: "/' + model_name + '.xml"\n')
        f.write("    }\n")
        f.write("}\n")
        f.write("GeomodellerTask {\n")
        f.write("    CloseProjectNoGUI {\n")
        f.write("    }\n")
        f.write("}\n")

        f.close()


# same same expect it builds a list that then gets written all at once (this version is slower!)
def loop2geomodeller2(
    model_name,
    test_data_path,
    tmp_path,
    output_path,
    dtm_file,
    bbox,
    save_faults,
    compute_etc,
    workflow,
):

    f = open(os.path.join(test_data_path, model_name, "m2l.taskfile"), "w")
    ostr = []

    ostr.append("#---------------------------------------------------------------\n")
    ostr.append("#-----------------------Project Header-----------------------\n")
    ostr.append("#---------------------------------------------------------------\n")
    ostr.append('name: "UWA_Intrepid"\n')
    ostr.append('description: "Automate_batch_Model"\n')
    ostr.append("    GeomodellerTask {\n")
    ostr.append("    CreateProject {\n")
    ostr.append('        name: "Hamersley"\n')
    ostr.append('        author: "Mark"\n')
    ostr.append('        date: "23/10/2019  0: 0: 0"\n')
    ostr.append('        projection { map_projection: "GDA94 / MGA50"}\n')
    ostr.append('        version: "2.0"\n')
    ostr.append("        units: meters\n")
    ostr.append("        precision: 1.0\n")
    ostr.append("        Extents {\n")
    ostr.append("            xmin: " + str(bbox[0]) + "\n")
    ostr.append("            ymin: " + str(bbox[1]) + "\n")
    ostr.append("            zmin: -7000\n")
    ostr.append("            xmax: " + str(bbox[2]) + "\n")
    ostr.append("            ymax: " + str(bbox[3]) + "\n")
    ostr.append("            zmax: 1200\n")
    ostr.append("        }\n")
    ostr.append("        deflection2d: 0.001\n")
    ostr.append("        deflection3d: 0.001\n")
    ostr.append("        discretisation: 10.0\n")
    ostr.append("        referenceTop: false\n")
    ostr.append("        CustomDTM {\n")
    ostr.append("            Extents {\n")
    ostr.append("            xmin: " + str(bbox[0]) + "\n")
    ostr.append("            ymin: " + str(bbox[1]) + "\n")
    ostr.append("            xmax: " + str(bbox[2]) + "\n")
    ostr.append("            ymax: " + str(bbox[3]) + "\n")
    ostr.append("            }\n")
    ostr.append('            name: "Topography"\n')
    ostr.append("            filename {\n")
    ostr.append('                Grid_Name: "' + dtm_file + '"\n')
    ostr.append("            }\n")
    ostr.append("            nx: 10\n")
    ostr.append("            ny: 10\n")
    ostr.append("        }\n")
    ostr.append("    }\n")
    ostr.append("}\n")

    orientations = pd.read_csv(
        os.path.join(output_path, "orientations_clean.csv"), sep=","
    )
    contacts = pd.read_csv(os.path.join(output_path, "contacts_clean.csv"), sep=",")
    all_sorts = pd.read_csv(os.path.join(tmp_path, "all_sorts_clean.csv"), sep=",")

    empty_fm = []

    for indx, afm in all_sorts.iterrows():
        foundcontact = False
        for indx2, acontact in contacts.iterrows():
            if acontact["formation"] in afm["code"]:
                foundcontact = True
                break
        foundorientation = False
        for indx3, ano in orientations.iterrows():
            if ano["formation"] in afm["code"]:
                foundorientation = True
                break
        if not foundcontact or not foundorientation:
            empty_fm.append(afm["code"])

    # print(empty_fm)

    all_sorts = np.genfromtxt(
        os.path.join(tmp_path, "all_sorts_clean.csv"), delimiter=",", dtype="U100"
    )
    nformations = len(all_sorts)

    ostr.append("#---------------------------------------------------------------\n")
    ostr.append("#-----------------------Create Formations-----------------------\n")
    ostr.append("#---------------------------------------------------------------\n")

    for i in range(1, nformations):
        if not all_sorts[i, 4] in empty_fm:
            ostr.append("GeomodellerTask {\n")
            ostr.append("CreateFormation {\n")

            ostr2 = '    name: "' + all_sorts[i, 4].replace("\n", "") + '"\n'
            ostr.append(ostr2)

            ostr2 = "    red: " + str(random.randint(1, 256) - 1) + "\n"
            ostr.append(ostr2)

            ostr2 = "    green: " + str(random.randint(1, 256) - 1) + "\n"
            ostr.append(ostr2)

            ostr2 = "    blue: " + str(random.randint(1, 256) - 1) + "\n"
            ostr.append(ostr2)

            ostr.append("    }\n")
            ostr.append("}\n")

    ostr.append("#---------------------------------------------------------------\n")
    ostr.append("#-----------------------Set Stratigraphic Pile------------------\n")
    ostr.append("#---------------------------------------------------------------\n")

    for i in range(1, nformations):
        # for i in range (nformations-1,0,-1):
        if all_sorts[i, 2] == str(1):
            ostr.append("GeomodellerTask {\n")
            ostr.append("SetSeries {\n")

            ostr2 = '    name: "' + all_sorts[i][5].replace("\n", "") + '"\n'
            ostr.append(ostr2)

            ostr2 = "    position: 1\n"
            ostr.append(ostr2)

            ostr2 = '    relation: "erode"\n'
            ostr.append(ostr2)

            ostr.append("    }\n")
            ostr.append("}\n")

            for j in range(nformations - 1, 0, -1):
                #        for j in range(1,nformations):
                if all_sorts[j, 1] == all_sorts[i, 1]:
                    if not all_sorts[j][4] in empty_fm:
                        ostr.append("GeomodellerTask {\n")
                        ostr.append("AddFormationToSeries {\n")

                        ostr2 = '    series: "' + all_sorts[j][5] + '"\n'
                        ostr.append(ostr2)

                        ostr2 = '    formation: "' + all_sorts[j][4] + '"\n'
                        ostr.append(ostr2)

                        ostr.append("    }\n")
                        ostr.append("}\n")

    if save_faults:
        output_path = os.path.join(test_data_path, "output")

        faults_len = pd.read_csv(os.path.join(output_path, "fault_dimensions.csv"))

        n_allfaults = len(faults_len)

        fcount = 0
        for i in range(0, n_allfaults):
            ostr.append("GeomodellerTask {\n")
            ostr.append("CreateFault {\n")
            ostr2 = '    name: "' + faults_len.iloc[i]["Fault"] + '"\n'
            ostr.append(ostr2)

            ostr2 = "    red: " + str(random.randint(1, 256) - 1) + "\n"
            ostr.append(ostr2)

            ostr2 = "    green: " + str(random.randint(1, 256) - 1) + "\n"
            ostr.append(ostr2)

            ostr2 = "    blue: " + str(random.randint(1, 256) - 1) + "\n"
            ostr.append(ostr2)

            ostr.append("    }\n")
            ostr.append("}\n")
            fcount = fcount + 1

            ostr.append("GeomodellerTask {\n")
            ostr.append("    Set3dFaultLimits {\n")
            ostr.append('        Fault_name: "' + faults_len.iloc[i]["Fault"] + '"\n')
            ostr.append(
                "        Horizontal: "
                + str(faults_len.iloc[i]["HorizontalRadius"])
                + "\n"
            )
            ostr.append(
                "        Vertical: " + str(faults_len.iloc[i]["VerticalRadius"]) + "\n"
            )
            ostr.append(
                "        InfluenceDistance: "
                + str(faults_len.iloc[i]["InfluenceDistance"])
                + "\n"
            )
            ostr.append("    }\n")
            ostr.append("}\n")

    ostr.append("#---------------------------------------------------------------\n")
    ostr.append("#-----------------------Import 3D contact data ---Base Model----\n")
    ostr.append("#---------------------------------------------------------------\n")

    contacts = pd.read_csv(os.path.join(output_path, "contacts_clean.csv"), sep=",")
    all_sorts = pd.read_csv(os.path.join(tmp_path, "all_sorts_clean.csv"), sep=",")
    # all_sorts.set_index('code',  inplace = True)
    # display(all_sorts)

    for inx, afm in all_sorts.iterrows():
        # print(afm[0])
        if not afm["code"] in empty_fm:
            ostr.append("GeomodellerTask {\n")
            ostr.append("    Add3DInterfacesToFormation {\n")
            ostr.append('          formation: "' + str(afm["code"]) + '"\n')

            for indx2, acontact in contacts.iterrows():
                if acontact["formation"] in afm["code"]:
                    ostr2 = (
                        "              point {x:"
                        + str(acontact["X"])
                        + "; y:"
                        + str(acontact["Y"])
                        + "; z:"
                        + str(acontact["Z"])
                        + "}\n"
                    )
                    ostr.append(ostr2)
            ostr.append("    }\n")
            ostr.append("}\n")
    ostr.append("#---------------------------------------------------------------\n")
    ostr.append("#------------------Import 3D orientation data ---Base Model-----\n")
    ostr.append("#---------------------------------------------------------------\n")

    orientations = pd.read_csv(
        os.path.join(output_path, "orientations_clean.csv"), sep=","
    )
    all_sorts = pd.read_csv(os.path.join(tmp_path, "all_sorts_clean.csv"), sep=",")
    # all_sorts.set_index('code',  inplace = True)
    # display(all_sorts)

    for inx, afm in all_sorts.iterrows():
        # print(groups[agp])
        if not afm["code"] in empty_fm:
            ostr.append("GeomodellerTask {\n")
            ostr.append("    Add3DFoliationToFormation {\n")
            ostr.append('          formation: "' + str(afm["code"]) + '"\n')
            for indx2, ano in orientations.iterrows():
                if ano["formation"] in afm["code"]:
                    ostr.append("           foliation {\n")
                    ostr2 = (
                        "                  Point3D {x:"
                        + str(ano["X"])
                        + "; y:"
                        + str(ano["Y"])
                        + "; z:"
                        + str(ano["Z"])
                        + "}\n"
                    )
                    ostr.append(ostr2)
                    ostr2 = "                  direction: " + str(ano["azimuth"]) + "\n"
                    ostr.append(ostr2)
                    ostr2 = "                  dip: " + str(ano["dip"]) + "\n"
                    ostr.append(ostr2)
                    if ano["polarity"] == 1:
                        ostr2 = "                  polarity: Normal_Polarity\n"
                    else:
                        ostr2 = "                  polarity: Reverse_Polarity\n"
                    ostr.append(ostr2)
                    ostr2 = "           }\n"
                    ostr.append(ostr2)
            ostr.append("    }\n")
            ostr.append("}\n")

    ostr.append("#---------------------------------------------------------------\n")
    ostr.append("#-----------------------Import 3D fault data ---Base Model------\n")
    ostr.append("#---------------------------------------------------------------\n")

    contacts = pd.read_csv(os.path.join(output_path, "faults.csv"), sep=",")
    faults = pd.read_csv(os.path.join(output_path, "fault_dimensions.csv"), sep=",")

    for indx, afault in faults.iterrows():
        ostr.append("GeomodellerTask {\n")
        ostr.append("    Add3DInterfacesToFormation {\n")
        ostr.append('          formation: "' + str(afault["Fault"]) + '"\n')
        for indx2, acontact in contacts.iterrows():
            if acontact["formation"] == afault["Fault"]:
                ostr2 = (
                    "              point {x:"
                    + str(acontact["X"])
                    + "; y:"
                    + str(acontact["Y"])
                    + "; z:"
                    + str(acontact["Z"])
                    + "}\n"
                )
                ostr.append(ostr2)
        ostr.append("    }\n")
        ostr.append("}\n")

    ostr.append("#---------------------------------------------------------------\n")
    ostr.append("#------------------Import 3D fault orientation data ------------\n")
    ostr.append("#---------------------------------------------------------------\n")

    orientations = pd.read_csv(
        os.path.join(output_path, "fault_orientations.csv"), sep=","
    )
    faults = pd.read_csv(os.path.join(output_path, "fault_dimensions.csv"), sep=",")

    for indx, afault in faults.iterrows():
        ostr.append("GeomodellerTask {\n")
        ostr.append("    Add3DFoliationToFormation {\n")
        ostr.append('          formation: "' + str(afault["Fault"]) + '"\n')
        for indx2, ano in orientations.iterrows():
            if ano["formation"] == afault["Fault"]:
                ostr.append("           foliation {\n")
                ostr2 = (
                    "                  Point3D {x:"
                    + str(ano["X"])
                    + "; y:"
                    + str(ano["Y"])
                    + "; z:"
                    + str(ano["Z"])
                    + "}\n"
                )
                ostr.append(ostr2)
                ostr2 = (
                    "                  direction: " + str(ano["DipDirection"]) + "\n"
                )
                ostr.append(ostr2)
                if ano["dip"] == -999:
                    ostr2 = (
                        "                  dip: " + str(random.randint(60, 90)) + "\n"
                    )
                else:
                    ostr2 = "                  dip: " + str(ano["dip"]) + "\n"
                ostr.append(ostr2)
                if ano["DipPolarity"] == 1:
                    ostr2 = "                  polarity: Normal_Polarity\n"
                else:
                    ostr2 = "                  polarity: Reverse_Polarity\n"
                ostr.append(ostr2)
                ostr2 = "           }\n"
                ostr.append(ostr2)
        ostr.append("    }\n")
        ostr.append("}\n")

    if save_faults:
        G = nx.read_gml(os.path.join(tmp_path, "fault_network.gml"), label="label")
        # nx.draw(G, with_labels=True, font_weight='bold')
        edges = list(G.edges)
        # for i in range(0,len(edges)):
        # print(edges[i][0],edges[i][1])
        cycles = list(nx.simple_cycles(G))
        # display(cycles)
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append(
            "#-----------------------Link faults with faults ----------------\n"
        )
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append("GeomodellerTask {\n")
        ostr.append("    LinkFaultsWithFaults {\n")

        for i in range(0, len(edges)):
            found = False
            for j in range(0, len(cycles)):
                if edges[i][0] == cycles[j][0] and edges[i][1] == cycles[j][1]:
                    found = True  # fault pair is first two elements in a cycle list so don't save to taskfile
            if not found:
                ostr2 = (
                    '        FaultStopsOnFaults{ fault: "'
                    + edges[i][1]
                    + '"; stopson: "'
                    + edges[i][0]
                    + '"}\n'
                )
                ostr.append(ostr2)

        ostr.append("    }\n")
        ostr.append("}\n")

    if save_faults:
        all_fault_group = np.genfromtxt(
            os.path.join(output_path, "group-fault-relationships.csv"),
            delimiter=",",
            dtype="U100",
        )
        ngroups = len(all_fault_group)
        all_fault_group = np.transpose(all_fault_group)
        nfaults = len(all_fault_group)

        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append(
            "#-----------------------Link series with faults ----------------\n"
        )
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append("GeomodellerTask {\n")
        ostr.append("    LinkFaultsWithSeries {\n")

        for i in range(1, nfaults):
            first = True
            for j in range(1, ngroups):
                if all_fault_group[i, j] == str(1):
                    if first:
                        ostr2 = (
                            '    FaultSeriesLinks{ fault: "'
                            + all_fault_group[i, 0]
                            + '"; series: ['
                        )
                        ostr.append(ostr2)
                        ostr2 = '"' + all_fault_group[0, j] + '"'
                        ostr.append(ostr2)
                        first = False
                    else:
                        ostr2 = ', "' + all_fault_group[0, j] + '"'
                        ostr.append(ostr2)
            if not first:
                ostr2 = "]}\n"
                ostr.append(ostr2)

        ostr.append("    }\n")
        ostr.append("}\n")

    ostr.append("GeomodellerTask {\n")
    ostr.append("    SaveProjectAs {\n")
    ostr.append('        filename: "./' + model_name + '.xml"\n')
    ostr.append("    }\n")
    ostr.append("}\n")
    f.close()

    if compute_etc:
        f = open(os.path.join(test_data_path, model_name, "m2l_compute.taskfile"), "w")
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append("#----------------------------Load Model----------------------\n")
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append("GeomodellerTask {\n")
        ostr.append("    OpenProjectNoGUI {\n")
        ostr.append('        filename: "./' + model_name + '.xml"\n')
        ostr.append("    }\n")
        ostr.append("}\n")

        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append(
            "#----------------------------Compute Model----------------------\n"
        )
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append("\n")
        ostr.append("GeomodellerTask {\n")
        ostr.append("    ComputeModel {\n")
        ostr.append("        SeriesList {\n")
        ostr.append('            node: "All" \n')
        ostr.append("        }\n")
        ostr.append("        SectionList {\n")
        ostr.append('            node: "All"\n')
        ostr.append("        }\n")
        ostr.append("        FaultList {\n")
        ostr.append('            node: "All"\n')
        ostr.append("        }\n")
        ostr.append("        radius: 10.0\n")
        ostr.append("    }\n")
        ostr.append("}\n")

        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append(
            "#-----------------------Add geophysical Properties--------------\n"
        )
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append("\n")
        ostr.append("\n")
        ostr.append("\n")
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append(
            "#--------------------------Export Lithology Voxet---------------\n"
        )
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append("GeomodellerTask {\n")
        ostr.append("    SaveLithologyVoxet {\n")
        ostr.append("        nx: 25\n")
        ostr.append("        ny: 25\n")
        ostr.append("        nz: 40\n")
        ostr.append('        LithologyVoxetFileStub: "./Litho_Voxet/LithoVoxet.vo"\n')
        ostr.append("    }\n")
        ostr.append("}\n")
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append(
            "#--------------------------Save As Model------------------------\n"
        )
        ostr.append(
            "#---------------------------------------------------------------\n"
        )
        ostr.append("\n")

        ostr.append("GeomodellerTask {\n")
        ostr.append("    SaveProjectAs {\n")
        ostr.append('        filename: "/' + model_name + '.xml"\n')
        ostr.append("    }\n")
        ostr.append("}\n")
        ostr.append("GeomodellerTask {\n")
        ostr.append("    CloseProjectNoGUI {\n")
        ostr.append("    }\n")
        ostr.append("}\n")
        f.writelines(ostr)
        f.close()


##########################################################################
# Import outputs from map2loop to LoopStructural and view with Lavavu
#
# loop2LoopStructural(thickness_file,orientation_file,contacts_file,bbox)
# Args:
# bbox model bounding box
#
# Calculates model and displays in LavaVu wthin notebook
##########################################################################
def loop2LoopStructural(m2l_directory):
    """create a model from a map2loop directory

    [extended_summary]

    Parameters
    ----------
    m2l_directory : string
        path to the map2loop directory
    """
    visualise = False
    # make sure everything is installed and can be imported
    try:
        from LoopStructural import GeologicalModel
        from LoopStructural.utils import process_map2loop
    except ImportError:
        print("Loop Structural not installed")
        return
    try:
        from LoopStructural.visualisation import LavaVuModelViewer

        visualise = True
    except ImportError:
        print(
            "Lavavu is not installed, try installing it with pip \n"
            "Model will be built but cannot be visualised"
        )

    m2l_data = process_map2loop(m2l_directory)
    boundary_points = np.zeros((2, 3))
    boundary_points[0, 0] = m2l_data["bounding_box"]["minx"]
    boundary_points[0, 1] = m2l_data["bounding_box"]["miny"]
    boundary_points[0, 2] = m2l_data["bounding_box"]["lower"]
    boundary_points[1, 0] = m2l_data["bounding_box"]["maxx"]
    boundary_points[1, 1] = m2l_data["bounding_box"]["maxy"]
    boundary_points[1, 2] = m2l_data["bounding_box"]["upper"]

    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.set_model_data(m2l_data["data"])

    faults = []
    for f in m2l_data["max_displacement"].keys():
        if model.data[model.data["type"] == f].shape[0] == 0:
            continue
        fault_id = f[6:]
        overprints = []
        try:
            overprint_id = m2l_data["fault_fault"][
                m2l_data["fault_fault"][fault_id] == 1
            ]["fault_id"].to_numpy()
            for i in overprint_id:
                overprints.append(["Fault_%i" % i])
        except:
            print("No entry for %s in fault_fault_relations" % f)
        #     continue
        faults.append(
            model.create_and_add_fault(
                f,
                -m2l_data["max_displacement"][f],
                faultfunction="BaseFault",
                interpolatortype="FDI",
                nelements=1e4,
                # data_region=.1, Not implemented currently
                #                                                  regularisation=[1,1,1],
                solver="pyamg",
                #                                                  damp=True,
                #                                                  buffer=0.1,
                #                                                  steps=1,
                overprints=overprints,
                cpw=10,
                npw=10,
            )
        )

    # loop through all of the groups and add them to the model in youngest to oldest.
    group_features = []
    for i in m2l_data["groups"]["group number"].unique():
        g = (
            m2l_data["groups"]
            .loc[m2l_data["groups"]["group number"] == i, "group"]
            .unique()[0]
        )
        group_features.append(
            model.create_and_add_foliation(
                g,
                interpolatortype="PLI",  # which interpolator to use
                nelements=1e5,  # how many tetras/voxels
                buffer=0.5,  # how much to extend nterpolation around box
                solver="pyamg",
                damp=True,
            )
        )
        # if the group was successfully added (not null) then lets add the base (0 to be unconformity)
        if group_features[-1]:
            model.add_unconformity(group_features[-1]["feature"], 0)
    model.set_stratigraphic_column(m2l_data["stratigraphic_column"])
    if visualise:
        viewer = LavaVuModelViewer(model)
        viewer.add_model(cmap="tab20")
        viewer.interactive()


##########################################################################
# Import outputs from map2loop to gempy and view with pyvtk
# loop2gempy(test_data_name,tmp_path,vtk_pth,orientations_file,contacts_file,groups_file,dtm_reproj_file,bbox,model_base, model_top,vtk)
# Args:
# test_data_name root name of project
# tmp_path path of temp files directory
# vtk_pth path of vtk output directory
# orientations_file path of orientations file
# contacts_file path of contacts file
# groups_file path of groups file
# dtm_reproj_file path of dtm file
# bbox model bounding box
# model_base z value ofbase of model
# model_top z value of top of model
# vtk flag as to wether to save out model to vtk
#
# Calculates model and displays in external vtk viewer
##########################################################################
def loop2gempy__(
    test_data_name: str,
    tmp_path: str,
    vtk_path: str,
    orientations_file: str,
    contacts_file: str,
    groups_file: str,
    bbox: tuple,
    model_base: float,
    model_top: float,
    vtk: bool,
    dtm_reproj_file: str = None,
    va=None,
    verbose: bool = False,
    compute: bool = True,
):
    """

    :param test_data_name:
    :param tmp_path:
    :param vtk_path:
    :param orientations_file:
    :param contacts_file:
    :param groups_file:
    :param bbox:
    :param model_base:
    :param model_top:
    :param vtk:
    :param dtm_reproj_file:
    :param va: vertical anisotropy. Factor by which all Z coordinates are multiplied by
    :param verbose:
    :param compute:
    :return:
    """
    import gempy as gp
    from gempy import plot

    geo_model = gp.create_model(test_data_name)

    # If depth coordinates are much smaller than XY the whole system of equations becomes very unstable. Until
    # I fix it properly in gempy this is a handcrafted hack
    if va is None:
        va = (float(bbox[0]) - float(bbox[2])) / (model_base - model_top) / 2

        if va < 3:
            va = 0
        else:
            print("The vertical exageration is: ", va)

    gp.init_data(
        geo_model,
        extent=[bbox[0], bbox[2], bbox[1], bbox[3], model_base * va, model_top * va],
        resolution=(50, 50, 50),
        path_o=orientations_file,
        path_i=contacts_file,
    )

    geo_model.modify_surface_points(
        geo_model.surface_points.df.index, Z=geo_model.surface_points.df["Z"] * va
    )

    if dtm_reproj_file is not None:
        # Load reprojected topography to model

        fp = dtm_reproj_file
        geo_model.set_topography(source="gdal", filepath=fp)

        # Rescaling topography:
        geo_model._grid.topography.values[:, 2] *= va
        geo_model._grid.update_grid_values()
        geo_model.update_from_grid()

    # Pile processing:
    contents = np.genfromtxt(groups_file, delimiter=",", dtype="U100")

    # Init dictionary Series:Surfaces
    map_series_to_surfaces = {}
    choice = 0
    for group in contents:
        # Reading surfaces groups
        surfaces_g = np.atleast_2d(
            np.genfromtxt(
                os.path.join(tmp_path, group + ".csv"), delimiter=",", dtype="U100"
            )
        )

        # Check if there are several choices
        if surfaces_g.shape[1] > 1:
            surfaces_g = surfaces_g[choice]
        # Deleting the first element since it is not a surface
        surfaces_g = surfaces_g[1:]
        # Creating the mapping dictionary
        map_series_to_surfaces[group] = surfaces_g.tolist()

    if verbose is True:
        print(map_series_to_surfaces)

    # Setting pile to model
    gp.map_series_to_surfaces(
        geo_model, map_series_to_surfaces, remove_unused_series=False
    )

    if "Default series" in map_series_to_surfaces:

        #    Removing related data
        del_surfaces = geo_model.surfaces.df.groupby("series").get_group(
            "Default series"
        )["surface"]
        geo_model.delete_surfaces(del_surfaces, remove_data=True)

        # Removing series that have not been mapped to any surface
        geo_model.delete_series("Default series")

    if compute is True:
        gp.set_interpolator(geo_model, theano_optimizer="fast_run", dtype="float64")
        gp.compute_model(geo_model)

    # Visualise Model
    # gp.plot.plot_3D(geo_model, render_data=False)
    p3d = gp.plot_3d(geo_model, plotter_type="background", notebook=False)

    p3d3 = gp.plot_3d(geo_model, notebook=True)

    # Save model as vtk
    if vtk:
        gp.plot.export_to_vtk(
            geo_model,
            path=vtk_path,
            name=test_data_name + ".vtk",
            voxels=False,
            block=None,
            surfaces=True,
        )

    return geo_model


def loop2gempy_(
    test_data_name,
    tmp_path,
    vtk_path,
    orientations_file,
    contacts_file,
    groups_file,
    dtm_reproj_file,
    bbox,
    model_base,
    model_top,
    vtk,
):
    import gempy as gp
    from gempy import plot

    geo_model = gp.create_model(test_data_name)

    # If depth coordinates are much smaller than XY the whole system of equations becomes very unstable. Until
    # I fix it properly in gempy this is a handcrafted hack
    ve = (bbox[0] - bbox[2]) / (model_base - model_top)

    if ve < 3:
        ve = 0
    else:
        print("The vertical exageration is: ", ve)

    gp.init_data(
        geo_model,
        extent=[bbox[0], bbox[2], bbox[1], bbox[3], model_base * ve, model_top * ve],
        resolution=(50, 50, 50),
        path_o=orientations_file,
        path_i=contacts_file,
        default_values=True,
    )

    # Show example lithological points
    # gp.get_data(geo_model, 'surface_points').head()

    # Show example orientations
    # gp.get_data(geo_model, 'orientations').head()

    # Plot some of this data
    # gp.plot.plot_data(geo_model, direction='z')

    geo_model.modify_surface_points(
        geo_model.surface_points.df.index, Z=geo_model.surface_points.df["Z"] * ve
    )

    # Load reprojected topgraphy to model

    fp = dtm_reproj_file
    geo_model.set_topography(source="gdal", filepath=fp)

    contents = np.genfromtxt(groups_file, delimiter=",", dtype="U100")
    ngroups = len(contents)

    faults = gp.Faults()
    series = gp.Series(faults)
    # series.df

    # display(ngroups,contents)
    groups = []

    for i in range(0, ngroups):
        groups.append(contents[i].replace("\n", ""))
        series.add_series(contents[i].replace("\n", ""))
        print(contents[i].replace("\n", ""))

    series.delete_series("Default series")

    # series

    # Load surfaces and assign to series
    surfaces = gp.Surfaces(series)

    print(ngroups, groups)
    for i in range(0, ngroups):
        contents = np.genfromtxt(
            os.path.join(tmp_path, groups[i] + ".csv"), delimiter=",", dtype="U100"
        )
        nformations = len(contents.shape)

        if nformations == 1:
            for j in range(1, len(contents)):
                surfaces.add_surface(str(contents[j]).replace("\n", ""))
                d = {groups[i]: str(contents[j]).replace("\n", "")}
                surfaces.map_series(
                    {groups[i]: (str(contents[j]).replace("\n", ""))}
                )  # working but no gps
        else:
            for j in range(1, len(contents[0])):
                surfaces.add_surface(str(contents[0][j]).replace("\n", ""))
                d = {groups[i]: str(contents[0][j]).replace("\n", "")}
                surfaces.map_series(
                    {groups[i]: (str(contents[0][j]).replace("\n", ""))}
                )  # working but no gps

    # Set Interpolation Data
    id_only_one_bool = geo_model.surface_points.df["id"].value_counts() == 1
    id_only_one = id_only_one_bool.index[id_only_one_bool]
    single_vals = geo_model.surface_points.df[
        geo_model.surface_points.df["id"].isin(id_only_one)
    ]
    for idx, vals in single_vals.iterrows():
        geo_model.add_surface_points(vals["X"], vals["Y"], vals["Z"], vals["surface"])

    geo_model.update_structure()

    gp.set_interpolation_data(
        geo_model, compile_theano=True, theano_optimizer="fast_compile", verbose=[]
    )

    # Provide summary data on model

    # geo_model.additional_data.structure_data

    # Calculate Model
    gp.compute_model(geo_model)

    # Extract surfaces to visualize in 3D renderers
    # gp.plot.plot_section(geo_model, 49, direction='z', show_data=False)

    ver, sim = gp.get_surfaces(geo_model)

    # import winsound
    # duration = 700  # milliseconds
    # freq = 1100  # Hz
    # winsound.Beep(freq, duration)
    # winsound.Beep(freq, duration)
    # winsound.Beep(freq, duration)

    # Visualise Model
    gp.plot.plot_3D(geo_model, render_data=False)

    # Save model as vtk
    if vtk:
        gp.plot.export_to_vtk(
            geo_model,
            path=vtk_path,
            name=test_data_name + ".vtk",
            voxels=False,
            block=None,
            surfaces=True,
        )

    return geo_model


#     Courtesy of https://gist.github.com/delestro/54d5a34676a8cef7477e


def rand_cmap(
    nlabels, type="bright", first_color_black=True, last_color_black=False, verbose=True
):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib

    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    if type not in ("bright", "soft"):
        print('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print("Number of labels: " + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == "bright":
        randHSVcolors = [
            (
                np.random.uniform(low=0.0, high=1),
                np.random.uniform(low=0.2, high=1),
                np.random.uniform(low=0.9, high=1),
            )
            for i in range(nlabels)
        ]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(
                colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2])
            )

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list(
            "new_map", randRGBcolors, N=nlabels
        )

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == "soft":
        low = 0.6
        high = 0.95
        randRGBcolors = [
            (
                np.random.uniform(low=low, high=high),
                np.random.uniform(low=low, high=high),
                np.random.uniform(low=low, high=high),
            )
            for i in range(nlabels)
        ]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list(
            "new_map", randRGBcolors, N=nlabels
        )

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(
            ax,
            cmap=random_colormap,
            norm=norm,
            spacing="proportional",
            ticks=None,
            boundaries=bounds,
            format="%1i",
            orientation="horizontal",
        )

    return random_colormap


def display_LS_map(
    model, dtm, geol_clip, faults_clip, dst_crs, use_cmap, cmap, use_topo, use_faults
):

    if not use_cmap:
        cmap = rand_cmap(
            100,
            type="soft",
            first_color_black=False,
            last_color_black=False,
            verbose=False,
        )

    dtm_val = dtm.read(1)

    grid = np.array((dtm_val.shape[0] * dtm_val.shape[1], 3))
    scale = (dtm.bounds[2] - dtm.bounds[0]) / dtm_val.shape[1]
    x = np.linspace(dtm.bounds[0], dtm.bounds[2], dtm_val.shape[1])
    y = np.linspace(dtm.bounds[3], dtm.bounds[1], dtm_val.shape[0])
    xx, yy = np.meshgrid(x, y, indexing="ij")

    if use_topo:
        zz = dtm_val.flatten()
    else:
        zz = np.zeros_like(xx)

    points = np.array(
        [xx.flatten(order="F"), yy.flatten(order="F"), zz.flatten(order="F")]
    ).T
    v = model.evaluate_model(model.scale(points), scale=False)
    transform = from_origin(dtm.bounds[0], dtm.bounds[3], scale, scale)

    memfile = MemoryFile()
    new_dataset = memfile.open(
        driver="GTiff",
        height=dtm.shape[0],
        width=dtm.shape[1],
        count=1,
        dtype="float64",
        crs=dst_crs,
        transform=transform,
    )
    new_dataset.write(
        v.astype("float64").reshape(dtm_val.shape[0], dtm_val.shape[1]), 1
    )

    fig, ax = matplotlib.pyplot.subplots(figsize=(15, 15))
    rasterio.plot.show(
        new_dataset.read(1), transform=new_dataset.transform, cmap=cmap, ax=ax
    )
    geol_clip.plot(ax=ax, facecolor="none", edgecolor="black", linewidth=0.4)
    if use_faults:
        faults_clip.plot(ax=ax, facecolor="none", edgecolor="red", linewidth=0.7)


@beartype.beartype
def export_to_projectfile(loopFilename, config: Config, overwrite: bool = False):
    if loopFilename is None or loopFilename == "":
        loopFilename = os.path.join(
            config.output_path, os.path.basename(config.output_path) + ".loop3d"
        )

    if not os.path.isdir(os.path.dirname(config.output_path)):
        print(
            "(ERROR) Output path "
            + os.path.dirname(config.output_path)
            + " does not exist aborting project file export"
        )
        return

    # Check that file exist and check version matches
    overwriteVersionMismatchedFile = False
    existingExtents = None
    if os.path.isfile(loopFilename):
        resp = LoopProjectFile.Get(loopFilename, "version")
        if resp["errorFlag"]:
            print(resp["errorString"])
        else:
            if resp["value"] != LoopProjectFile.LoopVersion():
                if overwrite:
                    print("Overwriting version mismatched loop project file")
                    overwriteVersionMismatchedFile = True
                    resp = LoopProjectFile.Get(loopFilename, "extents")
                    if not resp["errorFlag"]:
                        existingExtents = resp["value"]
                    os.remove(loopFilename)
                else:
                    print(
                        "(ERROR) LoopProjectFile version mismatch, export to loop project file aborted"
                    )
                    return

    # If file does not exist or overwriting a version mismatched file create one
    if not os.path.isfile(loopFilename) or overwriteVersionMismatchedFile:
        resp = LoopProjectFile.CreateBasic(loopFilename)
        if resp["errorFlag"]:
            print(resp["errorString"])
            print("(ERROR) Cannot continue with loop project file export")
            return
        else:
            if existingExtents == None:
                # Without meaningful geodesic and utm zone values set extents to a non-conflicting
                # location in the pacific
                LoopProjectFile.Set(
                    loopFilename,
                    "extents",
                    geodesic=[-180, -179, 0, 1],
                    utm=[
                        1,
                        1,
                        config.bbox_3d["maxx"],
                        config.bbox_3d["minx"],
                        config.bbox_3d["maxy"],
                        config.bbox_3d["miny"],
                    ],
                    depth=[config.bbox_3d["top"], config.bbox_3d["base"]],
                    spacing=[1000, 1000, 500],
                    preference="utm",
                )
            else:
                resp = LoopProjectFile.Set(loopFilename, "extents", **existingExtents)
                if resp["errorFlag"]:
                    print(resp["errorString"])

    form2supergroup = pd.read_csv(
        os.path.join(config.tmp_path, "all_sorts_clean.csv"), sep=","
    )[["code", "group", "colour"]].rename(
        columns={"code": "formation", "group": "supergroup"}
    )

    # Change to formation_summary_thicknesses to get all layers and also dont' need to calc thickness again
    stratigraphicLayers = pd.read_csv(
        os.path.join(config.output_path, "formation_thicknesses.csv")
    )

    stratAges = pd.read_csv(os.path.join(config.tmp_path, "age_sorted_groups.csv"))[
        ["group_", "min", "max"]
    ]
    stratAges.rename(
        columns={"group_": "supergroup", "min": "minAge", "max": "maxAge"}, inplace=True
    )
    stratLayers = pd.merge(form2supergroup, stratigraphicLayers, on=["formation"])
    stratLayers = pd.merge(stratLayers, stratAges, on=["supergroup"])
    stratLayers["colour1Red"] = [int(a[1:3], 16) for a in stratLayers["colour"]]
    stratLayers["colour1Green"] = [int(a[3:5], 16) for a in stratLayers["colour"]]
    stratLayers["colour1Blue"] = [int(a[5:7], 16) for a in stratLayers["colour"]]
    thickness = {}
    uniqueLayers = stratLayers[
        [
            "formation",
            "supergroup",
            "colour1Red",
            "colour1Green",
            "colour1Blue",
            "minAge",
            "maxAge",
        ]
    ].drop_duplicates(subset="formation")
    for f in uniqueLayers["formation"]:
        thickness[f] = np.mean(
            stratigraphicLayers[stratigraphicLayers["formation"] == f]["thickness"]
        )
    stratigraphicLogData = np.zeros(
        uniqueLayers.shape[0], LoopProjectFile.stratigraphicLayerType
    )
    stratigraphicLogData["layerId"] = range(uniqueLayers.shape[0])
    stratigraphicLogData["layerId"] += 1
    stratigraphicLogData["minAge"] = uniqueLayers["minAge"]
    stratigraphicLogData["maxAge"] = uniqueLayers["maxAge"]
    stratigraphicLogData["name"] = uniqueLayers["formation"]
    stratigraphicLogData["supergroup"] = uniqueLayers["supergroup"]
    stratigraphicLogData["enabled"] = 1
    stratigraphicLogData["rank"] = 0
    stratigraphicLogData["type"] = 4
    stratigraphicLogData["thickness"] = list(thickness.values())

    # Should check format of colour first for "#ffffff" perhaps with regex
    stratigraphicLogData["colour1Red"] = uniqueLayers["colour1Red"]
    stratigraphicLogData["colour1Green"] = uniqueLayers["colour1Green"]
    stratigraphicLogData["colour1Blue"] = uniqueLayers["colour1Blue"]
    stratigraphicLogData["colour2Red"] = [
        int(a * 0.95) for a in stratigraphicLogData["colour1Red"]
    ]
    stratigraphicLogData["colour2Green"] = [
        int(a * 0.95) for a in stratigraphicLogData["colour1Green"]
    ]
    stratigraphicLogData["colour2Blue"] = [
        int(a * 0.95) for a in stratigraphicLogData["colour1Blue"]
    ]
    resp = LoopProjectFile.Set(
        loopFilename, "stratigraphicLog", data=stratigraphicLogData, verbose=True
    )

    if resp["errorFlag"]:
        print(resp["errorString"])

    faults = pd.read_csv(os.path.join(config.output_path, "fault_orientations.csv"))
    faults["formation"] = [re.sub("\.0", "", s) for s in faults["formation"]]
    faultDims = pd.read_csv(os.path.join(config.output_path, "fault_dimensions.csv"))
    faultDims.rename(columns={"Fault": "formation"}, inplace=True)
    faultDims["formation"] = [re.sub("\.0", "", s) for s in faultDims["formation"]]
    faultDisplacements = pd.read_csv(
        os.path.join(config.output_path, "fault_displacements3.csv")
    )
    faultDisplacements.rename(columns={"fname": "formation"}, inplace=True)
    faultDisplacements["formation"] = [
        re.sub("\.0", "", s) for s in faultDisplacements["formation"]
    ]
    faults = faults.merge(faultDims, on="formation")

    minStratAge = np.nanmin(uniqueLayers["minAge"])
    maxStratAge = np.nanmax(uniqueLayers["maxAge"])
    faultObs = pd.read_csv(os.path.join(config.output_path, "faults.csv"))
    faultObs["formation"] = [re.sub("\.0", "", s) for s in faultObs["formation"]]
    faultObs["posOnly"] = 1
    faultsJoined = pd.concat([faults, faultObs])
    if len(faultDims) > 0:
        faultEvents = np.zeros(faultDims.shape[0], LoopProjectFile.faultEventType)
        # The fault eventId is called formation for some reason
        faultEvents["name"] = faultDims["formation"]
        faultEvents["enabled"] = 1
        faultEvents["rank"] = 0
        faultEvents["type"] = 0
        faultEvents["minAge"] = np.linspace(
            minStratAge, maxStratAge, faultDims.shape[0]
        )
        faultEvents["maxAge"] = faultEvents["minAge"]
        avgDisplacements = []
        avgDownthrowDir = []
        for formationName in faultDims["formation"].unique():
            avgDisplacements.append(
                np.average(
                    faultDisplacements[
                        faultDisplacements["formation"] == formationName
                    ]["vertical_displacement"]
                )
            )
            avgDownthrowDir.append(
                np.average(
                    faultDisplacements[
                        faultDisplacements["formation"] == formationName
                    ]["downthrow_dir"]
                )
            )
        faultEvents["avgDisplacement"] = avgDisplacements
        faultEvents["avgDownthrowDir"] = avgDownthrowDir
        faultEvents["influenceDistance"] = faultDims["InfluenceDistance"]
        faultEvents["verticalRadius"] = faultDims["VerticalRadius"]
        faultEvents["horizontalRadius"] = faultDims["HorizontalRadius"]
        faultEvents["colour"] = faultDims["colour"]

        faultEvents["eventId"] = [re.sub(".*_", "", s) for s in faultDims["formation"]]

        resp = LoopProjectFile.Set(
            loopFilename, "faultLog", data=faultEvents, verbose=False
        )
        if resp["errorFlag"]:
            print(resp["errorString"])

        faultsData = np.zeros(
            faultsJoined.shape[0], LoopProjectFile.faultObservationType
        )
        faultsData["eventId"] = [
            re.sub(".*_", "", s) for s in faultsJoined["formation"]
        ]
        faultsData["easting"] = faultsJoined["X"]
        faultsData["northing"] = faultsJoined["Y"]
        faultsData["altitude"] = faultsJoined["Z"]
        faultsData["dipDir"] = faultsJoined["DipDirection"]
        faultsData["dip"] = faultsJoined["dip"]
        faultsData["dipPolarity"] = faultsJoined["DipPolarity"]
        faultsData["displacement"] = 0
        faultsData["posOnly"] = faultsJoined["posOnly"]
        resp = LoopProjectFile.Set(
            loopFilename, "faultObservations", data=faultsData, verbose=True
        )
        if resp["errorFlag"]:
            print(resp["errorString"])

    # each contact contains a location and which formation it is on
    contacts = pd.read_csv(os.path.join(config.output_path, "contacts_clean.csv"))
    layerIds = []
    for form in contacts["formation"]:
        a = bytes(form, "ascii")
        if a in stratigraphicLogData["name"]:
            layerIds.append(
                int(stratigraphicLogData[stratigraphicLogData["name"] == a]["layerId"])
            )
        else:
            layerIds.append(0)
    contactsData = np.zeros(contacts.shape[0], LoopProjectFile.contactObservationType)
    contactsData["layerId"] = layerIds
    contactsData["easting"] = contacts["X"]
    contactsData["northing"] = contacts["Y"]
    contactsData["altitude"] = contacts["Z"]
    # contactsData['dipdir'] = contacts['']
    # contactsData['dip'] = contacts['']
    resp = LoopProjectFile.Set(
        loopFilename, "contacts", data=contactsData, verbose=True
    )
    if resp["errorFlag"]:
        print(resp["errorString"])

    observations = pd.read_csv(
        os.path.join(config.output_path, "orientations_clean.csv")
    )
    layerIds = []
    for form in observations["formation"]:
        a = bytes(form, "ascii")
        if a in stratigraphicLogData["name"]:
            layerIds.append(
                int(stratigraphicLogData[stratigraphicLogData["name"] == a]["layerId"])
            )
        else:
            layerIds.append(0)
    observations["layer"] = "s0"
    observationsData = np.zeros(
        observations.shape[0], LoopProjectFile.stratigraphicObservationType
    )
    observationsData["layerId"] = layerIds
    observationsData["easting"] = observations["X"]
    observationsData["northing"] = observations["Y"]
    observationsData["altitude"] = observations["Z"]
    observationsData["dipDir"] = observations["azimuth"]
    observationsData["dip"] = observations["dip"]
    observationsData["dipPolarity"] = observations["polarity"]
    observationsData["layer"] = observations["layer"]
    resp = LoopProjectFile.Set(
        loopFilename, "stratigraphicObservations", data=observationsData, verbose=True
    )
    if resp["errorFlag"]:
        print(resp["errorString"])

    # Check created file is valid
    if LoopProjectFile.CheckFileValid(loopFilename):
        return loopFilename
    else:
        return None


##########################################################################
# Import outputs from map2loop to gempy and view with pyvtk
# loop2gempy(test_data_name,tmp_path,vtk_pth,orientations_file,contacts_file,groups_file,dtm_reproj_file,bbox,model_base, model_top,vtk)
# Args:
# test_data_name root name of project
# tmp_path path of temp files directory
# vtk_pth path of vtk output directory
# orientations_file path of orientations file
# contacts_file path of contacts file
# groups_file path of groups file
# dtm_reproj_file path of dtm file
# bbox model bounding box
# model_base z value ofbase of model
# model_top z value of top of model
# vtk flag as to wether to save out model to vtk
#
# Calculates model and displays in external vtk viewer
##########################################################################
def loop2gempy(*args, **kwargs):
    """Calculate the model using gempy as backend.
    At the moment there is not support for finite faults since gempy does not
     accept passing the ellipsoid parameters directly.
    :param contacts_file (str): path of contacts file
    :param orientations_file: path of orientations file
    :param bbox: model bounding box
    :param groups_file: path of groups file
    :param model_base: z value ofbase of model
    :param model_top: z value of top of model
    :param dtm_reproj_file: path of dtm file
    :param faults_contact: path of contacts file with fault data
    :param faults_orientations: path of orientations file with fault data
    :param faults_rel_matrix: bool matrix describing the interaction between groups. Rows offset columns
    :param faults_groups_rel: bool matrix describing the interaction between faults and features
    :param faults_faults_rel: bool matrix describing the interaction between faults and faults
    :param model_name: name of the model
    :param compute (bool): Default True. Whether or not compute the model
    :param vtk (bool): Default False. Whether or not visualize the model
    :param vtk_path (str): Default None. Path of vtk output directory
    :param plot_3d_kwargs (dict): kwargs for `gempy.plot_3d`
    :return: gempy.Project
    """
    from gempy.addons.map2gempy import loop2gempy

    geo_model = loop2gempy(*args, **kwargs)
    return geo_model


def ElementToPandas(loopFilename, element, loopCompoundType):
    """
    **ElementToPandas** - Exports one element of the loop project file
    to a pandas dataframe

    Parameters
    ----------
    loopFilename: string
        The filename of the loop project file
    element: string
        The name of the element to extract
    loopCompoundType: numpy.compoundType
        The numpy data structure that the element is stored in

    Returns
        pandas dataframe
    -------

    """
    resp = LoopProjectFile.Get(loopFilename, element)
    if resp["errorFlag"]:
        print(resp["errorString"])
        df = pd.DataFrame()
    else:
        columns = list(loopCompoundType.names)
        df = pd.DataFrame.from_records(resp["value"], columns=columns)
        df = df.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
        df.set_index(columns[0], inplace=True)
    return df
