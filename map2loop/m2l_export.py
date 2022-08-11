import numpy as np
import pandas as pd
import os
import LoopProjectFile
import rasterio
import rasterio.transform
import rasterio.plot
import matplotlib
import beartype
from .config import Config
import re


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
        except Exception:
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
    transform = rasterio.transform.from_origin(dtm.bounds[0], dtm.bounds[3], scale, scale)

    memfile = rasterio.MemoryFile()
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
            if existingExtents is None:
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
    thickness = {}
    sum = 0
    num = 0
    for f in form2supergroup["formation"]:
        # Ignore cover or cover_up layers
        if f == "cover" or f == "cover_up":
            pass
        elif len(stratigraphicLayers[stratigraphicLayers["formation"] == f]):
            thickness[f] = np.median(
                stratigraphicLayers[stratigraphicLayers["formation"] == f]["thickness"]
            )
            sum += thickness[f]
            num += 1
        else:
            thickness[f] = -1
    # Note: thickness of unit without any thickness values is the mean thickness of those with values
    mean_thickness = sum / num
    for f in form2supergroup["formation"]:
        if f == "cover" or f == "cover_up":
            pass
        elif thickness[f] == -1:
            thickness[f] = mean_thickness

    stratAges = pd.read_csv(os.path.join(config.tmp_path, "age_sorted_groups.csv"))[
        ["group_", "min", "max"]
    ]
    stratAges.rename(
        columns={"group_": "supergroup", "min": "minAge", "max": "maxAge"}, inplace=True
    )
    stratLayers = pd.merge(form2supergroup, stratAges, on=["supergroup"])
    stratLayers["colour1Red"] = [int(a[1:3], 16) for a in stratLayers["colour"]]
    stratLayers["colour1Green"] = [int(a[3:5], 16) for a in stratLayers["colour"]]
    stratLayers["colour1Blue"] = [int(a[5:7], 16) for a in stratLayers["colour"]]
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
    faults["formation"] = [re.sub("\\.0", "", s) for s in faults["formation"]]
    faultDims = pd.read_csv(os.path.join(config.output_path, "fault_dimensions.csv"))
    faultDims.rename(columns={"Fault": "formation"}, inplace=True)
    faultDims["formation"] = [re.sub("\\.0", "", s) for s in faultDims["formation"]]
    faultDisplacements = pd.read_csv(
        os.path.join(config.output_path, "fault_displacements3.csv")
    )
    faultDisplacements.rename(columns={"fname": "formation"}, inplace=True)
    faultDisplacements["formation"] = [
        re.sub("\\.0", "", s) for s in faultDisplacements["formation"]
    ]
    faults = faults.merge(faultDims, on="formation")

    minStratAge = np.nanmin(uniqueLayers["minAge"])
    maxStratAge = np.nanmax(uniqueLayers["maxAge"])
    faultObs = pd.read_csv(os.path.join(config.output_path, "faults.csv"))
    faultObs["formation"] = [re.sub("\\.0", "", s) for s in faultObs["formation"]]
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
