import os
import geopandas as gpd
from shapely.geometry import LineString, Polygon, MultiLineString
import warnings
import numpy as np
import pandas as pd
from math import sqrt

# from .mapdata import MapData # Creates a circular dependency
from .m2l_enums import Datatype, Datastate, VerboseLevel
import beartype
from .config import Config


# explodes polylines and modifies objectid for exploded parts
def explode_polylines(indf, c_l, dst_crs):
    outdf = gpd.GeoDataFrame(columns=indf.columns, crs=dst_crs)
    for idx, row in indf.iterrows():
        if type(row.geometry) == LineString:
            rowdf = gpd.GeoDataFrame(data=dict(row), index=[0],crs=dst_crs)
            outdf = pd.concat([outdf, rowdf], ignore_index=True)
        if type(row.geometry) == MultiLineString:
            multdf = gpd.GeoDataFrame(columns=indf.columns, crs=dst_crs)
            recs = len(row.geometry.geoms)
            row_dict = dict(row)
            row_dict.pop("geometry")
            rowdf = gpd.GeoDataFrame(data=row_dict, index=[0])
            multdf = pd.concat([multdf, *[rowdf] * recs], ignore_index=True)
            i = 0
            for geom in range(recs):
                multdf.loc[geom, "geometry"] = row.geometry.geoms[geom]
                multdf.loc[geom, c_l["o"]] = (
                    str(multdf.loc[geom, c_l["o"]]) + "_" + str(i)
                )
                print(
                    "map2loop warning: Fault_" + multdf.loc[geom, c_l["o"]],
                    "is one of a set of duplicates, so renumbering",
                )
                i = i + 1
            multdf.set_crs(crs=dst_crs, inplace=True)
            outdf = pd.concat([outdf, multdf], ignore_index=True)
    return outdf


@beartype.beartype
def check_all_maps(
    mapdata,
    config: Config,
    use_roi_clip,
    roi_clip_path,
    verbose_level: VerboseLevel = VerboseLevel.ALL,
):

    m2l_errors = []
    m2l_warnings = []
    if (
        config.bbox[3] < config.bbox[1]
        or config.bbox[2] < config.bbox[0]
        or config.bbox_3d["top"] < config.bbox_3d["base"]
    ):
        m2l_errors.append("bounding box has negative range for x or y or z")

    f = open(os.path.join(config.tmp_path, "bbox.csv"), "w")
    f.write("minx,miny,maxx,maxy,lower,upper\n")
    ostr = "{},{},{},{},{},{}\n".format(
        config.bbox[0],
        config.bbox[1],
        config.bbox[2],
        config.bbox[3],
        config.bbox_3d["base"],
        config.bbox_3d["top"],
    )
    f.write(ostr)
    f.close()

    if use_roi_clip:
        polygo = gpd.read_file(roi_clip_path)
    else:
        y_point_list = [
            config.bbox[1],
            config.bbox[1],
            config.bbox[3],
            config.bbox[3],
            config.bbox[1],
        ]
        x_point_list = [
            config.bbox[0],
            config.bbox[2],
            config.bbox[2],
            config.bbox[0],
            config.bbox[0],
        ]
        bbox_geom = Polygon(zip(x_point_list, y_point_list))
        polygo = gpd.GeoDataFrame(
            index=[0], crs=mapdata.working_projection, geometry=[bbox_geom]
        )

    for i in [
        Datatype.GEOLOGY,
        Datatype.STRUCTURE,
        Datatype.FAULT,
        Datatype.FOLD,
        Datatype.MINERAL_DEPOSIT,
    ]:
        if not mapdata.get_filename(i).startswith("http") and not os.path.isfile(
            mapdata.get_filename(i)
        ):
            m2l_errors.append("file " + mapdata.get_filename(i) + " not found")

    orientations = check_structure_map(
        mapdata.get_map_data(Datatype.STRUCTURE),
        config.c_l,
        m2l_warnings,
        m2l_errors,
        verbose_level,
    )
    geology = check_geology_map(
        mapdata.get_map_data(Datatype.GEOLOGY),
        config.c_l,
        config.run_flags["ignore_codes"],
        m2l_warnings,
        m2l_errors,
        verbose_level=verbose_level,
    )
    folds = check_fold_map(
        mapdata.get_map_data(Datatype.FOLD),
        config.c_l,
        m2l_warnings,
        m2l_errors,
        verbose_level,
    )
    faults = check_fault_map(
        mapdata.get_map_data(Datatype.FAULT),
        config.c_l,
        m2l_warnings,
        m2l_errors,
        verbose_level,
    )
    mindeps = check_mindep_map(
        mapdata.get_map_data(Datatype.MINERAL_DEPOSIT),
        config.c_l,
        m2l_warnings,
        m2l_errors,
        verbose_level,
    )

    if len(m2l_warnings) > 0:
        print("\nWarnings:")
        warnings.warn("The warnings listed above were issued")
        for w in m2l_warnings:
            print("    ", w)
    if len(m2l_errors) > 0:
        print("\nErrors:")
        warnings.warn(
            "The errors listed above must be fixed prior to rerunning map2loop"
        )
        for e in m2l_errors:
            print("    ", e)
        sep = "\12"
        raise NameError(
            sep.join(m2l_errors) + "\n map2loop error: Fix errors before running again"
        )

    if len(m2l_errors) == 0:
        mapdata.data[Datatype.STRUCTURE] = orientations
        mapdata.data_states[Datatype.STRUCTURE] = Datastate.CONVERTED
        mapdata.data[Datatype.GEOLOGY] = geology
        mapdata.data_states[Datatype.GEOLOGY] = Datastate.CONVERTED
        mapdata.data[Datatype.FOLD] = folds
        mapdata.data_states[Datatype.FOLD] = Datastate.CONVERTED
        mapdata.data[Datatype.FAULT] = faults
        mapdata.data_states[Datatype.FAULT] = Datastate.CONVERTED
        mapdata.data[Datatype.MINERAL_DEPOSIT] = mindeps
        mapdata.data_states[Datatype.MINERAL_DEPOSIT] = Datastate.CONVERTED
        output_modified_maps(
            orientations,
            geology,
            folds,
            faults,
            mindeps,
            polygo,
            config.c_l,
            mapdata.working_projection,
            config.tmp_path,
        )


def rename_columns(
    dataframe,
    original_column_name,
    replacement_column_name,
    m2l_errors,
    verbose_level=VerboseLevel.ALL,
):
    if dataframe is not None and original_column_name != replacement_column_name:
        if original_column_name in dataframe.columns:
            if replacement_column_name in dataframe.columns:
                dataframe.rename(
                    columns={
                        replacement_column_name: "old_" + replacement_column_name,
                        original_column_name: replacement_column_name,
                    },
                    inplace=True,
                )
                if verbose_level != VerboseLevel.NONE:
                    print(
                        f'Moving column {replacement_column_name} to "old_{replacement_column_name}"'
                    )
            else:
                dataframe.rename(
                    columns={original_column_name: replacement_column_name},
                    inplace=True,
                )
        else:
            m2l_errors.append(
                f"ERROR: data does not contain column {original_column_name}"
            )
    return dataframe


def check_structure_map(
    orientations, c_l, m2l_warnings, m2l_errors, verbose_level=VerboseLevel.ALL
):
    # Check orientation points valid
    if orientations is None or type(orientations) != gpd.GeoDataFrame:
        m2l_errors.append("No structure map found")
        return None

    # Error if insufficient sturucture points
    if len(orientations) < 2:
        m2l_errors.append(
            "not enough orientations to complete calculations (need at least 2), projection may be inconsistent"
        )

    # Check for strike and convert to dip direction
    if c_l["otype"] == "strike" and c_l["dd"] in orientations.columns:
        if verbose_level != VerboseLevel.NONE:
            print("converting strike/dip orientation to dipdir/dip")
        orientations[c_l["dd"]] = orientations.apply(
            lambda row: row[c_l["dd"]] + 90.0, axis=1
        )

    # Check for duplicate column references
    if c_l["sf"] == c_l["ds"]:
        if c_l["sf"] in orientations.columns:
            orientations = rename_columns(
                orientations, c_l["sf"], "STRUCTURE_TYPE", m2l_errors, verbose_level
            )
            orientations[c_l["ds"]] = orientations["STRUCTURE_TYPE"].copy()
            c_l["sf"] = "STRUCTURE_TYPE"
    if c_l["bo"] == c_l["ds"]:
        if c_l["bo"] in orientations.columns:
            orientations = rename_columns(
                orientations, c_l["bo"], "POLARITY", m2l_errors, verbose_level
            )
            orientations[c_l["ds"]] = orientations["POLARITY"].copy()
            c_l["bo"] = "POLARITY"

    # Fill in missing data with default values
    if c_l["sf"] not in orientations.columns:
        orientations[c_l["sf"]] = "Bed"
        m2l_warnings.append(
            f"field named '{str(c_l['sf'])}' added with default value \"Bed\""
        )

    if c_l["ds"] not in orientations.columns:
        orientations[c_l["ds"]] = ""
        m2l_warnings.append(
            f"field named '{str(c_l['ds'])}' added with default empty value"
        )

    if c_l["gi"] not in orientations.columns:
        orientations[c_l["gi"]] = np.arange(len(orientations))
        m2l_warnings.append(f"field named '{str(c_l['gi'])}' added with default value")

    # Convert whitespace entries to nans
    orientations = orientations.replace(r"^\s+$", np.nan, regex=True)

    # Warn about number of nans present in data
    for code in ("sf", "d", "dd", "gi"):
        if c_l[code] in orientations.columns:
            nans = orientations[c_l[code]].isnull().sum()
            if nans > 0:
                m2l_warnings.append(
                    f"{str(nans)} NaN/blank found in column '{str(c_l[code])}' of orientations file, replacing with 0"
                )
                orientations[c_l[code]].fillna(0, inplace=True)

    # Remove orientation entries with unknown dip
    orientations = orientations[orientations[c_l["d"]] != -999]

    unique_o = set(orientations[c_l["gi"]])
    if not len(unique_o) == len(orientations):
        m2l_warnings.append("duplicate orientation point unique IDs")

    # Ensuring type of column is correct
    orientations[c_l["dd"]] = pd.to_numeric(orientations[c_l["dd"]])
    orientations[c_l["d"]] = pd.to_numeric(orientations[c_l["d"]])

    # convert c_l codes to meaningful column names
    orientations = rename_columns(
        orientations, c_l["d"], "DIP", m2l_errors, verbose_level
    )
    orientations = rename_columns(
        orientations, c_l["dd"], "DIPDIR", m2l_errors, verbose_level
    )
    orientations = rename_columns(
        orientations, c_l["sf"], "STRUCTURE_TYPE", m2l_errors, verbose_level
    )
    orientations = rename_columns(
        orientations, c_l["ds"], "DESCRIPTION", m2l_errors, verbose_level
    )
    orientations = rename_columns(
        orientations, c_l["bo"], "POLARITY", m2l_errors, verbose_level
    )
    orientations = rename_columns(
        orientations, c_l["gi"], "STRUCTURE_POINT_ID", m2l_errors, verbose_level
    )

    # Reset c_l labels to new labels (temporary unit all c_l references removed from non-mapchecker code)
    # c_l['d'] = 'DIP'
    # c_l['dd'] = 'DIPDIR'
    # c_l['sf'] = 'STRUCTURE_TYPE'
    # c_l['ds'] = 'DESCRIPTION'
    # c_l['bo'] = 'POLARITY'
    # c_l['go'] = 'STRUCTURE_POINT_ID'

    if verbose_level != VerboseLevel.NONE:
        show_metadata(orientations, "orientations layer")

    return orientations


def _explode_intrusives(geology, c_l):
    """Break multipart geometries of intrusions into individual features

    Parameters
    ----------
    geology : _type_
        _description_
    """
    geol_clip_tmp = geology.copy(deep=False)
    # print('raw',len(geol_clip_tmp))
    geol_clip_tmp.crs = geology.crs
    geol_clip_tmp = geol_clip_tmp.dissolve(by=c_l["c"], aggfunc="first")
    # Note: Can't use "INTRUSIVE" and "SILL" in place of c_l["..."] as they are str field
    # to use for look up not column names
    geol_clip_tmp = geol_clip_tmp[
        geol_clip_tmp[c_l["r1"]].str.contains(c_l["intrusive"])
        & ~geol_clip_tmp[c_l["ds"]].str.contains(c_l["sill"])
    ]
    # print('tmp',len(geol_clip_tmp))
    geol_clip_tmp.reset_index(inplace=True)
    geol_clip_tmp = geol_clip_tmp.explode(index_parts=False, ignore_index=True)
    if "level_0" in geol_clip_tmp.columns:
        geol_clip_tmp.drop(labels="level_0", axis=1, inplace=True)

    geol_clip_tmp.reset_index(inplace=True)
    geol_clip_not = geology[
        (~geology[c_l["r1"]].str.contains(c_l["intrusive"]))
        | geology[c_l["ds"]].str.contains(c_l["sill"])
    ]
    # print('not',len(geol_clip_not))
    geol_clip_not.index = pd.RangeIndex(
        start=10000, stop=10000 + len(geol_clip_not), step=1
    )
    if "level_0" in geol_clip_not.columns:
        geol_clip_not.drop(labels="level_0", axis=1, inplace=True)
    geol_clip_not.reset_index(inplace=True)
    geol_clip_tmp[c_l["c"]] = (
        geol_clip_tmp[c_l["c"]].astype(str) + "_" + geol_clip_tmp.index.astype(str)
    )
    geol_clip_tmp[c_l["o"]] = geol_clip_tmp.index
    geol_clip_tmp[c_l["g"]] = geol_clip_tmp[c_l["c"]]

    geology = pd.concat([geol_clip_tmp, geol_clip_not], sort=False)
    # print('geol',len(geology))
    # #TODO Lg- not sure what level_0/level_1 are for...
    if "level_0" in geology.columns:
        geology.drop(labels="level_0", axis=1, inplace=True)
    if "level_1" in geology.columns:
        geology.drop(labels="level_1", axis=1, inplace=True)
    return geology


def check_geology_map(
    geology,
    c_l={},
    ignore_codes=[],
    m2l_warnings=[],
    m2l_errors=[],
    explode_intrusives=True,
    verbose_level=VerboseLevel.ALL,
) -> gpd.GeoDataFrame:
    """Checks whether the geological map is valid or not, appends any errors to a string.

    Parameters
    ----------
    geology : GeoDataFrame
        _description_
    c_l : dict, optional
        dictionary for converting fields into m2l requirements, by default {}
    ignore_codes : list, optional
        Code or prefix of a code that are skipped by m2l, e.g. cover , by default []
    m2l_warnings : list, optional
        _description_, by default []
    m2l_error : list, optional
        _description_, by default []
    explode_intrusives : bool, optional
        whether to break each of the intrusions into separate features, by default True
    verbose_level : VerboseLevel, optional
        change level of m2l output, by default VerboseLevel.ALL

    Returns
    -------
    gpd.GeoDataFrame
        geology geodataframe with the corrections applied

    Raises
    ------
    ValueError
        Error when the geology file is none or not a geodataframe
    ValueError
        Error when the geology file is empty
    """
    # Check geology map
    if geology is None or type(geology) != gpd.GeoDataFrame:
        m2l_errors.append("No geology map found")
        return geology

    if len(geology) == 0:
        m2l_errors.append("Empty geology map found")
        return geology

    # Need to break these up into column renaming fields and other types e.g. (bedding, intrusive, sill)
    # Do we create a new version of c_l for these search fields???
    # c_l.search_strings = {'volcanic':'VOLCANIC','intrusive':'INTRUSIVE','sill':'SILL'}
    # temporary map between old c_l and new_mapping
    # new_mapping = {'c':'UNIT_NAME','u':'CODE','r1':'ROCKTYPE1','r2':'ROCKTYPE2','ds':'DESCRIPTION','r2':'ROCKTYPE2',
    # 'min':'MIN_AGE','max':'MAX_AGE','g':'GROUP','o':'GEOLOGY_FEATURE_ID','g2':'GROUP2'}
    # geology = geology.rename(columns=dict(zip(c_l.values(),[ new_mapping[item] if item in new_mapping else item for item in c_l.keys() ])))
    # Remapping done at the end of function

    # make each pluton its own formation and group
    if c_l["o"] not in geology.columns:  # object id
        geology = geology.reset_index()
        geology[c_l["o"]] = geology.index

    # make each pluton its own formation and group
    if c_l["r1"] not in geology.columns:
        m2l_warnings.append("No extra litho for geology polygons")
        geology[c_l["r1"]] = np.nan
    geology[c_l["r1"]].fillna("not known", inplace=True)

    if c_l["ds"] not in geology.columns:
        m2l_warnings.append("No extra litho for geology polygons")
        geology[c_l["ds"]] = np.nan
    geology[c_l["ds"]].fillna("not known", inplace=True)

    if explode_intrusives:
        geology = _explode_intrusives(geology, c_l)

    unique_c = list(set(geology[c_l["c"]]))
    if len(unique_c) < 2:
        m2l_errors.append(
            "The selected area only has one formation in it, so no model can be calculated"
        )

    unique_g = set(geology[c_l["o"]])
    if not len(unique_g) == len(geology):
        m2l_warnings.append("duplicate geology polygon unique IDs")

    # Check for nans
    nans = geology[c_l["c"]].isnull().sum()
    if nans > 0:
        m2l_errors.append(
            ""
            + str(nans)
            + ' NaN/blank found in column "'
            + str(c_l["c"])
            + '" of geology file, please fix'
        )

    if c_l["g"] == "No_col" or not c_l["g"] in geology.columns:
        m2l_warnings.append("No secondary strat coding for geology polygons")
        geology[c_l["g"]] = "Top"

    # Ensure that group column has some valid information
    geology = geology.replace(r"^\s+$", np.nan, regex=True)
    geology = geology.replace(",", " ", regex=True)
    geology[c_l["g"]].fillna(geology[c_l["g2"]], inplace=True)
    geology[c_l["g"]].fillna(geology[c_l["c"]], inplace=True)

    # Continue filling in missing columns
    if c_l["r2"] == "No_col" or not c_l["r2"] in geology.columns:
        m2l_warnings.append("No more extra litho for geology polygons")
        geology[c_l["r2"]] = "Nope"

    if c_l["min"] not in geology.columns:
        m2l_warnings.append("No min age for geology polygons")
        geology[c_l["min"]] = 0

    if c_l["max"] not in geology.columns:
        m2l_warnings.append("No max age for geology polygons")
        geology[c_l["max"]] = 100

    # Convert types
    geology[c_l["max"]] = geology[c_l["max"]].astype(np.float64)
    geology[c_l["min"]] = geology[c_l["min"]].astype(np.float64)

    if c_l["c"] not in geology.columns:
        m2l_errors.append("Must have primary strat coding field for geology polygons")

    for code in ("c", "g", "g2", "ds", "u", "r1"):
        if c_l[code] in geology.columns:
            geology[c_l[code]].str.replace(",", " ")
            if code == "c" or code == "g" or code == "g2":
                geology[c_l[code]] = geology[c_l[code]].str.replace(
                    "[ -/?]", "_", regex=True
                )

            nans = geology[c_l[code]].isnull().sum()
            if nans > 0:
                m2l_warnings.append(
                    ""
                    + str(nans)
                    + ' NaN/blank found in column "'
                    + str(c_l[code])
                    + '" of geology file, replacing with 0'
                )
                geology[c_l[code]].fillna("0", inplace=True)

    # Remove geology that has code starting with ignore_codes strings
    for code in ignore_codes:
        geology = geology[~geology[c_l["u"]].str.startswith(code)]

    if verbose_level != VerboseLevel.NONE:
        show_metadata(geology, "geology layer")

    # convert c_l codes to meaningful column names
    geology = rename_columns(geology, c_l["g"], "GROUP", m2l_errors, verbose_level)
    geology = rename_columns(geology, c_l["u"], "CODE", m2l_errors, verbose_level)
    geology = rename_columns(geology, c_l["g2"], "GROUP2", m2l_errors, verbose_level)
    geology = rename_columns(
        geology, c_l["ds"], "DESCRIPTION", m2l_errors, verbose_level
    )
    geology = rename_columns(geology, c_l["c"], "UNIT_NAME", m2l_errors, verbose_level)
    geology = rename_columns(geology, c_l["r1"], "ROCKTYPE1", m2l_errors, verbose_level)
    geology = rename_columns(geology, c_l["r2"], "ROCKTYPE2", m2l_errors, verbose_level)
    geology = rename_columns(geology, c_l["min"], "MIN_AGE", m2l_errors, verbose_level)
    geology = rename_columns(geology, c_l["max"], "MAX_AGE", m2l_errors, verbose_level)
    geology = rename_columns(
        geology, c_l["o"], "GEOMETRY_OBJECT_ID", m2l_errors, verbose_level
    )

    return geology


def check_fold_map(
    folds, c_l, m2l_warnings, m2l_errors, verbose_level=VerboseLevel.ALL
):
    # Process fold polylines
    if folds is not None:
        if len(folds) > 0:
            if not c_l["o"] in folds.columns:
                folds = folds.reset_index()
                folds[c_l["o"]] = folds.index
            unique_g = set(folds[c_l["o"]])

            if not len(unique_g) == len(folds):
                m2l_warnings.append("duplicate fold polyline unique IDs")

            folds = folds[folds[c_l["ff"]].str.contains(c_l["fold"], case=False)]

            folds = folds.replace(r"^\s+$", np.nan, regex=True)

            for code in ("ff", "t"):
                if c_l["ff"] == "No_col" or not c_l["ff"] in folds.columns:
                    m2l_warnings.append("No fold code for fold polylines")
                    c_l["ff"] = "ff"
                    folds[c_l["ff"]] = c_l["fold"]

                if c_l["t"] == "No_col" or not c_l["t"] in folds.columns:
                    m2l_warnings.append("No fold polarity for fold polylines")
                    c_l["t"] = "t"
                    folds[c_l["t"]] = "None"

                if c_l[code] in folds.columns:
                    folds[c_l[code]].str.replace(",", " ")

                    nans = folds[c_l[code]].isnull().sum()
                    if nans > 0:
                        m2l_warnings.append(
                            ""
                            + str(nans)
                            + ' NaN/blank found in column "'
                            + str(c_l[code])
                            + '" of folds file, replacing with 0'
                        )
                        folds[c_l[code]].fillna("0", inplace=True)

            if len(folds) > 0:
                folds_explode = explode_polylines(folds, c_l, folds.crs)
                if len(folds_explode) > len(folds):
                    m2l_warnings.append(
                        "some folds are MultiPolyLines, and have been split"
                    )
                # folds_explode.crs = dst_crs
                folds = folds_explode

            if verbose_level != VerboseLevel.NONE:
                show_metadata(folds, "fold layer")
        else:
            print("No folds in area, projection may be inconsistent")
    return folds


def check_fault_map(
    faults, c_l, m2l_warnings, m2l_errors, verbose_level=VerboseLevel.ALL
):
    # Check fault map
    if faults is None or type(faults) != gpd.GeoDataFrame:
        m2l_errors.append("No fault data found")
        return faults

    # Process fault polylines
    faults = faults[faults[c_l["f"]].str.contains(c_l["fault"], case=False)]
    faults = faults.replace(r"^\s+$", np.nan, regex=True)

    # Check for missing columns and fill with default values
    if not c_l["o"] in faults.columns:
        m2l_warnings.append(
            'field named "' + str(c_l["o"]) + '" added with default value'
        )
        faults[c_l["o"]] = np.arange(len(faults))

    if c_l["f"] == "No_col" or not c_l["f"] in faults.columns:
        m2l_warnings.append("No fault type for fault polylines")
        c_l["f"] = "ftype"
        faults[c_l["f"]] = c_l["fault"]

    if c_l["fdip"] == "No_col" or not c_l["fdip"] in faults.columns:
        m2l_warnings.append("No fault dip for fault polylines")
        c_l["fdip"] = "fdip"
        faults[c_l["fdip"]] = c_l["fdipnull"]

    if c_l["fdipdir"] == "No_col" or not c_l["fdipdir"] in faults.columns:
        m2l_warnings.append("No fault dip direction for fault polylines")
        c_l["fdipdir"] = "fdipdir"
        faults[c_l["fdipdir"]] = -999

    if c_l["fdipest"] == "No_col" or not c_l["fdipest"] in faults.columns:
        m2l_warnings.append("No fault dip estimate for fault polylines")
        c_l["fdipest"] = "fdipest"
        faults[c_l["fdipest"]] = "None"

    if c_l["fdipest_vals"] == "No_col" or not c_l["fdipest_vals"] in faults.columns:
        m2l_warnings.append("No fault dip estimate text for fault polylines")
        c_l["fdipest_vals"] = "fdipest_vals"
        faults[c_l["fdipest_vals"]] = "None"

    if c_l["n"] == "No_col" or not c_l["n"] in faults.columns:
        m2l_warnings.append("No fault name for fault polylines")
        c_l["n"] = "fname"
        faults[c_l["n"]] = "None"

    # Check for nans and warn
    for code in ("f", "o", "fdip", "fdipdir", "fdipest"):
        if c_l[code] in faults.columns:
            nans = faults[c_l[code]].isnull().sum()
            if nans > 0:
                m2l_warnings.append(
                    ""
                    + str(nans)
                    + ' NaN/blank found in column "'
                    + str(c_l[code])
                    + '" of fault file, replacing with -999'
                )
                faults[c_l[code]].fillna("-999", inplace=True)

    #
    faults[c_l["o"]] = faults[c_l["o"]].astype(int)

    # Check for duplicates and warn
    unique_f = set(faults[c_l["o"]])
    if not len(unique_f) == len(faults):
        m2l_errors.append("duplicate fault/fold polyline unique IDs")

    if len(faults) > 0:
        faults_explode = explode_polylines(faults, c_l, faults.crs)
        if len(faults_explode) > len(faults):
            m2l_warnings.append("some faults are MultiPolyLines, and have been split")
        # faults_explode.to_crs = dst_crs
        faults = faults_explode

        lengths = []

        for indx, flt in faults.iterrows():
            if c_l["fault"].lower() in flt[c_l["f"]].lower():
                if flt.geometry.type == "LineString":
                    flt_ls = LineString(flt.geometry)
                    if len(flt_ls.coords) > 0:
                        dlsx = (
                            flt_ls.coords[0][0]
                            - flt_ls.coords[len(flt_ls.coords) - 1][0]
                        )
                        dlsy = (
                            flt_ls.coords[0][1]
                            - flt_ls.coords[len(flt_ls.coords) - 1][1]
                        )
                        strike = sqrt((dlsx * dlsx) + (dlsy * dlsy))
                        lengths.append(strike)
                    else:
                        lengths.append(0)
                    # print(flt)

        faults["f_length"] = lengths
        faults = faults[faults["f_length"] > 500].copy()

        if verbose_level != VerboseLevel.NONE:
            show_metadata(faults, "fault layer")
    else:
        m2l_warnings.append("No faults in area, projection may be inconsistent")

    if c_l["o"] == c_l["n"]:
        if c_l["n"] in faults.columns:
            faults = rename_columns(faults, c_l["n"], "NAME", m2l_errors, verbose_level)
            faults[c_l["o"]] = faults["NAME"].copy()
            c_l["n"] = "NAME"

    # convert c_l codes to meaningful column names
    faults = rename_columns(faults, c_l["f"], "FEATURE", m2l_errors, verbose_level)
    faults = rename_columns(
        faults, c_l["o"], "GEOMETRY_OBJECT_ID", m2l_errors, verbose_level
    )
    faults = rename_columns(faults, c_l["n"], "NAME", m2l_errors, verbose_level)
    faults = rename_columns(faults, c_l["fdip"], "DIP", m2l_errors, verbose_level)
    faults = rename_columns(faults, c_l["fdipdir"], "DIPDIR", m2l_errors, verbose_level)
    faults = rename_columns(
        faults, c_l["fdipest"], "DIP_ESTIMATE", m2l_errors, verbose_level
    )
    return faults


def check_mindep_map(
    mindeps, c_l, m2l_warnings, m2l_errors, verbose_level=VerboseLevel.ALL
):
    # Process mindep points
    if mindeps is not None:
        try:
            if len(mindeps) == 0:
                m2l_warnings.append("no mindeps for analysis")
            else:
                mindeps = mindeps.replace(r"^\s+$", np.nan, regex=True)

                for code in ("msc", "msn", "mst", "mtc", "mscm", "mcom"):
                    if c_l[code] == "No_col":
                        mindeps[c_l[code]] = "No_col"
                    if not c_l[code] in mindeps.columns:
                        m2l_errors.append(
                            'field named "'
                            + str(c_l[code])
                            + '" not found in mineral deposits file'
                        )
                    else:
                        nans = mindeps[c_l[code]].isnull().sum()
                        if nans > 0:
                            m2l_warnings.append(
                                str(nans)
                                + " NaN/blank found in column "
                                + str(c_l[code])
                                + " of mindep file, replacing with 0"
                            )
                            mindeps[c_l[code]].fillna("0", inplace=True)
            if verbose_level != VerboseLevel.NONE:
                show_metadata(mindeps, "mindeps layer")

        except Exception:
            m2l_warnings.append("no mindeps for analysis")
            mindeps = 0

        else:
            mindeps = None
    return mindeps


def output_modified_maps(
    orientations, geology, folds, faults, mindeps, polygo, c_l, dst_crs, tmp_path
):
    structure_filename = ""
    geol_filename = ""
    fault_filename = ""
    mindep_filename = ""
    fold_filename = ""
    # explode fault/fold multipolylines
    # sometimes faults go off map and come back in again which after clipping creates multipolylines
    if folds is not None:
        fold_filename = os.path.join(tmp_path, "folds_clip.shp")
        if len(folds) > 0:
            folds = folds.dropna(subset=["geometry"])
            folds.to_file(fold_filename)
        else:
            print("\nFold layer metadata\n--------------------")
            print("No folds found, projection may be inconsistent")

    if faults is not None:
        fault_filename = os.path.join(tmp_path, "faults_clip.shp")
        if len(faults) > 0:
            faults.crs = dst_crs
            faults.to_file(fault_filename)
        else:
            print("\nFault layer metadata\n--------------------")
            print("No faults found, projection may be inconsistent")

    if geology is not None:
        geol_clip = gpd.overlay(geology, polygo, how="intersection")
        if len(geol_clip) > 0:
            geol_clip.crs = dst_crs
            geol_filename = os.path.join(tmp_path, "geol_clip.shp")
            geol_clip.to_file(geol_filename)

    if orientations is not None:
        if len(orientations) > 0:
            structure_filename = os.path.join(tmp_path, "structure_clip.shp")
            orientations.crs = dst_crs
            orientations.to_file(structure_filename)

    if mindeps is not None:
        try:
            if len(mindeps) > 0:
                mindep_filename = os.path.join(tmp_path, "mindeps_clip.shp")
                mindeps.crs = dst_crs
                mindeps.to_file(mindep_filename)
        except Exception:
            print("Not clipping mindeps, dataframe is empty")
    print("\nNo errors found, clipped and updated files saved to tmp")

    return (
        structure_filename,
        geol_filename,
        fault_filename,
        mindep_filename,
        fold_filename,
        c_l,
    )


def show_metadata(gdf, name):
    if len(gdf) > 0:
        print("\n", name, " metadata\n--------------------")
        print("    bbox", gdf.total_bounds)
        print("    CRS", gdf.crs)
        print("    # items", len(gdf))
        types = []
        for i, g in gdf.iterrows():
            if g.geometry.type not in types:
                types.append(g.geometry.type)

        print("    Data types", types)
    else:
        print("\n", name, " metadata\n--------------------")
        print("    empty file, check contents")
