import numpy as np
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from math import (
    atan2,
    asin,
    radians,
    degrees,
    sqrt,
    pow,
    acos,
    fabs,
    tan,
    isnan,
)
import geopandas as gpd
from geopandas import GeoDataFrame
import pandas as pd
import os
from shapely.geometry import LineString, Point
from . import m2l_utils
import rasterio
from .m2l_enums import Datatype, VerboseLevel

import beartype
from .config import Config

######################################
# inspired by https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
#
# Simple Inverse Distance Weighting interpolation of observations z at x,y locations returned at locations defined by xi,yi arrays. From...
# scipy_idw(x, y, z, xi, yi)
# Args:
# x,y coordinates of points to be interpolated
# z value to be interpolated
# xi,yi grid of points where interpolation of z will be calculated - sci_py version of Simple Inverse Distance Weighting interpolation of observations z at x,y locations returned at locations defined by xi,yi arrays
#
# simple Inverse Distance Weighting calculation
######################################


def simple_idw(x, y, z, xi, yi):
    dist = distance_matrix(x, y, xi, yi)

    # In IDW, weights are 1 / distance
    weights = 1.0 / (dist)

    # Make weights sum to one
    weights /= weights.sum(axis=0)

    # Multiply the weights for each interpolated point by all observed Z-values
    zi = np.dot(weights.T, z)
    return zi


######################################
# call scipy inverse distance weighting
#
# Simple Inverse Distance Weighting interpolation of observations z at x,y locations returned at locations defined by xi,yi arrays. From...
# scipy_idw(x, y, z, xi, yi)
# Args:
# x,y coordinates of points to be interpolated
# z value to be interpolated
# xi,yi grid of points where interpolation of z will be calculated - sci_py version of Simple Inverse Distance Weighting interpolation of observations z at x,y locations returned at locations defined by xi,yi arrays
#
######################################


def scipy_idw(x, y, z, xi, yi):
    interp = Rbf(x, y, z, function="linear")
    return interp(xi, yi)


######################################
# call scipy Radial basis function interpolation
#
# scipy_rbf(x, y, z, xi, yi)
# Args:
# x,y coordinates of points to be interpolated
# z value to be interpolated
# xi,yi grid of points where interpolation of z will be calculated
#
# sci_py version of Radial Basis Function interpolation of observations z at x,y locations returned at locations defined by xi,yi arraysplot(x,y,z,grid)
######################################


def euclidean_norm_numpy(x1, x2):
    return np.linalg.norm(x1 - x2, axis=0)


# def scipy_rbf_dask(x, y, z, xi, yi):
#     # import dask.array as da
#     rbf = Rbf(x, y, z, function='multiquadric',
#               smooth=.15, norm=euclidean_norm_numpy)
#     n1 = len(xi)
#     xi = xi[np.newaxis, :]
#     yi = yi[np.newaxis, :]
#     ix = da.from_array(xi, chunks=(1, n1))
#     iy = da.from_array(yi, chunks=(1, n1))
#     iz = da.map_blocks(rbf, ix, iy)
#     zz = iz.compute()
#     print(zz)
#     return zz


def scipy_rbf(x, y, z, xi, yi):
    interp = Rbf(x, y, z, function="multiquadric", smooth=0.15)
    return interp(xi, yi)


def scipy_LNDI(x, y, z, xi, yi):  # actually LinearNDInterpolator
    from scipy.interpolate import LinearNDInterpolator

    interp = LinearNDInterpolator(list(zip(x, y)), z)
    return interp(xi, yi)


def scipy_CT(x, y, z, xi, yi):  # actually CloughTocher2DInterpolator
    from scipy.interpolate import CloughTocher2DInterpolator

    interp = CloughTocher2DInterpolator(list(zip(x, y)), z, rescale=True)
    return interp(xi, yi)


######################################
# calculate all distances between to arrays of points
# Make a distance matrix between pairwise observations
# Note: from <http://stackoverflow.com/questions/1871536>
# (Yay for ufuncs!)
# distance_matrix(x0, y0, x1, y1)
# Args:
# x0,y0 array of point locations
# x1,y1 second array of point locations
#
# Returns array of distances between all points defined by arrays by x0,y0 and all points defined by arrays x1,y1 from http://stackoverflow.com/questions/1871536
######################################


def distance_matrix(x0, y0, x1, y1):
    obs = np.vstack((x0, y0)).T
    interp = np.vstack((x1, y1)).T

    d0 = np.subtract.outer(obs[:, 0], interp[:, 0])
    d1 = np.subtract.outer(obs[:, 1], interp[:, 1])

    return np.hypot(d0, d1)


######################################
# plot an array of data
######################################


def plot(x, y, z, grid):
    plt.figure()
    plt.imshow(grid, extent=(0, 100, 0, 100), origin="lower")
    # plt.hold(True)
    # plt.scatter(x,100-y,c=z)
    plt.colorbar()


######################################
# switch function to select which intepolator to call
#
# interpolator_switch(calc, x, y, z, xi, yi)
# Args:
# calc string naming the interpolator to use
# x,y coordinates of points to be interpolated
# z value to be interpolated
# xi,yi grid of points where interpolation of z will be calculated - sci_py version of Simple Inverse Distance Weighting interpolation of observations z at x,y locations returned at locations defined by xi,yi arrays
#
######################################
def interpolator_switch(calc, x, y, z, xi, yi):
    if calc == "simple_idw":
        val = simple_idw(x, y, z, xi, yi)
    elif calc == "scipy_idw":
        val = scipy_idw(x, y, z, xi, yi)
    elif calc == "scipy_LNDI":
        val = scipy_LNDI(x, y, z, xi, yi)
    elif calc == "scipy_CT":
        val = scipy_CT(x, y, z, xi, yi)
    else:
        val = scipy_rbf(x, y, z, xi, yi)
    return val


######################################
# interpolate three data arrays using various schemes
#
# call_interpolator(calc,x,y,l,m,n,xi,yi,nx,ny,fault_flag)
# Args:
# calc calculation mode, one of 'simple_idw', 'scipy_idw', 'scipy_rbf'
# l,m,n arrays of direction cosines of pole to plane
# xi,yi arrays of locations of interpolated locations (assumes a grid for plotting, otherwise doesn't matter)
# nx,ny number of x,y elemnts in grid
# fault_flag toggle whether calc for near-fault orientations or not
#
# Call interpolator defined by calc for arrays of arbitrary location x,y located observations as triple or double arrays of 3D or 2D direction cosine arrays (l,m,n) and returns grid of nx ,ny interpolated values for points defined by xi,yi locations. Inspired by https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
######################################


def call_interpolator(calc, x, y, l, m, n, xi, yi, nx, ny, fault_flag):
    # Calculate IDW or other interpolators

    ZIl = interpolator_switch(calc, x, y, l, xi, yi)
    if not fault_flag:
        ZIl = ZIl.reshape((ny, nx))

    ZIm = interpolator_switch(calc, x, y, m, xi, yi)
    if not fault_flag:
        ZIm = ZIm.reshape((ny, nx))

    if type(n) is not int:
        ZIn = interpolator_switch(calc, x, y, n, xi, yi)
        if not fault_flag:
            ZIn = ZIn.reshape((ny, nx))
    else:
        ZIn = 0
    return (ZIl, ZIm, ZIn)


######################################
# Interpolate dipd,dipdirection data from shapefile
#
# interpolate_orientations(structure_file,tmp_path,bbox,c_l,use_gcode,scheme,gridx,gridy,fault_flag)
# Args:
# structure_file path to orientation layer
# tmp_path directory of temporary outputs from m2l
# bbox bounding box of region of interest
# c_l dictionary of codes and labels specific to input geo information layers
# use_gcode list of groups whose orientation data will be interpolated
# scheme interpolation scheme one of 'simple_idw', 'scipy_idw', 'scipy_rbf'
# gridx,gridy number of cols & rows in interpolation grid
# fault_flag toggle whether calc for near-fault orientations or not
#
# Interpolate orientation layer to produce regular grid of l,m,n direction cosines
# Can choose between various RBF and IDW options
# The purpose of these interpolations and associated code is to help in three cases:
# -- Providing estimated dips and contacts in fault-bounded domains where no structural data are available
# -- Needed to estimate true thickness of formations
# -- Useful for populating parts of maps where little structural data is available
######################################


def interpolate_orientations(
    structure_file, output_path, bbox, c_l, this_gcode, calc, gridx, gridy, fault_flag
):
    structure = gpd.read_file(structure_file, bbox=bbox)

    if len(this_gcode) == 1:
        # subset orientations to just those with this group
        is_gp = structure["GROUP"] == this_gcode
        gp_structure = structure[is_gp]
        # print('single group')
        # display(gp_structure)
    else:
        # print('first code',this_gcode[0])
        # subset orientations to just those with this group
        is_gp = structure["GROUP"] == this_gcode[0]
        gp_structure = structure[is_gp]
        gp_structure_all = gp_structure.copy()
        # print('first group')
        # display(gp_structure)

        for i in range(1, len(this_gcode)):
            # print('next code',this_gcode[i])
            # subset orientations to just those with this group
            is_gp = structure["GROUP"] == this_gcode[i]
            temp_gp_structure = structure[is_gp]
            gp_structure_all = pd.concat(
                [gp_structure_all, temp_gp_structure], ignore_index=True
            )
            # print('next group')
            # display(gp_structure)

    npts = len(gp_structure_all)

    if fault_flag:
        nx, ny = len(gridx), len(gridy)
    else:
        nx, ny = gridx, gridy
        xi = np.linspace(bbox[0], bbox[2], nx)
        yi = np.linspace(bbox[1], bbox[3], ny)
        xi, yi = np.meshgrid(xi, yi)
        xi, yi = xi.flatten(), yi.flatten()

    x = np.zeros(npts)
    y = np.zeros(npts)
    dip = np.zeros(npts)
    dipdir = np.zeros(npts)

    i = 0
    for a_pt in gp_structure_all.iterrows():
        x[i] = a_pt[1]["geometry"].x + (np.random.ranf() * 0.01)
        y[i] = a_pt[1]["geometry"].y + (np.random.ranf() * 0.01)
        dip[i] = a_pt[1]["DIP"]

        # All orientation have been converted to dipdir
        dipdir[i] = a_pt[1]["DIPDIR"]
        i = i + 1

    l = np.zeros(npts)
    m = np.zeros(npts)
    n = np.zeros(npts)

    for i in range(0, npts):
        l[i], m[i], n[i] = m2l_utils.ddd2dircos(dip[i], dipdir[i])

    if fault_flag:
        ZIl, ZIm, ZIn = call_interpolator(
            calc, x, y, l, m, n, gridx, gridy, nx, ny, fault_flag
        )
    else:
        ZIl, ZIm, ZIn = call_interpolator(
            calc, x, y, l, m, n, xi, yi, nx, ny, fault_flag
        )

    # Comparisons...
    if not fault_flag:
        plot(x, -y, l, ZIl)
        plt.title("l")
        plot(x, -y, m, ZIm)
        plt.title("m")
        plot(x, -y, n, ZIn)
        plt.title("n")

        plt.show()

    if fault_flag:
        f = open(os.path.join(output_path, "f_input.csv"), "w")
        fi = open(os.path.join(output_path, "f_interpolation_" + calc + ".csv"), "w")
        fl = open(os.path.join(output_path, "f_interpolation_l.csv"), "w")
        fm = open(os.path.join(output_path, "f_interpolation_m.csv"), "w")
        fn = open(os.path.join(output_path, "f_interpolation_n.csv"), "w")
    else:
        f = open(os.path.join(output_path, "input.csv"), "w")
        fi = open(os.path.join(output_path, "interpolation_" + calc + ".csv"), "w")
        fl = open(os.path.join(output_path, "interpolation_l.csv"), "w")
        fm = open(os.path.join(output_path, "interpolation_m.csv"), "w")
        fn = open(os.path.join(output_path, "interpolation_n.csv"), "w")

    f.write("x,y,dip,dipdirection\n")
    fi.write("x,y,dip,dipdirection\n")
    fl.write("x,y,l\n")
    fm.write("x,y,m\n")
    fn.write("x,y,n\n")

    for i in range(0, npts):
        ostr = "{},{},{},{}\n".format(x[i], y[i], int(dip[i]), int(dipdir[i]))
        # ostr=str(x[i])+","+str(y[i])+","+str(int(dip[i]))+","+str(int(dipdir[i]))+'\n'
        f.write(ostr)

    if fault_flag:
        for i in range(0, len(gridx)):
            L = ZIl[i] / (
                sqrt((pow(ZIl[i], 2.0)) + (pow(ZIm[i], 2.0)) + (pow(ZIn[i], 2.0)))
            )
            M = ZIm[i] / (
                sqrt((pow(ZIl[i], 2.0)) + (pow(ZIm[i], 2.0)) + (pow(ZIn[i], 2.0)))
            )
            N = ZIn[i] / (
                sqrt((pow(ZIl[i], 2.0)) + (pow(ZIm[i], 2.0)) + (pow(ZIn[i], 2.0)))
            )

            dip, dipdir = m2l_utils.dircos2ddd(L, M, N)
            if isnan(dip) or isnan(dipdir):
                dip = dipdir = L = M = N = 0
                print("Warning, no interpolated value for element No. ", i)
            ostr = "{},{},{},{}\n".format(gridx[i], gridy[i], int(dip), int(dipdir))
            # ostr=str(gridx[i])+","+str(gridy[i])+","+str(int(dip))+","+str(int(dipdir))+'\n'
            fi.write(ostr)

            ostr = "{},{},{}\n".format(gridx[i], gridy[i], L)
            # ostr=str(gridx[i])+","+str(gridy[i])+","+str(L)+'\n'
            fl.write(ostr)
            ostr = "{},{},{}\n".format(gridx[i], gridy[i], M)
            # ostr=str(gridx[i])+","+str(gridy[i])+","+str(M)+'\n'
            fm.write(ostr)
            ostr = "{},{},{}\n".format(gridx[i], gridy[i], N)
            # ostr=str(gridx[i])+","+str(gridy[i])+","+str(N)+'\n'
            fn.write(ostr)
    else:
        for xx in range(0, gridx):
            for yy in range(0, gridy):
                yyy = xx
                xxx = gridy - 2 - yy
                L = ZIl[xxx, yyy] / (
                    sqrt(
                        (pow(ZIl[xxx, yyy], 2.0))
                        + (pow(ZIm[xxx, yyy], 2.0))
                        + (pow(ZIn[xxx, yyy], 2.0))
                    )
                )
                M = ZIm[xxx, yyy] / (
                    sqrt(
                        (pow(ZIl[xxx, yyy], 2.0))
                        + (pow(ZIm[xxx, yyy], 2.0))
                        + (pow(ZIn[xxx, yyy], 2.0))
                    )
                )
                N = ZIn[xxx, yyy] / (
                    sqrt(
                        (pow(ZIl[xxx, yyy], 2.0))
                        + (pow(ZIm[xxx, yyy], 2.0))
                        + (pow(ZIn[xxx, yyy], 2.0))
                    )
                )

                dip, dipdir = m2l_utils.dircos2ddd(L, M, N)
                if isnan(dip) or isnan(dipdir):
                    dip = dipdir = L = M = N = 0
                    print(
                        "Warning, no interpolated value for grid point No. ",
                        xx,
                        ",",
                        yy,
                    )
                ostr = "{},{},{},{}\n".format(
                    bbox[0] + (xx * ((bbox[2] - bbox[0]) / gridx)),
                    bbox[1] + ((gridy - 1 - yy) * ((bbox[3] - bbox[1]) / gridy)),
                    int(dip),
                    int(dipdir),
                )
                # ostr=str(bbox[0]+(xx*((bbox[2]-bbox[0])/gridx)))+","+str(bbox[1]+((gridy-1-yy)*((bbox[3]-bbox[1])/gridy)))+","+str(int(dip))+","+str(int(dipdir))+'\n'
                fi.write(ostr)

                ostr = "{},{},{}\n".format(xx, yy, L)
                # ostr=str(xx)+","+str(yy)+","+str(L)+'\n'
                fl.write(ostr)
                ostr = "{},{},{}\n".format(xx, yy, M)
                # ostr=str(xx)+","+str(yy)+","+str(M)+'\n'
                fm.write(ostr)
                ostr = "{},{},{}\n".format(xx, yy, N)
                # ostr=str(xx)+","+str(yy)+","+str(N)+'\n'
                fn.write(ostr)

    f.close()
    fi.close()
    fl.close()
    fm.close()
    fn.close()

    if fault_flag:
        print(
            "orientations interpolated as dip dip direction",
            os.path.join(output_path, "f_interpolation_" + calc + ".csv"),
        )
        print(
            "orientations interpolated as l,m,n dir cos",
            os.path.join(output_path, "f_interpolation_l.csv"),
            " etc.",
        )
    else:
        fig, ax = plt.subplots(
            figsize=(10, 10),
        )
        ax.quiver(xi, yi, -ZIm, ZIl, headwidth=0)
        plt.show()
        print(
            "orientations interpolated as dip dip direction",
            os.path.join(output_path, "interpolation_" + calc + ".csv"),
        )
        print(
            "orientations interpolated as l,m,n dir cos",
            os.path.join(output_path, "interpolation_l.csv"),
            " etc.",
        )


######################################
# Interpolate 2D contact data from shapefile
#
# interpolate_contacts(geology_file,tmp_path,dtm,bbox,c_l,use_gcode,scheme,gridx,gridy,fault_flag)
# Args:
# geology_file path to basal contacts layer
# tmp_path directory of temporary outputs from m2l
# dtm rasterio format elevation grid
# bbox bounding box of region of interest
# c_l dictionary of codes and labels specific to input geo information layers
# use_gcode list of groups whose contact data will be interpolated
# scheme interpolation scheme one of 'simple_idw', 'scipy_idw', 'scipy_rbf'
# gridx,gridy number of cols & rows in interpolation grid
# fault_flag toggle whether calc for near-fault orientations or not
#
# Interpolate basal contacts layer to produce regular grid of l,m direction cosines
######################################


def interpolate_contacts(
    geology_file,
    output_path,
    dtm,
    dtb,
    dtb_null,
    cover_map,
    bbox,
    c_l,
    use_gcode,
    calc,
    gridx,
    gridy,
    fault_flag,
):
    geol_file = gpd.read_file(geology_file, bbox=bbox)
    # print(len(geol_file))
    # geol_file.plot( color='black',edgecolor='black')

    # Setup: Generate data...
    npts = 0
    decimate = 1
    if fault_flag:
        nx, ny = len(gridx), len(gridy)
    else:
        nx, ny = gridx, gridy
        xi = np.linspace(bbox[0], bbox[2], nx)
        yi = np.linspace(bbox[1], bbox[3], ny)
        xi, yi = np.meshgrid(xi, yi)
        xi, yi = xi.flatten(), yi.flatten()

    x = np.zeros(20000)  # FUDGE ################
    # should go through geology file to see how many contact
    y = np.zeros(20000)
    l = np.zeros(20000)  # segments will be made then define arrays?
    m = np.zeros(20000)

    if fault_flag:
        f = open(os.path.join(output_path, "f_raw_contacts.csv"), "w")
    else:
        f = open(os.path.join(output_path, "raw_contacts.csv"), "w")

    f.write("X,Y,Z,angle,lsx,lsy,formation,group\n")
    j = 0
    i = 0
    for (
        indx,
        acontact,
    ) in geol_file.iterrows():  # loop through distinct linestrings in MultiLineString
        if acontact.geometry.type == "MultiLineString":
            # print(i)
            for line in acontact.geometry:  # loop through line segments
                # print(i,len(acontact.geometry))
                if i % decimate == 0 and acontact["GROUP"] in use_gcode:
                    # if(acontact['id']==170):
                    # display(npts,line.coords[0][0],line.coords[1][0])
                    dlsx = line.coords[0][0] - line.coords[1][0]
                    dlsy = line.coords[0][1] - line.coords[1][1]
                    if (
                        not line.coords[0][0] == line.coords[1][0]
                        or not line.coords[0][1] == line.coords[1][1]
                    ):
                        lsx = dlsx / sqrt((dlsx * dlsx) + (dlsy * dlsy))
                        lsy = dlsy / sqrt((dlsx * dlsx) + (dlsy * dlsy))
                        x[i] = line.coords[1][0] + (dlsx / 2)
                        y[i] = line.coords[1][1] + (dlsy / 2)
                        angle = degrees(atan2(lsx, lsy))
                        l[i] = lsx
                        m[i] = lsy
                        # doesn't like point right on edge?
                        locations = [(x[i], y[i])]
                        height = m2l_utils.value_from_dtm_dtb(
                            dtm, dtb, dtb_null, cover_map, locations
                        )
                        if str(acontact["GROUP"]) == "None":
                            ostr = (
                                str(x[i])
                                + ","
                                + str(y[i])
                                + ","
                                + str(height)
                                + ","
                                + str(angle % 180)
                                + ","
                                + str(lsx)
                                + ","
                                + str(lsy)
                                + ","
                                + acontact["UNIT_NAME"]
                                .replace(" ", "_")
                                .replace("-", "_")
                                + ","
                                + acontact["UNIT_NAME"]
                                .replace(" ", "_")
                                .replace("-", "_")
                                + "\n"
                            )
                        else:
                            ostr = (
                                str(x[i])
                                + ","
                                + str(y[i])
                                + ","
                                + str(height)
                                + ","
                                + str(angle % 180)
                                + ","
                                + str(lsx)
                                + ","
                                + str(lsy)
                                + ","
                                + acontact["UNIT_NAME"]
                                .replace(" ", "_")
                                .replace("-", "_")
                                + ","
                                + acontact["GROUP"].replace(" ", "_").replace("-", "_")
                                + "\n"
                            )
                        f.write(ostr)
                        npts = npts + 1
                i = i + 1
        else:
            # display(acontact.geometry,acontact.geometry.coords)
            # for line in acontact: # loop through line segments in LineString
            if i % decimate == 0 and acontact["GROUP"] in use_gcode:
                dlsx = acontact.geometry.coords[0][0] - acontact.geometry.coords[1][0]
                dlsy = acontact.geometry.coords[0][1] - acontact.geometry.coords[1][1]
                if (
                    not acontact.geometry.coords[0][0] == acontact.geometry.coords[1][0]
                    or not acontact.geometry.coords[0][1]
                    == acontact.geometry.coords[1][1]
                ):
                    lsx = dlsx / sqrt((dlsx * dlsx) + (dlsy * dlsy))
                    lsy = dlsy / sqrt((dlsx * dlsx) + (dlsy * dlsy))
                    x[i] = acontact.geometry.coords[1][0] + (dlsx / 2)
                    y[i] = acontact.geometry.coords[1][1] + (dlsy / 2)
                    angle = degrees(atan2(lsx, lsy))
                    l[i] = lsx
                    m[i] = lsy
                    # doesn't like point right on edge?
                    locations = [(x[i], y[i])]
                    height = m2l_utils.value_from_dtm_dtb(
                        dtm, dtb, dtb_null, cover_map, locations
                    )
                    if str(acontact["GROUP"]) == "None":
                        ostr = (
                            str(x[i])
                            + ","
                            + str(y[i])
                            + ","
                            + str(height)
                            + ","
                            + str(angle % 180)
                            + ","
                            + str(lsx)
                            + ","
                            + str(lsy)
                            + ","
                            + acontact["UNIT_NAME"].replace(" ", "_").replace("-", "_")
                            + ","
                            + acontact["UNIT_NAME"].replace(" ", "_").replace("-", "_")
                            + "\n"
                        )
                    else:
                        ostr = (
                            str(x[i])
                            + ","
                            + str(y[i])
                            + ","
                            + str(height)
                            + ","
                            + str(angle % 180)
                            + ","
                            + str(lsx)
                            + ","
                            + str(lsy)
                            + ","
                            + acontact["UNIT_NAME"].replace(" ", "_").replace("-", "_")
                            + ","
                            + acontact["GROUP"].replace(" ", "_").replace("-", "_")
                            + "\n"
                        )
                    # print(ostr)
                    f.write(ostr)
                    # print(npts,dlsx,dlsy)
                    npts = npts + 1
                i = i + 1
        j = j + 1
    f.close()
    # print("i",i,"npts",npts)

    for i in range(0, npts):
        x[i] = x[i] + (np.random.ranf() * 0.01)
        y[i] = y[i] + (np.random.ranf() * 0.01)

    if fault_flag:
        ZIl, ZIm, ZIn = call_interpolator(
            calc,
            x[:npts],
            y[:npts],
            l[:npts],
            m[:npts],
            0,
            gridx,
            gridy,
            nx,
            ny,
            fault_flag,
        )
    else:
        ZIl, ZIm, ZIn = call_interpolator(
            calc, x[:npts], y[:npts], l[:npts], m[:npts], 0, xi, yi, nx, ny, fault_flag
        )

    # Comparisons...
    if not fault_flag:
        plot(x, -y, l, ZIl)
        plt.title("l")
        plot(x, -y, m, ZIm)
        plt.title("m")

    if fault_flag:
        fi = open(
            os.path.join(output_path, "f_interpolation_contacts_" + calc + ".csv"), "w"
        )
        fl = open(os.path.join(output_path, "f_interpolation_contacts_l.csv"), "w")
        fm = open(os.path.join(output_path, "f_interpolation_contacts_m.csv"), "w")
    else:
        fi = open(
            os.path.join(output_path, "interpolation_contacts_" + calc + ".csv"), "w"
        )
        fl = open(os.path.join(output_path, "interpolation_contacts_l.csv"), "w")
        fm = open(os.path.join(output_path, "interpolation_contacts_m.csv"), "w")

    fi.write("x,y,angle\n")
    fl.write("x,y,l\n")
    fm.write("x,y,m\n")

    if fault_flag:
        for i in range(0, len(gridx)):
            L = ZIl[i] / (sqrt((pow(ZIl[i], 2.0)) + (pow(ZIm[i], 2.0))))
            M = ZIm[i] / (sqrt((pow(ZIl[i], 2.0)) + (pow(ZIm[i], 2.0))))
            S = degrees(atan2(L, M))

            if isnan(S):
                S = 0
                print("Warning, no interpolated value for element No. ", i)

            ostr = str(gridx[i]) + "," + str(gridy[i]) + "," + str(int(S)) + "\n"
            fi.write(ostr)

            ostr = str(gridx[i]) + "," + str(gridy[i]) + "," + str(L) + "\n"
            fl.write(ostr)
            ostr = str(gridx[i]) + "," + str(gridy[i]) + "," + str(M) + "\n"
            fm.write(ostr)
    else:
        for xx in range(0, gridx):
            for yy in range(0, gridy):
                yyy = xx
                xxx = gridy - 2 - yy
                L = ZIl[xxx, yyy] / (
                    sqrt((pow(ZIl[xxx, yyy], 2.0)) + (pow(ZIm[xxx, yyy], 2.0)))
                )
                M = ZIm[xxx, yyy] / (
                    sqrt((pow(ZIl[xxx, yyy], 2.0)) + (pow(ZIm[xxx, yyy], 2.0)))
                )
                S = degrees(atan2(L, M))

                ostr = (
                    str(bbox[0] + (xx * ((bbox[2] - bbox[0]) / (gridx))))
                    + ","
                    + str(
                        bbox[1] + ((gridy - 2 - yy) * ((bbox[3] - bbox[1]) / (gridy)))
                    )
                    + ","
                    + str(int(S))
                    + "\n"
                )
                fi.write(ostr)

                ostr = str(xx) + "," + str(yy) + "," + str(L) + "\n"
                fl.write(ostr)
                ostr = str(xx) + "," + str(yy) + "," + str(M) + "\n"
                fm.write(ostr)

    fi.close()
    fl.close()
    fm.close()
    if fault_flag:
        print(
            "contacts interpolated as strike",
            os.path.join(output_path, "f_interpolation_contacts_" + calc + ".csv"),
        )
        print(
            "contacts interpolated as l,m dir cos",
            os.path.join(output_path, "f_interpolation_contacts_l.csv"),
            "etc.",
        )
    else:
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.quiver(xi, yi, ZIl, ZIm, headwidth=0)
        plt.show()
        print(
            "contacts interpolated as strike",
            os.path.join(output_path, "interpolation_contacts_" + calc + ".csv"),
        )
        print(
            "contacts interpolated as l,m dir cos",
            os.path.join(output_path, "interpolation_contacts_l.csv"),
            "etc.",
        )


######################################
# save all contacts as vectors (used for debugging)
#
# save_contact_vectors(geology_file,tmp_path,dtm,bbox,c_l,calc,decimate)
# Args:
# geology_file file path to geology polygons
# tmp_path directory of temporary outputs from m2l
# dtm rasterio format dtm raster
# bbox bounding box of model
# c_l dictionary of codes and labels specific to input geo information layers
# calc NOT USED
# decimate simple decimation factor for saving vectors
######################################


@beartype.beartype
def save_contact_vectors(config: Config, map_data, workflow: dict):
    geol_file = map_data.basal_contacts_no_faults
    dtm = map_data.get_map_data(Datatype.DTM).open()
    # geol_file = gpd.read_file(geology_file, bbox=config.bbox)
    # print(len(geol_file))
    # geol_file.plot( color='black',edgecolor='black')

    npts = 0
    i = 0
    for (
        indx,
        acontact,
    ) in geol_file.iterrows():  # loop through distinct linestrings in MultiLineString
        if acontact.geometry and acontact.geometry.type == "MultiLineString":
            for line in acontact.geometry.geoms:  # loop through line segments
                # if i % config.run_flags["contact_decimate"] == 0:
                #     npoint = 1
                i = i + 1
        else:
            # if i % config.run_flags["contact_decimate"] == 0:
            #     npoint = 1
            i = i + 1

    x = np.zeros(i + 1)
    y = np.zeros(i + 1)
    l = np.zeros(i + 1)
    m = np.zeros(i + 1)

    f = open(os.path.join(config.tmp_path, "raw_contacts.csv"), "w")
    f.write("X,Y,Z,angle,lsx,lsy,formation,group\n")
    j = 0
    i = 0
    for (
        indx,
        acontact,
    ) in geol_file.iterrows():  # loop through distinct linestrings in MultiLineString
        if acontact.geometry and acontact.geometry.type == "MultiLineString":
            # print(i)
            for line in acontact.geometry.geoms:  # loop through line segments
                # print(i,len(acontact.geometry))
                if i % config.run_flags["contact_decimate"] == 0:
                    # if(acontact['id']==170):
                    # display(npts,line.coords[0][0],line.coords[1][0])
                    dlsx = line.coords[0][0] - line.coords[1][0]
                    dlsy = line.coords[0][1] - line.coords[1][1]
                    if (
                        not line.coords[0][0] == line.coords[1][0]
                        or not line.coords[0][1] == line.coords[1][1]
                    ):
                        lsx = dlsx / sqrt((dlsx * dlsx) + (dlsy * dlsy))
                        lsy = dlsy / sqrt((dlsx * dlsx) + (dlsy * dlsy))
                        x[i] = line.coords[1][0] + (dlsx / 2)
                        y[i] = line.coords[1][1] + (dlsy / 2)
                        angle = degrees(atan2(lsx, lsy))
                        l[i] = lsx
                        m[i] = lsy
                        # doesn't like point right on edge?
                        locations = [(x[i], y[i])]
                        height = m2l_utils.value_from_dtm_dtb(
                            dtm,
                            map_data.dtb,
                            map_data.dtb_null,
                            workflow["cover_map"],
                            locations,
                        )

                        if str(acontact["GROUP"]) == "None":
                            ostr = "{},{},{},{},{},{},{},{}\n".format(
                                x[i],
                                y[i],
                                height,
                                angle % 180,
                                lsx,
                                lsy,
                                acontact["UNIT_NAME"]
                                .replace(" ", "_")
                                .replace("-", "_"),
                                acontact["UNIT_NAME"]
                                .replace(" ", "_")
                                .replace("-", "_"),
                            )
                            # ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact["UNIT_NAME"].replace(" ","_").replace("-","_")+","+acontact["UNIT_NAME"].replace(" ","_").replace("-","_")+"\n"
                        else:
                            ostr = "{},{},{},{},{},{},{},{}\n".format(
                                x[i],
                                y[i],
                                height,
                                angle % 180,
                                lsx,
                                lsy,
                                acontact["UNIT_NAME"]
                                .replace(" ", "_")
                                .replace("-", "_"),
                                acontact["GROUP"].replace(" ", "_").replace("-", "_"),
                            )
                            # ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact["UNIT_NAME"].replace(" ","_").replace("-","_")+","+acontact["GROUP"].replace(" ","_").replace("-","_")+"\n"
                        f.write(ostr)
                        npts = npts + 1
                i = i + 1
        else:
            # display(acontact.geometry,acontact.geometry.coords)
            # for line in acontact: # loop through line segments in LineString
            if acontact.geometry and i % config.run_flags["contact_decimate"] == 0:
                dlsx = acontact.geometry.coords[0][0] - acontact.geometry.coords[1][0]
                dlsy = acontact.geometry.coords[0][1] - acontact.geometry.coords[1][1]
                if (
                    not acontact.geometry.coords[0][0] == acontact.geometry.coords[1][0]
                    or not acontact.geometry.coords[0][1]
                    == acontact.geometry.coords[1][1]
                ):
                    lsx = dlsx / sqrt((dlsx * dlsx) + (dlsy * dlsy))
                    lsy = dlsy / sqrt((dlsx * dlsx) + (dlsy * dlsy))
                    x[i] = acontact.geometry.coords[1][0] + (dlsx / 2)
                    y[i] = acontact.geometry.coords[1][1] + (dlsy / 2)
                    angle = degrees(atan2(lsx, lsy))
                    l[i] = lsx
                    m[i] = lsy
                    # doesn't like point right on edge?
                    locations = [(x[i], y[i])]
                    height = m2l_utils.value_from_dtm_dtb(
                        dtm,
                        map_data.dtb,
                        map_data.dtb_null,
                        workflow["cover_map"],
                        locations,
                    )
                    if str(acontact["GROUP"]) == "None":
                        ostr = "{},{},{},{},{},{},{},{}\n".format(
                            x[i],
                            y[i],
                            height,
                            angle % 180,
                            lsx,
                            lsy,
                            acontact["UNIT_NAME"].replace(" ", "_").replace("-", "_"),
                            acontact["UNIT_NAME"].replace(" ", "_").replace("-", "_"),
                        )
                        # ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact["UNIT_NAME"].replace(" ","_").replace("-","_")+","+acontact["UNIT_NAME"].replace(" ","_").replace("-","_")+"\n"
                    else:
                        ostr = "{},{},{},{},{},{},{},{}\n".format(
                            x[i],
                            y[i],
                            height,
                            angle % 180,
                            lsx,
                            lsy,
                            acontact["UNIT_NAME"].replace(" ", "_").replace("-", "_"),
                            acontact["GROUP"].replace(" ", "_").replace("-", "_"),
                        )
                        # ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact["UNIT_NAME"].replace(" ","_").replace("-","_")+","+acontact["GROUP"].replace(" ","_").replace("-","_")+"\n"
                    # print(ostr)
                    f.write(ostr)
                    # print(npts,dlsx,dlsy)
                    npts = npts + 1
                i = i + 1
        j = j + 1
    f.close()
    if config.verbose_level != VerboseLevel.NONE:
        print(
            npts, "points saved to", os.path.join(config.tmp_path, "raw_contacts.csv")
        )


####################################
# combine interpolated contact information (to provide l,m with interpolated dip,dipdirection data (to provide n)
#
# join_contacts_and_orientations(combo_file,geology_file,tmp_path,dtm_reproj_file,c_l,lo,mo,no,lc,mc,xy,dst_crs,bbox,fault_flag)
# combo_file path to temporary combined information geology_file path to basal contacts layer
# tmp_path directory of temporary outputs from m2l
# dtm_reproj_file path to reprojected dtm file
# c_l dictionary of codes and labels specific to input geo information layers
# lo,mo,no 3D direction cosines of interpolated orientations
# lc,mc 2D direction cosines of interpolated contacts
# xy interpolated orientations (used to get x,y locations only) dst_crs Coordinate Reference System of destination geotif (any length-based projection)
# bbox bounding box of region of interest
# fault_flag toggle whether calc for near-fault orientations or not
#
# Combine interpolation orientations with interpolated basal contacts layers to produce regular grid of interpolated dip, dip direction estimates
# Uses normalised direction cosines (l,m,n):
# -- l,m from RBF of basal contact orientations -- signs of l & m from misorientation with RBF of orientation data and -- n from RBF of orientation data
#
# Useful for adding data where no orientations are available (e.g. in fault bounded domains) and for calculating true thickness of layers. Assumes a 2D plane of data, but if 3D RBF was calulated and projected contact info was used it should apply with topography too.
####################################
def join_contacts_and_orientations(
    combo_file,
    geology_file,
    output_path,
    dtm_reproj_file,
    dtb,
    dtb_null,
    cover_map,
    c_l,
    lo,
    mo,
    no,
    lc,
    mc,
    xy,
    dst_crs,
    bbox,
    fault_flag,
):
    f = open(combo_file, "w")
    f.write("x,y,dip,dipdirection,misorientation,dotproduct\n")

    for i in range(0, len(lc)):
        # print(mc[i,2],lc[i,2],lo[i,2],mo[i,2],no[i,2])
        # scaling contact dircos to *include* dip info
        scale = sqrt(1 - pow(no[i, 2], 2))
        # includes 90 rotation to account for orthogonality of contact and dip direction
        lcscaled = scale * -mc[i, 2]
        mcscaled = scale * lc[i, 2]
        # scaling dip dipdir dircos to *exclude* dip info
        scale2 = sqrt(pow(lo[i, 2], 2) + pow(mo[i, 2], 2))
        if scale2 > 0.0:
            loscaled = lo[i, 2] / scale2
            moscaled = mo[i, 2] / scale2
        else:
            loscaled = 0
            moscaled = 0
        # includes 90 rotation to account for orthogonality of contact and dip direction
        dotproduct = (-mc[i, 2] * loscaled) + (lc[i, 2] * moscaled)
        if dotproduct < 0:
            lcscaled = -lcscaled
            mcscaled = -mcscaled

        misorientation = degrees(acos(dotproduct))
        dip, dipdir = m2l_utils.dircos2ddd(lcscaled, mcscaled, no[i, 2])
        ostr = "{},{},{},{},{},{}\n".format(
            xy[i, 0], xy[i, 1], int(dip), int(dipdir), int(misorientation), dotproduct
        )
        # ostr=str(xy[i,0])+','+str(xy[i,1])+','+str(int(dip))+','+str(int(dipdir))+','+str(int(misorientation))+','+str(dotproduct)+'\n'
        f.write(ostr)
    f.close()

    geology = gpd.read_file(geology_file, bbox=bbox)
    geology.crs = dst_crs
    # geology = m2l_utils.explode(geology)
    geology = geology.explode(ignore_index=True)

    data = pd.read_csv(combo_file)

    geometry = [Point(xy) for xy in zip(data["x"], data["y"])]

    gdf = GeoDataFrame(data, crs=dst_crs, geometry=geometry)

    gdf.crs = dst_crs
    # print(gdf.crs, geology.crs)
    structure_code = gpd.sjoin(gdf, geology, how="left", predicate="within")
    dtm = rasterio.open(dtm_reproj_file)
    if fault_flag:
        f = open(os.path.join(output_path, "f_combo_full.csv"), "w")
    else:
        f = open(os.path.join(output_path, "combo_full.csv"), "w")
    f.write("X,Y,Z,azimuth,dip,polarity,formation\n")
    # last_code = ""
    for indx, a_point in structure_code.iterrows():

        locations = [(a_point["x"], a_point["y"])]
        height = m2l_utils.value_from_dtm_dtb(dtm, dtb, dtb_null, cover_map, locations)
        ostr = str(a_point["x"]) + ","
        ostr = ostr + str(a_point["y"]) + ","
        ostr = ostr + str(height) + "," + str(int(a_point["dipdirection"])) + ","
        ostr = ostr + str(int(a_point["dip"])) + ",1,"
        ostr = (
            ostr + str(a_point["UNIT_NAME"]).replace("-", "_").replace(" ", "_") + "\n"
        )

        if not str(a_point["UNIT_NAME"]) == "nan":
            f.write(ostr)
        # last_code = a_point["UNIT_NAME"]
    f.close()
    if fault_flag:
        print(
            "contacts and orientations interpolated as dip dip direction",
            os.path.join(output_path, "f_combo_full.csv"),
        )
    else:
        print(
            "contacts and orientations interpolated as dip dip direction",
            os.path.join(output_path, "combo_full.csv"),
        )


######################################
# Interpolate dipd,dipdirection data from shapefile usin fold axial traces as additional constraints
# interpolate_orientations_with_fat(structure_file,output_path,bbox,c_l,this_gcode,calc,gridx,gridy)
# structure_file path to orientation layer
# output_path directory for outputs from m2l
# bbox bounding box of region of interest
# c_l dictionary of codes and labels specific to input geo information layers
# this_gcode list of groups whose orientation data will be interpolated
# calc interpolation scheme one of 'simple_idw', 'scipy_idw', 'scipy_rbf'
# gridx,gridy number of cols & rows in interpolation grid
#
# Interpolate orientation layer to produce regular grid of l,m,n direction cosines
# Can choose between various RBF and IDW options
# The purpose of these interpolations and associated code is to help in three cases:
# -- Providing estimated dips and contacts in fault-bounded domains where no structural data are available
# -- Needed to estimate true thickness of formations
# -- Useful for poulating parts of maps where little structural data is available
######################################


def interpolate_orientations_with_fat(
    structure_file, output_path, bbox, c_l, this_gcode, calc, gridx, gridy
):
    structure = gpd.read_file(structure_file, bbox=bbox)
    fat_orientations = pd.read_csv(
        os.path.join(output_path, "fold_axial_trace_orientations2.csv"), sep=","
    )

    if len(this_gcode) == 1:
        # subset orientations to just those with this group
        is_gp = structure["GROUP"] == this_gcode
        gp_structure = structure[is_gp]
        print("single group")
        # display(gp_structure)
    else:
        print("first code", this_gcode[0])
        # subset orientations to just those with this group
        is_gp = structure["GROUP"] == this_gcode[0]
        gp_structure = structure[is_gp]
        gp_structure_all = gp_structure.copy()
        print("first group")
        # display(gp_structure)

        for i in range(1, len(this_gcode)):
            print("next code", this_gcode[i])
            # subset orientations to just those with this group
            is_gp = structure["GROUP"] == this_gcode[i]
            temp_gp_structure = structure[is_gp]
            gp_structure_all = pd.concat(
                [gp_structure_all, temp_gp_structure], ignore_index=True
            )
            print("next group")
            # display(gp_structure)

    npts = len(gp_structure_all) + len(fat_orientations)

    nx, ny = gridx, gridy

    xi = np.linspace(bbox[0], bbox[2], nx)
    yi = np.linspace(bbox[1], bbox[3], ny)
    xi, yi = np.meshgrid(xi, yi)
    xi, yi = xi.flatten(), yi.flatten()
    x = np.zeros(npts)
    y = np.zeros(npts)
    dip = np.zeros(npts)
    dipdir = np.zeros(npts)

    i = 0
    for ind, a_pt in gp_structure_all.iterrows():
        x[i] = a_pt["geometry"].x
        y[i] = a_pt["geometry"].y
        dip[i] = a_pt["DIP"]
        # All orientation have been converted to dipdir
        dipdir[i] = a_pt["DIPDIR"]
        i = i + 1

    for ind, a_pt in fat_orientations.iterrows():
        x[i] = a_pt["X"]
        y[i] = a_pt["Y"]
        dip[i] = a_pt["dip"]
        dipdir[i] = a_pt["azimuth"]
        i = i + 1

    l = np.zeros(npts)
    m = np.zeros(npts)
    n = np.zeros(npts)

    for i in range(0, npts):
        l[i], m[i], n[i] = m2l_utils.ddd2dircos(dip[i], dipdir[i])

    ZIl, ZIm, ZIn = call_interpolator(calc, x, y, l, m, n, xi, yi, nx, ny)

    # Comparisons...
    plot(x, -y, l, ZIl)
    plt.title("l")
    plot(x, -y, m, ZIm)
    plt.title("m")
    plot(x, -y, n, ZIn)
    plt.title("n")

    plt.show()

    f = open(os.path.join(output_path, "input.csv"), "w")
    fi = open(os.path.join(output_path, "interpolation_" + calc + ".csv"), "w")
    fl = open(os.path.join(output_path, "interpolation_l.csv"), "w")
    fm = open(os.path.join(output_path, "interpolation_m.csv"), "w")
    fn = open(os.path.join(output_path, "interpolation_n.csv"), "w")

    f.write("x,y,dip,dipdirection\n")
    fi.write("x,y,dip,dipdirection\n")
    fl.write("x,y,l\n")
    fm.write("x,y,m\n")
    fn.write("x,y,n\n")

    for i in range(0, npts):
        ostr = "{},{},{},{}\n".format(x[i], y[i], int(dip[i]), int(dipdir[i]))
        # ostr=str(x[i])+","+str(y[i])+","+str(int(dip[i]))+","+str(int(dipdir[i]))+'\n'
        f.write(ostr)

    for xx in range(0, gridx):
        for yy in range(0, gridy):
            yyy = xx
            xxx = gridy - 2 - yy
            L = ZIl[xxx, yyy] / (
                sqrt(
                    (pow(ZIl[xxx, yyy], 2.0))
                    + (pow(ZIm[xxx, yyy], 2.0))
                    + (pow(ZIn[xxx, yyy], 2.0))
                )
            )
            M = ZIm[xxx, yyy] / (
                sqrt(
                    (pow(ZIl[xxx, yyy], 2.0))
                    + (pow(ZIm[xxx, yyy], 2.0))
                    + (pow(ZIn[xxx, yyy], 2.0))
                )
            )
            N = ZIn[xxx, yyy] / (
                sqrt(
                    (pow(ZIl[xxx, yyy], 2.0))
                    + (pow(ZIm[xxx, yyy], 2.0))
                    + (pow(ZIn[xxx, yyy], 2.0))
                )
            )

            dip, dipdir = m2l_utils.dircos2ddd(L, M, N)
            ostr = "{},{},{},{}\n".format(
                bbox[0] + (xx * ((bbox[2] - bbox[0]) / gridx)),
                bbox[1] + ((gridy - 1 - yy) * ((bbox[3] - bbox[1]) / gridy)),
                int(dip),
                int(dipdir),
            )
            # ostr=str(bbox[0]+(xx*((bbox[2]-bbox[0])/gridx)))+","+str(bbox[1]+((gridy-1-yy)*((bbox[3]-bbox[1])/gridy)))+","+str(int(dip))+","+str(int(dipdir))+'\n'
            fi.write(ostr)

            ostr = "{},{},{}\n".format(xx, yy, L)
            # ostr=str(xx)+","+str(yy)+","+str(L)+'\n'
            fl.write(ostr)
            ostr = "{},{},{}\n".format(xx, yy, M)
            # ostr=str(xx)+","+str(yy)+","+str(M)+'\n'
            fm.write(ostr)
            ostr = "{},{},{}\n".format(xx, yy, N)
            # ostr=str(xx)+","+str(yy)+","+str(N)+'\n'
            fn.write(ostr)

    f.close()
    fi.close()
    fl.close()
    fm.close()
    fn.close()

    fig, ax = plt.subplots(
        figsize=(10, 10),
    )
    ax.quiver(xi, yi, -ZIm, ZIl, headwidth=0)
    plt.show()
    print(
        "orientations interpolated as dip dip direction",
        os.path.join(output_path, "interpolation_" + calc + ".csv"),
    )
    print(
        "orientations interpolated as l,m,n dir cos",
        os.path.join(output_path, "interpolation_l.csv"),
        "etc.",
    )


####################################################
# For each fault string:
# process_fault_throw_and_near_orientations(tmp_path,output_path,dtm_reproj_file,c_l,use_gcode,use_gcode2,dst_crs,bbox,scheme)
# Args:
#
#    incementally advance along polyline every at each inter-node (no point in doing more?)
#    find local stratigraphy 10m to left and right of fault
# Once full fault has been traversed:
#
#    Find list of contacts left
#    Find equivalent contacts on right
#    use interpolated orientations to estimate minimum true offset assuming vertical displacement and store
#    if no equivalent found, flag as domain fault and find min strat offset for contact, use cumulative minimum thickness estimate and store with flag (not implemented)
#    estimate median & sd of minimum fault offset and store with flag (not implemented)
# Local Orientations: Since much of the code is the same, we benefit by calculating local orientation data either side of
# fault so that geomodeller/gempy have satisfied fault compartment orientation data
###################################################


def process_fault_throw_and_near_orientations(
    tmp_path,
    output_path,
    dtm_reproj_file,
    dtb,
    dtb_null,
    cover_map,
    c_l,
    use_gcode,
    use_gcode2,
    dst_crs,
    bbox,
    scheme,
):
    fault_file = os.path.join(tmp_path, "faults_clip.shp")
    geology_file = os.path.join(tmp_path, "geol_clip.shp")

    faults = gpd.read_file(fault_file)
    geology = gpd.read_file(geology_file)
    dtm = rasterio.open(dtm_reproj_file)

    all_long_faults = np.genfromtxt(
        os.path.join(output_path, "fault_dimensions.csv"), delimiter=",", dtype="U100"
    )
    fault_names = all_long_faults[1:, :1]
    m_step = 10.0  # outstep from fault
    # if we have two faults at low angle this can be a problem as the point crosses the next fault
    xi = []
    yi = []
    fdc = []
    all_coordsdist = []
    all_coords_x = []
    all_coords_y = []

    fftc = open(os.path.join(output_path, "fault_tip_contacts.csv"), "w")
    fftc.write("X,Y,Z,formation\n")

    # loop through all faults

    for index, fault in faults.iterrows():
        # if(fault["GEOMETRY_OBJECT_ID"]!=1071):
        # continue
        if fault.geometry.type == "LineString":
            # in is dangerous as Fault_1 is in Fault_10
            if "Fault_" + str(fault["GEOMETRY_OBJECT_ID"]) in fault_names:
                # print('LineString','Fault_'+str(fault["GEOMETRY_OBJECT_ID"]),len(fault.geometry.coords))
                lcoords = []
                rcoords = []
                index = []

                # make a list of points just offset from mid-points between fault nodes and join
                # geology polygon information to points
                j = 0
                for i in range(0, len(fault.geometry.coords) - 1):
                    for inc in np.arange(0.01, 1, 0.01):
                        midx = fault.geometry.coords[i][0] + (
                            (
                                fault.geometry.coords[i + 1][0]
                                - fault.geometry.coords[i][0]
                            )
                            * inc
                        )
                        midy = fault.geometry.coords[i][1] + (
                            (
                                fault.geometry.coords[i + 1][1]
                                - fault.geometry.coords[i][1]
                            )
                            * inc
                        )
                        l, m = m2l_utils.pts2dircos(
                            fault.geometry.coords[i][0],
                            fault.geometry.coords[i][1],
                            fault.geometry.coords[i + 1][0],
                            fault.geometry.coords[i + 1][1],
                        )
                        lcoords.append([(midx + (m_step * m), midy - (m_step * l))])
                        rcoords.append([(midx - (m_step * m), midy + (m_step * l))])
                        index.append([(j)])
                        j = j + 1
                lgeom = [Point(xy) for xy in lcoords]
                rgeom = [Point(xy) for xy in rcoords]

                lgdf = GeoDataFrame(index, crs=dst_crs, geometry=lgeom)
                rgdf = GeoDataFrame(index, crs=dst_crs, geometry=rgeom)
                lcode = gpd.sjoin(lgdf, geology, how="left", predicate="within")
                rcode = gpd.sjoin(rgdf, geology, how="left", predicate="within")
                # display(lcode)

                # add points to list if they have different geology code than previous node on left side

                first = True
                lcontact = []
                lastlcode = ""

                for ind, indl in lcode.iterrows():

                    if ind < len(lcode) and not isnan(indl["index_right"]):
                        ntest1 = str(indl["DESCRIPTION"])
                        ntest2 = str(indl["ROCKTYPE1"])

                        if not ntest1 == "None" and not ntest2 == "None":
                            if ind == 1 or (
                                not lastlcode == indl["UNIT_NAME"]
                                and (
                                    (not c_l["sill"] in indl["DESCRIPTION"])
                                    or (not c_l["intrusive"] in indl["ROCKTYPE1"])
                                )
                            ):
                                all_coords_x.append(indl.geometry.x)
                                all_coords_y.append(indl.geometry.y)
                                if first:
                                    first = False
                                    firstlx = indl.geometry.x
                                    firstly = indl.geometry.y
                                    firstlc = (
                                        indl["UNIT_NAME"]
                                        .replace(" ", "_")
                                        .replace("-", "_")
                                    )
                                lastlx = indl.geometry.x
                                lastly = indl.geometry.y
                                lastlc = (
                                    indl["UNIT_NAME"]
                                    .replace(" ", "_")
                                    .replace("-", "_")
                                )

                            if lastlcode == "" and (
                                (not c_l["sill"] in indl["DESCRIPTION"])
                                or (not c_l["intrusive"] in indl["ROCKTYPE1"])
                            ):
                                lastlcode = indl["UNIT_NAME"]

                            # print('l',ind,indl["UNIT_NAME"],indl["DESCRIPTION"],indl["ROCKTYPE1"])
                            if (
                                not ntest1 == "None"
                                and not ntest2 == "None"
                                and not str(indl["UNIT_NAME"]) == "nan"
                            ):
                                if (not indl["UNIT_NAME"] == lastlcode) and (
                                    (not c_l["sill"] in indl["DESCRIPTION"])
                                    or (not c_l["intrusive"] in indl["ROCKTYPE1"])
                                ):
                                    lcontact.append(
                                        [
                                            (
                                                ind,
                                                lastlcode,
                                                indl["UNIT_NAME"],
                                                indl.geometry,
                                            )
                                        ]
                                    )
                                    lastlcode = indl["UNIT_NAME"]

                # add points to list if they have different geology code than previous node on right side

                first = True
                rcontact = []
                lastrcode = ""
                for ind, indr in rcode.iterrows():
                    if ind < len(rcode) and not isnan(indr["index_right"]):
                        ntest1 = str(indr["DESCRIPTION"])
                        ntest2 = str(indr["ROCKTYPE1"])

                        if not ntest1 == "None" and not ntest2 == "None":
                            if ind == 1 or (
                                not lastrcode == indr["UNIT_NAME"]
                                and (
                                    (not c_l["sill"] in indr["DESCRIPTION"])
                                    or (not c_l["intrusive"] in indr["ROCKTYPE1"])
                                )
                            ):
                                all_coords_x.append(indr.geometry.x)
                                all_coords_y.append(indr.geometry.y)
                                if first:
                                    first = False
                                    firstrx = indr.geometry.x
                                    firstry = indr.geometry.y
                                    firstrc = (
                                        indr["UNIT_NAME"]
                                        .replace(" ", "_")
                                        .replace("-", "_")
                                    )
                                lastrx = indr.geometry.x
                                lastry = indr.geometry.y
                                lastrc = (
                                    indr["UNIT_NAME"]
                                    .replace(" ", "_")
                                    .replace("-", "_")
                                )

                            if lastrcode == "" and (
                                (not c_l["sill"] in indr["DESCRIPTION"])
                                or (not c_l["intrusive"] in indr["ROCKTYPE1"])
                            ):
                                lastrcode = indr["UNIT_NAME"]
                            # print('r',ind,indr["UNIT_NAME"],indr["DESCRIPTION"],indr["ROCKTYPE1"])

                            # print(lastrcode,ntest1,ntest2,str(indr["UNIT_NAME"]),indr["UNIT_NAME"],c_l['sill'],c_l['intrusive'])

                            if (
                                not ntest1 == "None"
                                and not ntest2 == "None"
                                and not str(indr["UNIT_NAME"]) == "nan"
                            ):
                                if (not indr["UNIT_NAME"] == lastrcode) and (
                                    (not c_l["sill"] in indr["DESCRIPTION"])
                                    or (not c_l["intrusive"] in indr["ROCKTYPE1"])
                                ):
                                    rcontact.append(
                                        [
                                            (
                                                ind,
                                                lastrcode,
                                                indr["UNIT_NAME"],
                                                indr.geometry,
                                            )
                                        ]
                                    )
                                    lastrcode = indr["UNIT_NAME"]

                locations = [(firstlx, firstly)]
                first_height_l = m2l_utils.value_from_dtm_dtb(
                    dtm, dtb, dtb_null, cover_map, locations
                )
                ostr = "{},{},{},{}\n".format(firstlx, firstly, first_height_l, firstlc)
                fftc.write(ostr)
                locations = [(firstrx, firstry)]
                first_height_r = m2l_utils.value_from_dtm_dtb(
                    dtm, dtb, dtb_null, cover_map, locations
                )
                ostr = "{},{},{},{}\n".format(firstrx, firstry, first_height_r, firstrc)
                fftc.write(ostr)
                locations = [(lastlx, lastly)]
                last_height_l = m2l_utils.value_from_dtm_dtb(
                    dtm, dtb, dtb_null, cover_map, locations
                )
                ostr = "{},{},{},{}\n".format(lastlx, lastly, last_height_l, lastlc)
                fftc.write(ostr)
                locations = [(lastrx, lastry)]
                last_height_r = m2l_utils.value_from_dtm_dtb(
                    dtm, dtb, dtb_null, cover_map, locations
                )
                ostr = "{},{},{},{}\n".format(lastrx, lastry, last_height_r, lastrc)
                fftc.write(ostr)

                # loop through left and right sides to find equivalent contact pairs along fault

                if len(lcontact) > 0 and len(rcontact) > 0:
                    for lc in lcontact:
                        for rc in rcontact:
                            # display('l',lc[0][3].x,'r',rc[0][3].x)
                            if (
                                lc[0][1] == rc[0][1]
                                and lc[0][2] == rc[0][2]
                                and not lc[0][1] == ""
                            ):
                                dist = m2l_utils.ptsdist(
                                    lc[0][3].x, lc[0][3].y, rc[0][3].x, rc[0][3].y
                                )
                                if lc[0][0] < rc[0][0]:
                                    dist = -dist
                                # print('***',lc,rc)

                                xi.append((lc[0][3].x))
                                yi.append((lc[0][3].y))
                                l, m = m2l_utils.pts2dircos(
                                    lc[0][3].x, lc[0][3].y, rc[0][3].x, rc[0][3].y
                                )
                                if not (l == 0.0 and m == 0.0):
                                    fdc.append(
                                        (
                                            l,
                                            m,
                                            "Fault_" + str(fault["GEOMETRY_OBJECT_ID"]),
                                        )
                                    )
                                    all_coordsdist.append((dist))
        else:
            # in is dangerous as Fault_1 is in Fault_10
            if "Fault_" + str(fault["GEOMETRY_OBJECT_ID"]) in fault_names:
                for fls in fault.geometry:
                    fault_ls = LineString(fls)
                    lcoords = []
                    rcoords = []
                    index = []
                    # display("MLS DEBUG",fault.geometry.type)

                    j = 0
                    for i in range(0, len(fault_ls.coords) - 1):
                        for inc in np.arange(0.01, 1, 0.01):
                            midx = fault_ls.coords[i][0] + (
                                (fault_ls.coords[i + 1][0] - fault_ls.coords[i][0])
                                * inc
                            )
                            midy = fault_ls.coords[i][1] + (
                                (fault_ls.coords[i + 1][1] - fault_ls.coords[i][1])
                                * inc
                            )
                            l, m = m2l_utils.pts2dircos(
                                fault_ls.coords[i][0],
                                fault_ls.coords[i][1],
                                fault_ls.coords[i + 1][0],
                                fault_ls.coords[i + 1][1],
                            )
                            lcoords.append([(midx + (m_step * m), midy - (m_step * l))])
                            rcoords.append([(midx - (m_step * m), midy + (m_step * l))])
                            index.append([(j)])
                            j = j + 1

                    lgeom = [Point(xy) for xy in lcoords]
                    rgeom = [Point(xy) for xy in rcoords]
                    lgdf = GeoDataFrame(index, crs=dst_crs, geometry=lgeom)
                    rgdf = GeoDataFrame(index, crs=dst_crs, geometry=rgeom)
                    lcode = gpd.sjoin(lgdf, geology, how="left", predicate="within")
                    rcode = gpd.sjoin(rgdf, geology, how="left", predicate="within")

                    # add points to list if they have different geology code than previous node on left side

                    first = True
                    lcontact = []
                    lastlcode = ""
                    for ind, indl in lcode.iterrows():

                        if ind < len(lcode) and not isnan(indl["index_right"]):
                            ntest1 = str(indl["DESCRIPTION"])
                            ntest2 = str(indl["ROCKTYPE1"])

                            if not ntest1 == "None" and not ntest2 == "None":
                                if ind == 1 or (
                                    not lastlcode == indl["UNIT_NAME"]
                                    and (
                                        (not c_l["sill"] in indl["DESCRIPTION"])
                                        or (not c_l["intrusive"] in indl["ROCKTYPE1"])
                                    )
                                ):
                                    all_coords_x.append(indl.geometry.x)
                                    all_coords_y.append(indl.geometry.y)
                                    if first:
                                        first = False
                                        firstlx = indl.geometry.x
                                        firstly = indl.geometry.y
                                        firstlc = (
                                            indl["UNIT_NAME"]
                                            .replace(" ", "_")
                                            .replace("-", "_")
                                        )
                                    lastlx = indl.geometry.x
                                    lastly = indl.geometry.y
                                    lastlc = (
                                        indl["UNIT_NAME"]
                                        .replace(" ", "_")
                                        .replace("-", "_")
                                    )

                                if lastlcode == "" and (
                                    (not c_l["sill"] in indl["DESCRIPTION"])
                                    or (not c_l["intrusive"] in indl["ROCKTYPE1"])
                                ):
                                    lastlcode = indl["UNIT_NAME"]

                                if (
                                    not ntest1 == "None"
                                    and not ntest2 == "None"
                                    and not str(indl["UNIT_NAME"]) == "nan"
                                ):
                                    if (not indl["UNIT_NAME"] == lastlcode) and (
                                        (not c_l["sill"] in indl["DESCRIPTION"])
                                        or (not c_l["intrusive"] in indl["ROCKTYPE1"])
                                    ):
                                        lcontact.append(
                                            [
                                                (
                                                    ind,
                                                    lastlcode,
                                                    indl["UNIT_NAME"],
                                                    indl.geometry,
                                                )
                                            ]
                                        )
                                        lastlcode = indl["UNIT_NAME"]

                    # add points to list if they have different geology code than previous node on right side

                    first = True
                    rcontact = []
                    lastrcode = ""
                    for ind, indr in rcode.iterrows():
                        if ind < len(rcode) and not isnan(indr["index_right"]):
                            ntest1 = str(indr["DESCRIPTION"])
                            ntest2 = str(indr["ROCKTYPE1"])

                            if not ntest1 == "None" and not ntest2 == "None":
                                if ind == 1 or (
                                    not lastrcode == indr["UNIT_NAME"]
                                    and (
                                        (not c_l["sill"] in indr["DESCRIPTION"])
                                        or (not c_l["intrusive"] in indr["ROCKTYPE1"])
                                    )
                                ):
                                    all_coords_x.append(indr.geometry.x)
                                    all_coords_y.append(indr.geometry.y)
                                    if first:
                                        first = False
                                        firstrx = indr.geometry.x
                                        firstry = indr.geometry.y
                                        firstrc = (
                                            indr["UNIT_NAME"]
                                            .replace(" ", "_")
                                            .replace("-", "_")
                                        )
                                    lastrx = indr.geometry.x
                                    lastry = indr.geometry.y
                                    lastrc = (
                                        indr["UNIT_NAME"]
                                        .replace(" ", "_")
                                        .replace("-", "_")
                                    )

                                if lastrcode == "" and (
                                    (not c_l["sill"] in indr["DESCRIPTION"])
                                    or (not c_l["intrusive"] in indr["ROCKTYPE1"])
                                ):
                                    lastrcode = indr["UNIT_NAME"]

                                if (
                                    not ntest1 == "None"
                                    and not ntest2 == "None"
                                    and not str(indr["UNIT_NAME"]) == "nan"
                                ):
                                    if (not indr["UNIT_NAME"] == lastrcode) and (
                                        (not c_l["sill"] in indr["DESCRIPTION"])
                                        or (not c_l["intrusive"] in indr["ROCKTYPE1"])
                                    ):
                                        rcontact.append(
                                            [
                                                (
                                                    ind,
                                                    lastrcode,
                                                    indr["UNIT_NAME"],
                                                    indr.geometry,
                                                )
                                            ]
                                        )
                                        lastrcode = indr["UNIT_NAME"]

                    locations = [(firstlx, firstly)]
                    first_height_l = m2l_utils.value_from_dtm_dtb(
                        dtm, dtb, dtb_null, cover_map, locations
                    )
                    ostr = "{},{},{},{}\n".format(
                        firstlx, firstly, first_height_l, firstlc
                    )
                    fftc.write(ostr)
                    locations = [(firstrx, firstry)]
                    first_height_r = m2l_utils.value_from_dtm_dtb(
                        dtm, dtb, dtb_null, cover_map, locations
                    )
                    ostr = "{},{},{},{}\n".format(
                        firstrx, firstry, first_height_r, firstrc
                    )
                    fftc.write(ostr)
                    locations = [(lastlx, lastly)]
                    last_height_l = m2l_utils.value_from_dtm_dtb(
                        dtm, dtb, dtb_null, cover_map, locations
                    )
                    ostr = "{},{},{},{}\n".format(lastlx, lastly, last_height_l, lastlc)
                    fftc.write(ostr)
                    locations = [(lastrx, lastry)]
                    last_height_r = m2l_utils.value_from_dtm_dtb(
                        dtm, dtb, dtb_null, cover_map, locations
                    )
                    ostr = "{},{},{},{}\n".format(lastrx, lastry, last_height_r, lastrc)
                    fftc.write(ostr)

                    # loop through left and right sides to find equivalent contact pairs along fault

                    if len(lcontact) > 0 and len(rcontact) > 0:
                        for lc in lcontact:
                            for rc in rcontact:
                                # display('l',lc[0][3].x,'r',rc[0][3].x)
                                if (
                                    lc[0][1] == rc[0][1]
                                    and lc[0][2] == rc[0][2]
                                    and not lc[0][1] == ""
                                ):
                                    dist = m2l_utils.ptsdist(
                                        lc[0][3].x, lc[0][3].y, rc[0][3].x, rc[0][3].y
                                    )
                                    if lc[0][0] < rc[0][0]:
                                        dist = -dist
                                    # print('***',lc,rc)

                                    xi.append((lc[0][3].x))
                                    yi.append((lc[0][3].y))
                                    l, m = m2l_utils.pts2dircos(
                                        lc[0][3].x, lc[0][3].y, rc[0][3].x, rc[0][3].y
                                    )
                                    if not (l == 0.0 and m == 0.0):
                                        fdc.append(
                                            (
                                                l,
                                                m,
                                                "Fault_"
                                                + str(fault["GEOMETRY_OBJECT_ID"]),
                                            )
                                        )
                                        all_coordsdist.append((dist))

    fftc.close()
    structure_file = os.path.join(tmp_path, "structure_clip.shp")

    # first calculate interpolation for fault displacement calcs
    interpolate_orientations(
        structure_file, tmp_path, bbox, c_l, use_gcode, scheme, xi, yi, True
    )
    # then for near-fault calcs
    interpolate_orientations(
        structure_file,
        os.path.join(tmp_path, "ex_"),
        bbox,
        c_l,
        use_gcode,
        scheme,
        all_coords_x,
        all_coords_y,
        True,
    )

    basal_contacts_file = os.path.join(tmp_path, "basal_contacts.shp")

    interpolate_contacts(
        basal_contacts_file,
        tmp_path,
        dtm,
        dtb,
        dtb_null,
        cover_map,
        bbox,
        c_l,
        use_gcode2,
        scheme,
        xi,
        yi,
        True,
    )
    interpolate_contacts(
        basal_contacts_file,
        os.path.join(tmp_path, "ex_"),
        dtm,
        dtb,
        dtb_null,
        cover_map,
        bbox,
        c_l,
        use_gcode2,
        scheme,
        all_coords_x,
        all_coords_y,
        True,
    )

    combo_file = os.path.join(tmp_path, "f_combo.csv")
    ex_combo_file = os.path.join(tmp_path, "ex_f_combo.csv")

    lc = np.loadtxt(
        os.path.join(tmp_path, "f_interpolation_contacts_l.csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )
    mc = np.loadtxt(
        os.path.join(tmp_path, "f_interpolation_contacts_m.csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )
    lo = np.loadtxt(
        os.path.join(tmp_path, "f_interpolation_l.csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )
    mo = np.loadtxt(
        os.path.join(tmp_path, "f_interpolation_m.csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )
    no = np.loadtxt(
        os.path.join(tmp_path, "f_interpolation_n.csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )
    xy = np.loadtxt(
        os.path.join(tmp_path, "f_interpolation_" + scheme + ".csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )

    ex_lc = np.loadtxt(
        os.path.join(tmp_path, "ex_f_interpolation_contacts_l.csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )
    ex_mc = np.loadtxt(
        os.path.join(tmp_path, "ex_f_interpolation_contacts_m.csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )
    ex_lo = np.loadtxt(
        os.path.join(tmp_path, "ex_f_interpolation_l.csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )
    ex_mo = np.loadtxt(
        os.path.join(tmp_path, "ex_f_interpolation_m.csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )
    ex_no = np.loadtxt(
        os.path.join(tmp_path, "ex_f_interpolation_n.csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )
    ex_xy = np.loadtxt(
        os.path.join(tmp_path, "ex_f_interpolation_" + scheme + ".csv"),
        skiprows=1,
        delimiter=",",
        dtype=float,
    )

    join_contacts_and_orientations(
        combo_file,
        geology_file,
        tmp_path,
        dtm_reproj_file,
        dtb,
        dtb_null,
        cover_map,
        c_l,
        lo,
        mo,
        no,
        lc,
        mc,
        xy,
        dst_crs,
        bbox,
        True,
    )
    join_contacts_and_orientations(
        ex_combo_file,
        geology_file,
        os.path.join(tmp_path, "ex_"),
        dtm_reproj_file,
        dtb,
        dtb_null,
        cover_map,
        c_l,
        ex_lo,
        ex_mo,
        ex_no,
        ex_lc,
        ex_mc,
        ex_xy,
        dst_crs,
        bbox,
        True,
    )

    # loop though all nodes and calulate true displacement based on local estimated dip

    ddd = pd.read_csv(os.path.join(tmp_path, "f_combo_full.csv"))
    f = open(os.path.join(output_path, "fault_displacements3.csv"), "w")
    f.write("X,Y,fname,apparent_displacement,vertical_displacement\n")

    for i in range(len(fdc)):
        l, m, n = m2l_utils.ddd2dircos(ddd.iloc[i]["dip"], ddd.iloc[i]["azimuth"])
        lnorm = l / sqrt(pow(l, 2) + pow(m, 2))
        mnorm = m / sqrt(pow(l, 2) + pow(m, 2))
        dotproduct = fabs((fdc[i][0] * lnorm) + (fdc[i][1] * mnorm))
        ostr = (
            str(xi[i])
            + ","
            + str(yi[i])
            + ","
            + str(fdc[i][2])
            + ","
            + str(int(all_coordsdist[i]))
            + ","
            + str(
                abs(
                    int(
                        all_coordsdist[i]
                        * tan(radians(dotproduct * ddd.iloc[i]["dip"]))
                    )
                )
            )
            + "\n"
        )
        f.write(ostr)

    f.close()
    print(
        "fault displacement estimates saved as",
        os.path.join(output_path, "fault_displacements3.csv"),
    )
    print(
        "near-fault orientations saved as",
        os.path.join(tmp_path, "ex_f_combo_full.csv"),
    )


def call_interpolator_grid(calc, x, y, l, m, n, xi, yi):
    # Calculate IDW or other interpolators

    ZIl = interpolator_switch(calc, x, y, l, xi, yi)
    ZIm = interpolator_switch(calc, x, y, m, xi, yi)
    if type(n) is not int:
        ZIn = interpolator_switch(calc, x, y, n, xi, yi)
    else:
        ZIn = 0

    return (ZIl, ZIm, ZIn)


def interpolate_orientation_grid(structures, calc, xcoords, ycoords, c_l):

    npts = len(structures)
    x = np.zeros(npts)
    y = np.zeros(npts)
    dip = np.zeros(npts)
    dipdir = np.zeros(npts)

    i = 0
    l = np.zeros(npts)
    m = np.zeros(npts)
    n = np.zeros(npts)
    for a_pt in structures.iterrows():
        x[i] = a_pt[1]["geometry"].x + (np.random.ranf())
        y[i] = a_pt[1]["geometry"].y + (np.random.ranf())
        dip[i] = a_pt[1]["DIP"]

        # All orientation have been converted to dipdir
        dipdir[i] = a_pt[1]["DIPDIR"]

        i = i + 1

    for i in range(0, npts):
        l[i], m[i], n[i] = m2l_utils.ddd2dircos(dip[i], dipdir[i])
        # this code is now in the right place?
        if structures.iloc[i]["POLARITY"] == c_l["btype"]:
            l[i] = -l[i]
            m[i] = -m[i]
            n[i] = -n[i]

    ZIl, ZIm, ZIn = call_interpolator_grid(calc, x, y, l, m, n, xcoords, ycoords)

    l2 = ZIl / np.sqrt(ZIl**2 + ZIm**2 + ZIn**2)
    m2 = ZIm / np.sqrt(ZIl**2 + ZIm**2 + ZIn**2)
    n2 = ZIn / np.sqrt(ZIl**2 + ZIm**2 + ZIn**2)

    dip = 90.0 - np.degrees(np.arcsin(n2))
    dip_direction = np.where(m2 > 0, (360 + np.degrees(np.arctan(l2 / m2))) % 360, 1)
    dip_direction = np.where(
        m2 < 0, (540 + np.degrees(np.arctan(l2 / m2))) % 360, dip_direction
    )
    dip_direction = np.where(m2 == 0, 90, dip_direction)

    return (l2, m2, n2, dip, dip_direction)


def pts2dircos_arr(p1x, p1y, p2x, p2y):
    dx = p1x - p2x
    dy = p1y - p2y
    return (dx, dy)


def interpolate_contacts_grid(contacts, calc, xcoords_group, ycoords_group):
    decimate = 1
    i = 0
    listarray = []
    for (
        indx,
        acontact,
    ) in contacts.iterrows():  # loop through distinct linestrings in MultiLineString
        if acontact.geometry.type == "MultiLineString":
            for line in acontact.geometry.geoms:
                if i % decimate == 0:
                    listarray.append([line.coords[0][0], line.coords[0][1]])
                i = i + 1
        else:
            if 1 % decimate == 0:
                listarray.append(
                    [acontact.geometry.coords[0][0], acontact.geometry.coords[0][1]]
                )
            i = i + 1

    coords = np.array(listarray)
    if len(coords) & 0x1:
        coords2 = coords[: len(coords) - 1, :].reshape((int(len(coords) / 2), 4))
    else:
        coords2 = coords[: len(coords), :].reshape((int(len(coords) / 2), 4))

    dx = coords2[:, 0:1] - coords2[:, 2:3]
    dy = coords2[:, 1:2] - coords2[:, 3:4]
    x = coords2[:, 0:1] + (dx / 2)
    y = coords2[:, 1:2] + (dy / 2)

    dx2 = dx**2
    dy2 = dy**2
    mask = np.logical_and(dx2 > 0, dy2 > 0)

    dx = dx[mask]
    dy = dy[mask]
    dx2 = dx2[mask]
    dy2 = dy2[mask]
    x = x[mask]
    y = y[mask]
    scale = np.sqrt(dx2 + dy2)
    l = dx / scale
    m = dy / scale
    # m=np.where(l<0, -m, m)
    # l=np.where(l<0, -l, l)

    if len(x) > 2:
        ZIl, ZIm, ZIn = call_interpolator_grid(
            calc, x, y, l, m, 0, xcoords_group, ycoords_group
        )
        l2 = ZIl / np.sqrt(ZIl**2 + ZIm**2)
        m2 = ZIm / np.sqrt(ZIl**2 + ZIm**2)
        S = np.degrees(np.arctan2(l2, m2))

        return (l2, m2, S)
    else:
        return (0, 0, 0)


@beartype.beartype
def interpolation_grids(
    config: Config, map_data, basal_contacts_filename: str, super_groups: list
):  # -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Note: Hint type returns don't work properly for python 3.8 or older
    geology = map_data.get_map_data(Datatype.GEOLOGY).copy()
    orientations = map_data.get_map_data(Datatype.STRUCTURE).copy()

    is_bed = orientations["STRUCTURE_TYPE"].str.contains(
        config.c_l["bedding"], regex=False
    )

    orientations = orientations[is_bed]
    contacts = gpd.read_file(basal_contacts_filename, bbox=config.bbox)

    geology["GROUP"].fillna(geology["GROUP2"], inplace=True)
    geology["GROUP"].fillna(geology["UNIT_NAME"], inplace=True)

    spacing = config.run_flags["interpolation_spacing"]
    scheme = config.run_flags["interpolation_scheme"]
    if spacing < 0:
        spacing = -(config.bbox[2] - config.bbox[0]) / spacing

    # x = (config.bbox[2] - config.bbox[0]) / spacing
    # y = (config.bbox[3] - config.bbox[1]) / spacing

    xcoords = np.arange(config.bbox[0], config.bbox[2], spacing)
    ycoords = np.arange(config.bbox[1], config.bbox[3], spacing)
    xcoords, ycoords = np.meshgrid(xcoords, ycoords)
    bsh = xcoords.shape
    xcoords, ycoords = xcoords.flatten(), ycoords.flatten()

    xycoords = np.vstack((xcoords, ycoords)).transpose()
    xycoords = xycoords.reshape(len(xcoords), 2)

    nodes = gpd.GeoDataFrame(
        xycoords, geometry=[Point(xy) for xy in zip(xcoords, ycoords)]
    )
    nodes.crs = map_data.working_projection
    nodes_code = gpd.sjoin(nodes, geology, how="left", predicate="within")
    # orientations = gpd.sjoin(structures, geology, how="left", predicate="within")
    orientations = orientations[orientations["DIP"] != 0]
    first_supergroup = True
    # avoids massive memory rbf calcs by splitting calc into (non-threaded) chunks, maybe try dask + masks??
    split = 25000
    # print('split',split,'bsh[0]',bsh[0],'bsh[1]',bsh[1],'spacing',spacing,'xdim',config.bbox[2]-config.bbox[0],'ydim',config.bbox[3]-config.bbox[1])
    for i in np.arange(0, bsh[0] * bsh[1], split):
        min = i
        max = i + split
        if max > bsh[0] * bsh[1]:
            max = bsh[0] * bsh[1] + 1
        # print('rbf_split', min, max)
        nodes_code_row = nodes_code[min:max]
        for groups in super_groups:
            if config.verbose_level != VerboseLevel.NONE:
                print(groups)
            first = True
            for group in groups:

                if first:
                    all_nodes = nodes_code_row[nodes_code_row["GROUP"] == group]
                    all_structures = orientations[orientations["GROUP"] == group]
                    all_contacts = contacts[
                        contacts["GROUP"] == group.replace(" ", "_").replace("-", "_")
                    ]
                    first = False
                else:
                    another_node = nodes_code_row[nodes_code_row["GROUP"] == group]
                    all_nodes = pd.concat([all_nodes, another_node], sort=False)

                another_contact = contacts[
                    contacts["GROUP"] == group.replace(" ", "_").replace("-", "_")
                ]
                all_contacts = pd.concat([all_contacts, another_contact], sort=False)

                another_structure = orientations[orientations["GROUP"] == group]
                all_structures = pd.concat(
                    [all_structures, another_structure], sort=False
                )

            xcoords_group = all_nodes.geometry.x
            ycoords_group = all_nodes.geometry.y

            if len(xcoords_group) > 0:
                if len(all_structures) > 2:
                    l, m, n, d, dd = interpolate_orientation_grid(
                        all_structures, scheme, xcoords_group, ycoords_group, config.c_l
                    )
                    xy_lmn = np.vstack(
                        (xcoords_group, ycoords_group, l, m, n, d, dd)
                    ).transpose()
                    xy_lmn = xy_lmn.reshape(len(l), 7)
                else:
                    if config.verbose_level != VerboseLevel.NONE:
                        print(groups, "has no structures")

                    xy_lmn = np.zeros((5, len(xcoords_group)))
                    xy_lmn = np.vstack(
                        (xcoords_group, ycoords_group, xy_lmn)
                    ).transpose()
                    xy_lmn = xy_lmn.reshape(len(xcoords_group), 7)

                if len(all_contacts) > 0:
                    l, m, S = interpolate_contacts_grid(
                        all_contacts, scheme, xcoords_group, ycoords_group
                    )
                    if type(l) is not int:
                        xy_lm_contacts = np.vstack(
                            (xcoords_group, ycoords_group, l, m, S)
                        ).transpose()
                        xy_lm_contacts = xy_lm_contacts.reshape(len(l), 5)
                    else:
                        xy_lm_contacts = np.zeros((3, len(xcoords_group)))
                        xy_lm_contacts = np.vstack(
                            (xcoords_group, ycoords_group, xy_lm_contacts)
                        ).transpose()
                        xy_lm_contacts = xy_lm_contacts.reshape(len(xcoords_group), 5)
                else:
                    if config.verbose_level != VerboseLevel.NONE:
                        print(groups, "has no contacts")

                    xy_lm_contacts = np.zeros((3, len(xcoords_group)))
                    xy_lm_contacts = np.vstack(
                        (xcoords_group, ycoords_group, xy_lm_contacts)
                    ).transpose()
                    xy_lm_contacts = xy_lm_contacts.reshape(len(xcoords_group), 5)

                if first_supergroup:
                    first_supergroup = False
                    xy_lmn_all = np.copy(xy_lmn)
                    xy_lm_contacts_all = np.copy(xy_lm_contacts)
                else:
                    xy_lmn_all = np.vstack((xy_lmn_all, xy_lmn))
                    xy_lm_contacts_all = np.vstack((xy_lm_contacts_all, xy_lm_contacts))

    # sort to get back to x,y grid ordering
    dt = [
        ("X", xy_lmn_all.dtype),
        ("Y", xy_lmn_all.dtype),
        ("l", xy_lmn_all.dtype),
        ("m", xy_lmn_all.dtype),
        ("n", xy_lmn_all.dtype),
        ("dip", xy_lmn_all.dtype),
        ("dip_dir", xy_lmn_all.dtype),
    ]
    # assert xy_lmn_all.flags['C_CONTIGUOUS']
    orientation_interp = xy_lmn_all.ravel().view(dt)
    orientation_interp.sort(order=["X", "Y", "l", "m", "n", "dip", "dip_dir"])

    dt = [
        ("X", xy_lm_contacts_all.dtype),
        ("Y", xy_lm_contacts_all.dtype),
        ("l", xy_lm_contacts_all.dtype),
        ("m", xy_lm_contacts_all.dtype),
        ("angle", xy_lm_contacts_all.dtype),
    ]
    # assert xy_lm_contacts_all.flags['C_CONTIGUOUS']
    contact_interp = xy_lm_contacts_all.ravel().view(dt)
    contact_interp.sort(order=["X", "Y", "l", "m", "angle"])

    scale = np.sqrt(1 - (orientation_interp["n"] ** 2))
    lscaled = -scale * contact_interp["m"]
    mscaled = scale * contact_interp["l"]
    scale2 = np.sqrt(orientation_interp["l"] ** 2 + orientation_interp["m"] ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        loscaled = np.where(scale2 > 0, contact_interp["l"] / scale2, 0)
        moscaled = np.where(scale2 > 0, contact_interp["m"] / scale2, 0)
    dotproduct = (-contact_interp["m"] * loscaled) + (contact_interp["l"] * moscaled)

    # lscaled=np.where(dotproduct<0, -lscaled,lscaled)
    # mscaled=np.where(dotproduct<0, -mscaled,mscaled)
    dotproduct = np.where(dotproduct > 1, 1, dotproduct)
    dotproduct = np.where(dotproduct < -1, -1, dotproduct)
    # misorientation = np.degrees(np.arccos(dotproduct))

    dip = 90.0 - np.degrees(np.arcsin(orientation_interp["n"]))
    mscaled = np.where(mscaled == 0, 1e-5, mscaled)
    # dip_direction=np.where(mscaled>0, (360+np.degrees(np.arctan2(lscaled,mscaled)))%360, 1)
    dip_direction = (360 + np.degrees(np.arctan2(lscaled, mscaled))) % 360
    # dip_direction=np.where(mscaled<0, (540+np.degrees(np.arctan2(lscaled,mscaled)))%360, dip_direction)
    # dip_direction=np.where(mscaled==0, 90, dip_direction)

    combo_interp = np.vstack(
        (
            contact_interp["X"],
            contact_interp["Y"],
            lscaled,
            mscaled,
            orientation_interp["n"],
            dip,
            dip_direction,
        )
    ).transpose()

    orientation_interp = pd.DataFrame(orientation_interp)
    contact_interp = pd.DataFrame(contact_interp)
    combo_interp = pd.DataFrame(combo_interp)
    return (contact_interp, combo_interp)


@beartype.beartype
def process_fault_throw_and_near_faults_from_grid(
    config: Config, map_data, workflow, dip_grid, dip_dir_grid
):

    local_faults = map_data.get_map_data(Datatype.FAULT)
    local_faults = local_faults.dropna(subset=["geometry"])
    geology = map_data.get_map_data(Datatype.GEOLOGY)
    dtm = map_data.get_map_data(Datatype.DTM).open()

    all_long_faults = np.genfromtxt(
        os.path.join(config.output_path, "fault_dimensions.csv"),
        delimiter=",",
        dtype="U100",
    )
    fault_names = all_long_faults[1:, :1]
    m_step = 5  # outstep from fault
    decimate_near = 50
    xi = []
    yi = []
    fdc = []
    all_coordsdist = []
    all_coords_x = []
    all_coords_y = []

    fftc = open(os.path.join(config.output_path, "fault_tip_contacts.csv"), "w")
    fftc.write("X,Y,Z,formation\n")

    # loop through all faults

    # Looping through each fault
    for index, fault in local_faults.iterrows():
        # if(fault["GEOMETRY_OBJECT_ID"]!=1071):
        # continue
        if not str(fault.geometry.type) == "None":
            if fault.geometry.type == "LineString":
                # in is dangerous as Fault_1 is in Fault_10
                if "Fault_" + str(fault["GEOMETRY_OBJECT_ID"]) in fault_names:
                    # print('LineString','Fault_'+str(fault["GEOMETRY_OBJECT_ID"]),len(fault.geometry.coords))
                    lcoords = []
                    rcoords = []
                    index = []

                    # make a list of points just offset from mid-points between fault nodes and join
                    # geology polygon information to points
                    j = 0
                    for i in range(0, len(fault.geometry.coords) - 1):
                        l, m = m2l_utils.pts2dircos(
                            fault.geometry.coords[i][0],
                            fault.geometry.coords[i][1],
                            fault.geometry.coords[i + 1][0],
                            fault.geometry.coords[i + 1][1],
                        )
                        dx = m_step * m
                        dy = m_step * l
                        for inc in np.arange(0.1, 1, 0.01):
                            midx = fault.geometry.coords[i][0] + (
                                (
                                    fault.geometry.coords[i + 1][0]
                                    - fault.geometry.coords[i][0]
                                )
                                * inc
                            )
                            midy = fault.geometry.coords[i][1] + (
                                (
                                    fault.geometry.coords[i + 1][1]
                                    - fault.geometry.coords[i][1]
                                )
                                * inc
                            )
                            lcoords.append([(midx + dx, midy - dy)])
                            rcoords.append([(midx - dx, midy + dy)])
                            index.append([(j)])
                            j = j + 1
                    lgeom = [Point(xy) for xy in lcoords]
                    rgeom = [Point(xy) for xy in rcoords]

                    lgdf = GeoDataFrame(
                        index, crs=map_data.working_projection, geometry=lgeom
                    )
                    rgdf = GeoDataFrame(
                        index, crs=map_data.working_projection, geometry=rgeom
                    )
                    lcode = gpd.sjoin(lgdf, geology, how="left", predicate="within")
                    rcode = gpd.sjoin(rgdf, geology, how="left", predicate="within")
                    # display(lcode)

                    # add lots of points left and right of fault to make sure ellipse range is happy in geomodeller
                    lgroups = []
                    for ind, indl in lcode.iterrows():
                        if (
                            not str(indl["UNIT_NAME"]) == "nan"
                            and not str(indl["ROCKTYPE1"]) == "nan"
                        ):
                            if ind % decimate_near == 0 or ind == len(lcode) - 1:
                                if (not config.c_l["sill"] in indl["DESCRIPTION"]) or (
                                    not config.c_l["intrusive"] in indl["ROCKTYPE1"]
                                ):
                                    locations = [(indl.geometry.x, indl.geometry.y)]
                                    last_height_l = m2l_utils.value_from_dtm_dtb(
                                        dtm,
                                        map_data.dtb,
                                        map_data.dtb_null,
                                        workflow["cover_map"],
                                        locations,
                                    )
                                    if not indl["GROUP"] in lgroups:
                                        ostr = "{},{},{},{}\n".format(
                                            indl.geometry.x,
                                            indl.geometry.y,
                                            last_height_l,
                                            indl["UNIT_NAME"]
                                            .replace(" ", "_")
                                            .replace("-", "_"),
                                        )
                                        fftc.write(ostr)
                                        lgroups.append(indl["GROUP"])
                    rgroups = []
                    for ind, indr in rcode.iterrows():
                        if (
                            not str(indr["UNIT_NAME"]) == "nan"
                            and not str(indr["ROCKTYPE1"]) == "nan"
                        ):
                            if (not config.c_l["sill"] in indr["DESCRIPTION"]) or (
                                not config.c_l["intrusive"] in indr["ROCKTYPE1"]
                            ):
                                if ind % decimate_near == 0 or ind == len(rcode) - 1:
                                    locations = [(indr.geometry.x, indr.geometry.y)]
                                    last_height_r = m2l_utils.value_from_dtm_dtb(
                                        dtm,
                                        map_data.dtb,
                                        map_data.dtb_null,
                                        workflow["cover_map"],
                                        locations,
                                    )
                                    if not indr["GROUP"] in rgroups:
                                        ostr = "{},{},{},{}\n".format(
                                            indr.geometry.x,
                                            indr.geometry.y,
                                            last_height_r,
                                            indr["UNIT_NAME"]
                                            .replace(" ", "_")
                                            .replace("-", "_"),
                                        )
                                        fftc.write(ostr)
                                        rgroups.append(indr["GROUP"])

                    # add points to list if they have different geology code than previous node on left side

                    first = True
                    lcontact = []
                    lastlcode = ""

                    for ind, indl in lcode.iterrows():
                        if ind < len(lcode) and not isnan(indl["index_right"]):
                            ntest1 = str(indl["DESCRIPTION"])
                            ntest2 = str(indl["ROCKTYPE1"])

                            if not ntest1 == "None" and not ntest2 == "None":
                                if ind == 1 or (
                                    not lastlcode == indl["UNIT_NAME"]
                                    and (
                                        (not config.c_l["sill"] in indl["DESCRIPTION"])
                                        or (
                                            not config.c_l["intrusive"]
                                            in indl["ROCKTYPE1"]
                                        )
                                    )
                                ):
                                    all_coords_x.append(indl.geometry.x)
                                    all_coords_y.append(indl.geometry.y)
                                    if first:
                                        first = False
                                        firstlx = indl.geometry.x
                                        firstly = indl.geometry.y
                                        firstlc = (
                                            indl["UNIT_NAME"]
                                            .replace(" ", "_")
                                            .replace("-", "_")
                                        )
                                    lastlx = indl.geometry.x
                                    lastly = indl.geometry.y
                                    lastlc = (
                                        indl["UNIT_NAME"]
                                        .replace(" ", "_")
                                        .replace("-", "_")
                                    )

                                if lastlcode == "" and (
                                    (not config.c_l["sill"] in indl["DESCRIPTION"])
                                    or (
                                        not config.c_l["intrusive"] in indl["ROCKTYPE1"]
                                    )
                                ):
                                    lastlcode = indl["UNIT_NAME"]

                                # print('l',ind,indl["UNIT_NAME"],indl["DESCRIPTION"],indl["ROCKTYPE1"])
                                if (
                                    not ntest1 == "None"
                                    and not ntest2 == "None"
                                    and not str(indl["UNIT_NAME"]) == "nan"
                                ):
                                    if (not indl["UNIT_NAME"] == lastlcode) and (
                                        (not config.c_l["sill"] in indl["DESCRIPTION"])
                                        or (
                                            not config.c_l["intrusive"]
                                            in indl["ROCKTYPE1"]
                                        )
                                    ):
                                        lcontact.append(
                                            [
                                                (
                                                    ind,
                                                    lastlcode,
                                                    indl["UNIT_NAME"],
                                                    indl.geometry,
                                                )
                                            ]
                                        )
                                        lastlcode = indl["UNIT_NAME"]
                                        locations = [(lastlx, lastly)]
                                        last_height_l = m2l_utils.value_from_dtm_dtb(
                                            dtm,
                                            map_data.dtb,
                                            map_data.dtb_null,
                                            workflow["cover_map"],
                                            locations,
                                        )
                                        ostr = "{},{},{},{}\n".format(
                                            lastlx,
                                            lastly,
                                            last_height_l,
                                            indl["UNIT_NAME"]
                                            .replace(" ", "_")
                                            .replace("-", "_"),
                                        )
                                        fftc.write(ostr)

                    # add points to list if they have different geology code than previous node on right side

                    first = True
                    rcontact = []
                    lastrcode = ""
                    for ind, indr in rcode.iterrows():
                        if ind < len(rcode) and not isnan(indr["index_right"]):
                            ntest1 = str(indr["DESCRIPTION"])
                            ntest2 = str(indr["ROCKTYPE1"])

                            if not ntest1 == "None" and not ntest2 == "None":
                                if ind == 1 or (
                                    not lastrcode == indr["UNIT_NAME"]
                                    and (
                                        (not config.c_l["sill"] in indr["DESCRIPTION"])
                                        or (
                                            not config.c_l["intrusive"]
                                            in indr["ROCKTYPE1"]
                                        )
                                    )
                                ):
                                    all_coords_x.append(indr.geometry.x)
                                    all_coords_y.append(indr.geometry.y)
                                    if first:
                                        first = False
                                        firstrx = indr.geometry.x
                                        firstry = indr.geometry.y
                                        firstrc = (
                                            indr["UNIT_NAME"]
                                            .replace(" ", "_")
                                            .replace("-", "_")
                                        )
                                    lastrx = indr.geometry.x
                                    lastry = indr.geometry.y
                                    lastrc = (
                                        indr["UNIT_NAME"]
                                        .replace(" ", "_")
                                        .replace("-", "_")
                                    )

                                if lastrcode == "" and (
                                    (not config.c_l["sill"] in indr["DESCRIPTION"])
                                    or (
                                        not config.c_l["intrusive"] in indr["ROCKTYPE1"]
                                    )
                                ):
                                    lastrcode = indr["UNIT_NAME"]
                                # print('r',ind,indr["UNIT_NAME"],indr["DESCRIPTION"],indr["ROCKTYPE1"])

                                # print(lastrcode,ntest1,ntest2,str(indr["UNIT_NAME"]),indr["UNIT_NAME"],config.c_l['sill'],config.c_l['intrusive'])

                                if (
                                    not ntest1 == "None"
                                    and not ntest2 == "None"
                                    and not str(indr["UNIT_NAME"]) == "nan"
                                ):
                                    if (not indr["UNIT_NAME"] == lastrcode) and (
                                        (not config.c_l["sill"] in indr["DESCRIPTION"])
                                        or (
                                            not config.c_l["intrusive"]
                                            in indr["ROCKTYPE1"]
                                        )
                                    ):
                                        rcontact.append(
                                            [
                                                (
                                                    ind,
                                                    lastrcode,
                                                    indr["UNIT_NAME"],
                                                    indr.geometry,
                                                )
                                            ]
                                        )
                                        lastrcode = indr["UNIT_NAME"]
                                        locations = [(lastrx, lastry)]
                                        last_height_r = m2l_utils.value_from_dtm_dtb(
                                            dtm,
                                            map_data.dtb,
                                            map_data.dtb_null,
                                            workflow["cover_map"],
                                            locations,
                                        )
                                        ostr = "{},{},{},{}\n".format(
                                            lastrx,
                                            lastry,
                                            last_height_r,
                                            indr["UNIT_NAME"]
                                            .replace(" ", "_")
                                            .replace("-", "_"),
                                        )
                                        fftc.write(ostr)

                    # loop through left and right sides to find equivalent contact pairs along fault

                    if len(lcontact) > 0 and len(rcontact) > 0:
                        for lc in lcontact:
                            for rc in rcontact:
                                # display('l',lc[0][3].x,'r',rc[0][3].x)
                                if (
                                    lc[0][1] == rc[0][1]
                                    and lc[0][2] == rc[0][2]
                                    and not lc[0][1] == ""
                                ):
                                    dist = m2l_utils.ptsdist(
                                        lc[0][3].x, lc[0][3].y, rc[0][3].x, rc[0][3].y
                                    )
                                    if lc[0][0] < rc[0][0]:
                                        dist = -dist
                                    # print('***',lc,rc)

                                    xi.append((lc[0][3].x))
                                    yi.append((lc[0][3].y))
                                    l, m = m2l_utils.pts2dircos(
                                        lc[0][3].x, lc[0][3].y, rc[0][3].x, rc[0][3].y
                                    )
                                    if not (l == 0.0 and m == 0.0):
                                        fdc.append(
                                            (
                                                l,
                                                m,
                                                "Fault_"
                                                + str(fault["GEOMETRY_OBJECT_ID"]),
                                            )
                                        )
                                        all_coordsdist.append((dist))
            ##############################################################################################
            # Between these dividers shouldn't be neccessary anymore.
            else:
                # in is dangerous as Fault_1 is in Fault_10
                if "Fault_" + str(fault["GEOMETRY_OBJECT_ID"]) in fault_names:
                    for fls in fault.geometry:
                        fault_ls = LineString(fls)
                        lcoords = []
                        rcoords = []
                        index = []
                        # display("MLS DEBUG",fault.geometry.type)

                        j = 0
                        for i in range(0, len(fault_ls.coords) - 1):
                            for inc in np.arange(0.01, 1, 0.01):
                                midx = fault_ls.coords[i][0] + (
                                    (fault_ls.coords[i + 1][0] - fault_ls.coords[i][0])
                                    * inc
                                )
                                midy = fault_ls.coords[i][1] + (
                                    (fault_ls.coords[i + 1][1] - fault_ls.coords[i][1])
                                    * inc
                                )
                                l, m = m2l_utils.pts2dircos(
                                    fault_ls.coords[i][0],
                                    fault_ls.coords[i][1],
                                    fault_ls.coords[i + 1][0],
                                    fault_ls.coords[i + 1][1],
                                )
                                lcoords.append(
                                    [(midx + (m_step * m), midy - (m_step * l))]
                                )
                                rcoords.append(
                                    [(midx - (m_step * m), midy + (m_step * l))]
                                )
                                index.append([(j)])
                                j = j + 1

                        lgeom = [Point(xy) for xy in lcoords]
                        rgeom = [Point(xy) for xy in rcoords]
                        lgdf = GeoDataFrame(
                            index, crs=map_data.working_projection, geometry=lgeom
                        )
                        rgdf = GeoDataFrame(
                            index, crs=map_data.working_projection, geometry=rgeom
                        )
                        lcode = gpd.sjoin(lgdf, geology, how="left", predicate="within")
                        rcode = gpd.sjoin(rgdf, geology, how="left", predicate="within")

                        # add points to list if they have different geology code than previous node on left side

                        first = True
                        lcontact = []
                        lastlcode = ""
                        for ind, indl in lcode.iterrows():

                            if ind < len(lcode) and not isnan(indl["index_right"]):
                                ntest1 = str(indl["DESCRIPTION"])
                                ntest2 = str(indl["ROCKTYPE1"])

                                if not ntest1 == "None" and not ntest2 == "None":
                                    if ind == 1 or (
                                        not lastlcode == indl["UNIT_NAME"]
                                        and (
                                            (
                                                not config.c_l["sill"]
                                                in indl["DESCRIPTION"]
                                            )
                                            or (
                                                not config.c_l["intrusive"]
                                                in indl["ROCKTYPE1"]
                                            )
                                        )
                                    ):
                                        all_coords_x.append(indl.geometry.x)
                                        all_coords_y.append(indl.geometry.y)
                                        if first:
                                            first = False
                                            firstlx = indl.geometry.x
                                            firstly = indl.geometry.y
                                            firstlc = (
                                                indl["UNIT_NAME"]
                                                .replace(" ", "_")
                                                .replace("-", "_")
                                            )
                                        lastlx = indl.geometry.x
                                        lastly = indl.geometry.y
                                        lastlc = (
                                            indl["UNIT_NAME"]
                                            .replace(" ", "_")
                                            .replace("-", "_")
                                        )

                                    if lastlcode == "" and (
                                        (not config.c_l["sill"] in indl["DESCRIPTION"])
                                        or (
                                            not config.c_l["intrusive"]
                                            in indl["ROCKTYPE1"]
                                        )
                                    ):
                                        lastlcode = indl["UNIT_NAME"]

                                    if (
                                        not ntest1 == "None"
                                        and not ntest2 == "None"
                                        and not str(indl["UNIT_NAME"]) == "nan"
                                    ):
                                        if (not indl["UNIT_NAME"] == lastlcode) and (
                                            (
                                                not config.c_l["sill"]
                                                in indl["DESCRIPTION"]
                                            )
                                            or (
                                                not config.c_l["intrusive"]
                                                in indl["ROCKTYPE1"]
                                            )
                                        ):
                                            lcontact.append(
                                                [
                                                    (
                                                        ind,
                                                        lastlcode,
                                                        indl["UNIT_NAME"],
                                                        indl.geometry,
                                                    )
                                                ]
                                            )
                                            lastlcode = indl["UNIT_NAME"]
                                            locations = [(lastlx, lastly)]
                                            last_height_l = (
                                                m2l_utils.value_from_dtm_dtb(
                                                    dtm,
                                                    map_data.dtb,
                                                    map_data.dtb_null,
                                                    workflow["cover_map"],
                                                    locations,
                                                )
                                            )
                                            ostr = "{},{},{},{}\n".format(
                                                lastlx,
                                                lastly,
                                                last_height_l,
                                                indl["UNIT_NAME"]
                                                .replace(" ", "_")
                                                .replace("-", "_"),
                                            )
                                            fftc.write(ostr)

                        # add points to list if they have different geology code than previous node on right side

                        first = True
                        rcontact = []
                        lastrcode = ""
                        for ind, indr in rcode.iterrows():
                            if ind < len(rcode) and not isnan(indr["index_right"]):
                                ntest1 = str(indr["DESCRIPTION"])
                                ntest2 = str(indr["ROCKTYPE1"])

                                if not ntest1 == "None" and not ntest2 == "None":
                                    if ind == 1 or (
                                        not lastrcode == indr["UNIT_NAME"]
                                        and (
                                            (
                                                not config.c_l["sill"]
                                                in indr["DESCRIPTION"]
                                            )
                                            or (
                                                not config.c_l["intrusive"]
                                                in indr["ROCKTYPE1"]
                                            )
                                        )
                                    ):
                                        all_coords_x.append(indr.geometry.x)
                                        all_coords_y.append(indr.geometry.y)
                                        if first:
                                            first = False
                                            firstrx = indr.geometry.x
                                            firstry = indr.geometry.y
                                            firstrc = (
                                                indr["UNIT_NAME"]
                                                .replace(" ", "_")
                                                .replace("-", "_")
                                            )
                                        lastrx = indr.geometry.x
                                        lastry = indr.geometry.y
                                        lastrc = (
                                            indr["UNIT_NAME"]
                                            .replace(" ", "_")
                                            .replace("-", "_")
                                        )

                                    if lastrcode == "" and (
                                        (not config.c_l["sill"] in indr["DESCRIPTION"])
                                        or (
                                            not config.c_l["intrusive"]
                                            in indr["ROCKTYPE1"]
                                        )
                                    ):
                                        lastrcode = indr["UNIT_NAME"]

                                    if (
                                        not ntest1 == "None"
                                        and not ntest2 == "None"
                                        and not str(indr["UNIT_NAME"]) == "nan"
                                    ):
                                        if (not indr["UNIT_NAME"] == lastrcode) and (
                                            (
                                                not config.c_l["sill"]
                                                in indr["DESCRIPTION"]
                                            )
                                            or (
                                                not config.c_l["intrusive"]
                                                in indr["ROCKTYPE1"]
                                            )
                                        ):
                                            rcontact.append(
                                                [
                                                    (
                                                        ind,
                                                        lastrcode,
                                                        indr["UNIT_NAME"],
                                                        indr.geometry,
                                                    )
                                                ]
                                            )
                                            lastrcode = indr["UNIT_NAME"]
                                            locations = [(lastrx, lastry)]
                                            last_height_r = (
                                                m2l_utils.value_from_dtm_dtb(
                                                    dtm,
                                                    map_data.dtb,
                                                    map_data.dtb_null,
                                                    workflow["cover_map"],
                                                    locations,
                                                )
                                            )
                                            ostr = "{},{},{},{}\n".format(
                                                lastrx,
                                                lastry,
                                                last_height_r,
                                                indr["UNIT_NAME"]
                                                .replace(" ", "_")
                                                .replace("-", "_"),
                                            )
                                            fftc.write(ostr)

                        locations = [(firstlx, firstly)]
                        first_height_l = m2l_utils.value_from_dtm_dtb(
                            dtm,
                            map_data.dtb,
                            map_data.dtb_null,
                            workflow["cover_map"],
                            locations,
                        )
                        ostr = "{},{},{},{}\n".format(
                            firstlx,
                            firstly,
                            first_height_l,
                            firstlc.replace(" ", "_").replace("-", "_"),
                        )
                        fftc.write(ostr)
                        locations = [(firstrx, firstry)]
                        first_height_r = m2l_utils.value_from_dtm_dtb(
                            dtm,
                            map_data.dtb,
                            map_data.dtb_null,
                            workflow["cover_map"],
                            locations,
                        )
                        ostr = "{},{},{},{}\n".format(
                            firstrx,
                            firstry,
                            first_height_r,
                            firstrc.replace(" ", "_").replace("-", "_"),
                        )
                        fftc.write(ostr)
                        locations = [(lastlx, lastly)]
                        last_height_l = m2l_utils.value_from_dtm_dtb(
                            dtm,
                            map_data.dtb,
                            map_data.dtb_null,
                            workflow["cover_map"],
                            locations,
                        )
                        ostr = "{},{},{},{}\n".format(
                            lastlx,
                            lastly,
                            last_height_l,
                            lastlc.replace(" ", "_").replace("-", "_"),
                        )
                        fftc.write(ostr)
                        locations = [(lastrx, lastry)]
                        last_height_r = m2l_utils.value_from_dtm_dtb(
                            dtm,
                            map_data.dtb,
                            map_data.dtb_null,
                            workflow["cover_map"],
                            locations,
                        )
                        ostr = "{},{},{},{}\n".format(
                            lastrx,
                            lastry,
                            last_height_r,
                            lastrc.replace(" ", "_").replace("-", "_"),
                        )
                        fftc.write(ostr)
                        if len(lcode) > 5:
                            locations = [
                                (
                                    lcode.iloc[len(lcode) - 3].geometry.x,
                                    lcode.iloc[len(lcode) - 3].geometry.y,
                                )
                            ]
                            last_height_l = m2l_utils.value_from_dtm_dtb(
                                dtm,
                                map_data.dtb,
                                map_data.dtb_null,
                                workflow["cover_map"],
                                locations,
                            )
                            ostr = "{},{},{},{}\n".format(
                                lcode.iloc[len(lcode) - 3].geometry.x,
                                lcode.iloc[len(lcode) - 3].geometry.y,
                                last_height_l,
                                str(lcode.iloc[len(lcode) - 3]["UNIT_NAME"])
                                .replace(" ", "_")
                                .replace("-", "_"),
                            )
                            if (
                                not str(lcode.iloc[len(lcode) - 3]["UNIT_NAME"])
                                == "nan"
                            ):
                                fftc.write(ostr)
                        if len(rcode) > 5:
                            locations = [
                                (
                                    rcode.iloc[len(rcode) - 3].geometry.x,
                                    rcode.iloc[len(rcode) - 3].geometry.y,
                                )
                            ]
                            last_height_r = m2l_utils.value_from_dtm_dtb(
                                dtm,
                                map_data.dtb,
                                map_data.dtb_null,
                                workflow["cover_map"],
                                locations,
                            )
                            ostr = "{},{},{},{}\n".format(
                                lcode.iloc[len(rcode) - 3].geometry.x,
                                rcode.iloc[len(rcode) - 3].geometry.y,
                                last_height_r,
                                str(rcode.iloc[len(rcode) - 3]["UNIT_NAME"])
                                .replace(" ", "_")
                                .replace("-", "_"),
                            )
                            if (
                                not str(rcode.iloc[len(rcode) - 3]["UNIT_NAME"])
                                == "nan"
                            ):
                                fftc.write(ostr)

                        # loop through left and right sides to find equivalent contact pairs along fault

                        if len(lcontact) > 0 and len(rcontact) > 0:
                            for lc in lcontact:
                                for rc in rcontact:
                                    # display('l',lc[0][3].x,'r',rc[0][3].x)
                                    if (
                                        lc[0][1] == rc[0][1]
                                        and lc[0][2] == rc[0][2]
                                        and not lc[0][1] == ""
                                    ):
                                        dist = m2l_utils.ptsdist(
                                            lc[0][3].x,
                                            lc[0][3].y,
                                            rc[0][3].x,
                                            rc[0][3].y,
                                        )
                                        if lc[0][0] < rc[0][0]:
                                            dist = -dist
                                        # print('***',lc,rc)

                                        xi.append((lc[0][3].x))
                                        yi.append((lc[0][3].y))
                                        l, m = m2l_utils.pts2dircos(
                                            lc[0][3].x,
                                            lc[0][3].y,
                                            rc[0][3].x,
                                            rc[0][3].y,
                                        )
                                        if not (l == 0.0 and m == 0.0):
                                            fdc.append(
                                                (
                                                    l,
                                                    m,
                                                    "Fault_"
                                                    + str(fault["GEOMETRY_OBJECT_ID"]),
                                                )
                                            )
                                            all_coordsdist.append((dist))

                ##############################################################################################

    fftc.close()
    fault_dim = pd.read_csv(
        os.path.join(config.output_path, "fault_dimensions.csv"), sep=","
    )
    fault_dim.set_index("Fault", inplace=True)

    fault_orien = pd.read_csv(
        os.path.join(config.output_path, "fault_orientations.csv"), sep=","
    )
    fault_orien = fault_orien.drop_duplicates(subset=["formation"])
    fault_orien.set_index("formation", inplace=True)
    f = open(os.path.join(config.output_path, "fault_displacements3.csv"), "w")
    f.write("X,Y,fname,apparent_displacement,vertical_displacement,downthrow_dir\n")

    for i in range(len(fdc)):
        r = int((yi[i] - config.bbox[1]) / config.run_flags["interpolation_spacing"])
        c = int((xi[i] - config.bbox[0]) / config.run_flags["interpolation_spacing"])
        if not np.isnan(dip_dir_grid[r, c]) and not np.isnan(dip_grid[r, c]):
            l, m, n = m2l_utils.ddd2dircos(dip_grid[r, c], dip_dir_grid[r, c])
            lnorm = l / sqrt(pow(l, 2) + pow(m, 2))
            mnorm = m / sqrt(pow(l, 2) + pow(m, 2))

            # product of (rotation sign between bed dip direction and fault dip direction)
            # and appararent throw sign gives downthrown side relative to dip direction of fault
            lf, mf, nf = m2l_utils.ddd2dircos(
                fault_orien.loc[fdc[i][2]]["dip"],
                fault_orien.loc[fdc[i][2]]["DipDirection"],
            )
            lfnorm = lf / sqrt(pow(lf, 2) + pow(mf, 2))
            mfnorm = mf / sqrt(pow(lf, 2) + pow(mf, 2))

            dotproduct = fabs((fdc[i][0] * lnorm) + (fdc[i][1] * mnorm))
            crossproduct = asin(lnorm * mfnorm - mnorm * lfnorm)
            if crossproduct * all_coordsdist[i] < 0:
                down_dipdir = fault_orien.loc[fdc[i][2]]["DipDirection"]
            else:
                down_dipdir = (fault_orien.loc[fdc[i][2]]["DipDirection"] + 180.0) % 360
            v_d = all_coordsdist[i] * tan(radians(dotproduct * dip_grid[r, c]))
            if not np.isnan(v_d):
                vert_displacement = abs(
                    int(all_coordsdist[i] * tan(radians(dotproduct * dip_grid[r, c])))
                )
                # cap displacements more than 20% of fault length
                if (
                    vert_displacement
                    > fault_dim.loc[fdc[i][2]]["HorizontalRadius"] / 100.0
                ):
                    if config.verbose_level != VerboseLevel.NONE:
                        print(
                            "Fault",
                            fdc[i][2],
                            "with displacement of",
                            vert_displacement,
                            "capped to",
                            fault_dim.loc[fdc[i][2]]["HorizontalRadius"] / 100.0,
                        )
                    vert_displacement = (
                        fault_dim.loc[fdc[i][2]]["HorizontalRadius"] / 100.0
                    )

                ostr = (
                    str(xi[i])
                    + ","
                    + str(yi[i])
                    + ","
                    + str(fdc[i][2])
                    + ","
                    + str(int(all_coordsdist[i]))
                    + ","
                    + str(vert_displacement)
                    + ","
                    + str(down_dipdir)
                    + "\n"
                )
        f.write(ostr)

    f.close()
    if config.verbose_level != VerboseLevel.NONE:
        print(
            "fault displacement estimates saved as",
            os.path.join(config.output_path, "fault_displacements3.csv"),
        )
        print(
            "near-fault orientations saved as",
            os.path.join(config.tmp_path, "ex_f_combo_full.csv"),
        )
        print(
            "near-fault orientations saved as",
            os.path.join(config.tmp_path, "ex_f_combo_full.csv"),
        )

    # when no fault displacment data are available, set to 1m displacment so they still are calculated

    if local_faults is not None:
        fault_ori = pd.read_csv(
            os.path.join(config.output_path, "fault_orientations.csv")
        )
        fault_disp = pd.read_csv(
            os.path.join(config.output_path, "fault_displacements3.csv")
        )

        fault_disp_found = fault_disp["fname"].unique()

        f = open(os.path.join(config.output_path, "fault_displacements3.csv"), "a+")

        for ind, fault in fault_ori.iterrows():
            if not fault["formation"] in fault_disp_found:
                ostr = (
                    str(fault["X"])
                    + ","
                    + str(fault["Y"])
                    + ","
                    + str(fault["formation"])
                    + ","
                    + str(1)
                    + ","
                    + str(1)
                    + ","
                    + str(1)
                    + "\n"
                )
                f.write(ostr)

        f.close()
