import os

# ROI
# # padding arounf dtm to ensure reprojected dtm covers target area (in degrees)
# step_out = 0.1
# inset = 0  # unused??

# # region of interest coordinates in metre-based system (or non-degree system)
# minx = 500057
# maxx = 603028
# miny = 7455348
# maxy = 7567953
# model_top = 1200
# model_base = -8200

# PATHS

# flag to use local or WFS source for data inputs (True = local)
# local_paths = True

# test_data_path = './test_data/'

# if(not os.path.isdir(test_data_path)):
#     os.mkdir(test_data_path)

geology_file = 'hams2_geol.shp'  # input geology file (if local)
# input fault file (if local)
fault_file = 'GEOS_GEOLOGY_LINEARSTRUCTURE_500K_GSD.shp'
# input bedding orientation file (if local)
structure_file = 'hams2_structure.shp'
mindep_file = 'mindeps_2018.shp'  # input mineral deposit file (if local)

# CRS

# coordinate reference system for imported dtms (geodetic lat/long WGS84)
src_crs = {'init': 'EPSG:4326'}
dst_crs = {'init': 'EPSG:28350'}  # coordinate system for data

# CODES AND LABELS
# these refer to specific fields (codes) in GIS layer or database that contain the info needed for these calcs and text substrings (labels) in the contents of these fields
c_l = {
    # Orientations
    "d": "DIP",  # field that contains dip information
    "dd": "DIP_DIR",  # field that contains dip direction information
    "sf": 'FEATURE',  # field that contains information on type of structure
    # text to search for in field defined by sf code to show that this is a bedding measurement
    "bedding": 'Bed',
    # flag to determine measurement convention (currently 'strike' or 'dip direction')
    "otype": 'dip direction',
    "bo": "TYPE",  # field that contains type of foliation
    # text to search for in field defined by bo code to show that this is an overturned bedding measurement
    "btype": 'overturned',
    # Stratigraphy
    "g": 'GROUP_',  # field that contains coarser stratigraphic coding
    # field that contains alternate coarser stratigraphic coding if 'g' is blank
    "g2": 'SUPERSUITE',
    "c": 'UNITNAME',  # field that contains finer stratigraphic coding
    "ds": 'DESCRIPTN',  # field that contains information about lithology
    # field that contains alternate stratigraphic coding (not used??)
    "u": 'CODE',
    "r1": 'ROCKTYPE1',  # field that contains  extra lithology information
    "r2": 'ROCKTYPE2',  # field that contains even more lithology information
    "sill": 'sill',  # text to search for in field defined by ds code to show that this is a sill
    # text to search for in field defined by r1 code to show that this is an intrusion
    "intrusive": 'intrusive',
    # text to search for in field defined by ds code to show that this is an volv=canic (not intrusion)
    "volcanic": 'volcanic',
    # Mineral Deposits
    "msc": 'SITE_CODE',  # field that contains site code of deposit
    "msn": 'SHORT_NAME',  # field that contains short name of deposit
    "mst": 'SITE_TYPE_',  # field that contains site type of deposit
    "mtc": 'TARGET_COM',  # field that contains target commodity of deposit
    "mscm": 'SITE_COMMO',  # field that contains site commodity of deposit
    "mcom": 'COMMODITY_',  # field that contains commodity group of deposit
    # text to search for in field defined by mst code that shows site to ignore
    "minf": 'Infrastructure',
    # Timing
    "min": 'MIN_AGE_MA',  # field that contains minimum age of unit defined by ccode
    "max": 'MAX_AGE_MA',  # field that contains maximum age of unit defined by ccode
    #faults and folds
    "f": 'FEATURE',  # field that contains information on type of structure
    # text to search for in field defined by f code to show that this is a fault
    "fault": 'Fault',
    "ff": 'FEATURE',  # field that contains information on type of structure
    # text to search for in field defined by f code to show that this is a fold axial trace
    "fold": 'Fold axial trace',
    "fdip": 'DIP',               # field for numeric fault dip value
    # text to search for in field defined by fdip to show that this has no known dip
    "fdipnull": '0',
    "fdipdir": 'DIP_DIR',        # field for text fault dip direction value
    # flag for text fault dip direction type num e.g. 045 or alpha e.g. southeast
    "fdipdir_flag": 'alpha',
    "fdipest": 'DIP_EST',        # field for text fault dip estimate value
    # text to search for in field defined by fdipest to give fault dip estimate in increasing steepness
    "fdipest_vals": 'gentle,moderate,steep',
    # field that contains information on name of fault (not used??)
    "n": 'NAME',
    "t": 'TYPE',  # field that contains information on type of fold
    # text to search for in field defined by t to show that this is a syncline
    "syn": 'syncline',
    # ids
    "o": 'OBJECTID',  # field that contains unique id of geometry object
    "gi": 'GEOPNT_ID'  # field that contains unique id of structure point
}

# DECIMATION

# store every nth orientation (in object order) 0 = save all
orientation_decimate = 0
# store every nth contact point (in object order) 0 = save all
contact_decimate = 10
# store every nth fault point (in object order) 0 = save all
fault_decimate = 5
# store every nth fold axial trace point (in object order) 0 = save all
fold_decimate = 5

# INTERPOLATION

gridx = 50  # x grid dimensions (no of points, not distance) for interpolations
gridy = 50  # x grid dimensions (no of points, not distance) for interpolations
scheme = 'scipy_rbf'  # interpolation scheme
# buffer distance for clipping points by faults (in metres or same units as dst_crs)
dist_buffer = 5
intrusion_mode = 0        # 1 all intrusions exluded from basal contacts, 0 only sills
use_interpolations = False    # flag to use interpolated orientations or not.

# ASSUMPTIONS

pluton_dip = 45  # surface dip of pluton contacts
# saucers: \__+_+__/  batholith: +/       \+   domes: /  +  + \  pendant: +\_____/+
pluton_form = 'domes'
fault_dip = 90  # surface dip of faults
