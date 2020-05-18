#!/bin/sh

# Define paths to the shapefiles here.
PATH_IN_GEOLOGY="../data/test_wa_wgs84/T1_500k_geol.shp"
PATH_IN_FAULTS="../data/test_wa_wgs84/T1_500k_linear.shp"
PATH_IN_POINTS="../data/test_wa_wgs84/warox_wgs84.shp"

#PATH_IN_GEOLOGY="../data/full_state_wa/500k_geol.shp"
#PATH_IN_FAULTS="../data/full_state_wa/500k_faults.shp"

# Define paths to the output CSV files here. (They should point to different locations.)
PATH_OUT_GEOLOGY="../data/out_geology"
PATH_OUT_FAULTS="../data/out_faults"
PATH_OUT_POINTS="../data/out_points"

# Remove old data.
rm -r $PATH_OUT_GEOLOGY
rm -r $PATH_OUT_FAULTS
rm -r $PATH_OUT_POINTS

# Set here format parameters for conversion.
# (EPSG:4326 is the srs code for lat/long; and EPSG:32750 is the srs code for UTM 50S)
FORMAT="-lco GEOMETRY=AS_WKT -lco SEPARATOR=TAB -s_srs EPSG:4326 -t_srs EPSG:32750 -f CSV"

#NOTE: if there is an error about PROJ.4 library then try:
#sudo ln -s /usr/lib/libproj.so.0 /usr/lib/libproj.so

# Geology data.
FIELDS="OBJECTID,UNITNAME,GROUP_,MIN_AGE_MA,MAX_AGE_MA,CODE,DESCRIPTN,ROCKTYPE1,ROCKTYPE2" 
ogr2ogr -select $FIELDS $FORMAT $PATH_OUT_GEOLOGY $PATH_IN_GEOLOGY

# Faults data.
FIELDS="OBJECTID,FEATURE"
ogr2ogr -select $FIELDS $FORMAT $PATH_OUT_FAULTS $PATH_IN_FAULTS

# Points data.
FIELDS="GEOPNT_ID,DIP,DIP_DIR"
ogr2ogr -select $FIELDS $FORMAT $PATH_OUT_POINTS $PATH_IN_POINTS
