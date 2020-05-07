## To build map2model execute:

mkdir build  
cd build  
cmake ..  
make  

## To run:

./map2model ../Parfile  

## Visualization of output

To visalize output images (.svg file) use e.g. a util 'display' (Linux), or Firefox.  
To visalize output graphs (.gml files) use e.g. yEd Graph Editor (cross-platform).  
In eEd program, to see a graph, choose e.g. Layout --> Organic.  

## Input data format convertion

To convert shape files data to WKT data in UTM coordinates, use ogr2ogr util. Run script shp2csv.sh  
To convert a single point in lat/long coordinates to UTM use gdaltransform util as:  
gdaltransform -s_srs EPSG:4326 -t_srs EPSG:32750  
and then type the point coordinates separated by a space, and press Enter.  

## Outputs description

An output image 'polygons_read.svg' shows all input polygons that were read, with marked clipping window.  
An output image 'polygons_clipped.svg' shows the result of clipping of the input polygons with the clipping window.  
An output graph 'graph_all.gml' shows topology for units, i.e., nodes correspond to units, and edges to unit contacts.  
An output graph 'graph_no_sills.gml' shows same info as 'graph_all.gml', but contacts with units of type 'sills' are excluded.  




