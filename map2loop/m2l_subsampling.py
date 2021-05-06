import geopandas as gpd
import pandas as pd
import shapely
import matplotlib.pyplot as plt
import sys
import copy
import fiona
import shapely
from shapely.geometry import shape, mapping, LineString, Polygon, MultiLineString, MultiPolygon, MultiPoint
from shapely.geometry.polygon import LinearRing
from shapely.ops import unary_union
import heapq
from itertools import chain

def check_overlaps (gpdataframe):
    """
    Check for overlaps between the features (polygons) within the same geodataframe.

    Args:
        gpdataframe (geopandas dataframe): geopandas dataframe containing features/polygons
        
    Returns:
        geopandas dataframe: geodataframe containing the overlaps
    """
    
    data_temp=gpdataframe
    data_temp['id'] = data_temp.index
    data_overlaps=gpd.GeoDataFrame(crs=data_temp.crs)
    for index, row in data_temp.iterrows():
        data_temp1=data_temp.loc[data_temp.id!=row.id,]
        # check if intersection occured
        overlaps=data_temp1[data_temp1.geometry.overlaps(row.geometry)]['id'].tolist()
        if len(overlaps)>0:
            temp_list=[]
            # compare the area with threshold
            for y in overlaps:
                temp_area=gpd.overlay(data_temp.loc[data_temp.id==y,],data_temp.loc[data_temp.id==row.id,],how='intersection')
                temp_area=temp_area.loc[temp_area.geometry.area>=9e-9]
                if temp_area.shape[0]>0:
                    data_overlaps=gpd.GeoDataFrame(pd.concat([temp_area,data_overlaps],ignore_index=True),crs=data_temp.crs)
    # get unique of list id
    data_overlaps['sorted']=data_overlaps.apply(lambda y: sorted([y['id_1'],y['id_2']]),axis=1)
    data_overlaps['sorted']=data_overlaps.sorted.apply(lambda y: ''.join(str(y)))
    data_overlaps=data_overlaps.drop_duplicates('sorted')
    data_overlaps=data_overlaps.reset_index()[['id_1','id_2','geometry']]
    
    return (data_overlaps)

def check_gaps (gpdataframe):
    """
    Check for gaps between the features (polygons) within the same geodataframe.

    Args:
        gpdataframe (geopandas dataframe): geopandas dataframe containing features/polygons

    Returns:
        geopandas dataframe: geodataframe containing the gaps
    """
    
    data_temp['diss_id']=100
    data_temp_diss=data_temp.dissolve(by='diss_id')
    interior=data_temp_diss.interiors.values.tolist()[0]
    gap_list=[]
    for i in interior:
        gap_list.append(LineString(i))
    data_gaps=gpd.GeoDataFrame(geometry=gap_list,crs=data_temp.crs)
    data_gaps['feature_touches']=data_gaps.geometry.apply(lambda y: data_temp.loc[data_temp.touches(y)]['id'].tolist())
    
    return (data_gaps)