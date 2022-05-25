import geopandas as gpd
import numpy as np
import pandas as pd

def extract_basal_contacts(geology_polygons,column_names):
    """Create a geodataframe with all of the contact lines
    between features in the geology_polygons dataset
    """
    # break multipart features into separate features
    all_geom = geology_polygons.explode(ignore_index=True).reset_index()
    # A feature has an exterior line and can have any number of interior lines
    # we need to extract all interior lines into a geoseries
    indexes = []
    data = []
    for i in all_geom.interiors.index:
        for geom in all_geom.interiors.loc[i]:
            indexes.append(i)
            data.append(geom)
    series = gpd.GeoSeries(data)
    # we want to create a final dataframe with all of the lines but we want to maintain the
    # link to the original features, so copy the indexes from the original feature into the 
    # new dataframe
    interior = gpd.GeoDataFrame(geometry=series)
    interior['indexes'] = indexes
    lines = gpd.GeoDataFrame(geometry=all_geom.exterior)
    lines['indexes'] = all_geom.index
    # merge the interior lines into the exterior lines
    lines = pd.concat([lines,interior]).reset_index()
    # now we want to join the lines feature to the columns of the original shapefile we
    # are interested in
    lines.set_index('indexes').join(all_geom[column_names])
    # now change the index back to the original index
    lines = lines.set_index('index')
    # calculate the intersection of the lines to determine
    # which lines are in contact with each other
    intersection_matrix = np.zeros((len(lines),len(lines)),dtype=bool)
    for i in range(intersection_matrix.shape[0]):
        intersection_matrix[i,:] = lines.intersects(lines.loc[i,'geometry'])
    # now we want to create a dataframe with the lines that are in contact
    linestrings = []
    data = []
    for i in range(intersection_matrix.shape[0]):
        intersection_matrix[i,i]=False
        for j in lines.loc[intersection_matrix[i,:]].index:
            #if lines[i,['age']]
            l = lines.loc[j,'geometry']
            linestrings.append(lines.loc[i,:].geometry.intersection(
                                            l)) 
            data.append([])

    basal_contacts = gpd.gpd.GeoSeries(linestrings)
    return basal_contacts