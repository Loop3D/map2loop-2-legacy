from scipy.interpolate import RegularGridInterpolator
import numpy as np
import geopandas
from shapely.geometry import Point

class MapUtil:
    
    def __init__(self,bounding_box,geology=None,dtm=None):
        """Wrapper for evaluating map data on xy points
        """
        self.bounding_box = bounding_box
        self.geology = geology
        self.dtm = dtm

    def evaluate_dtm_at_points(self,xy,nodataval=np.nan):
        """interpolate the dtm on new locations
        """
        dtm_val = self.dtm.read(1)
        x=np.linspace(self.dtm.bounds[0],self.dtm.bounds[2],dtm_val.shape[1])
        y=np.linspace(self.dtm.bounds[1],self.dtm.bounds[3],dtm_val.shape[0])

        dtm_interpolator = RegularGridInterpolator((x,y),dtm_val.T)
#         inside = self._is_inside(xy)
        interpolated_dtm =np.zeros(xy.shape[0])
        interpolated_dtm[:] = nodataval
        interpolated_dtm[:] = dtm_interpolator(xy[:,:])
        return interpolated_dtm
    def evaluate_geology_at_points(self,xy,column='colour_index'):
        """Extract a numerical column from a shape file, default
        is colour index
        """
        nodes = geopandas.GeoDataFrame(
            xy[:,:1], geometry=[Point(xy) for xy in zip(xy[:,0], xy[:,1])])
        points_geology = geopandas.sjoin(nodes, self.geology, how="left", op="within")
        return points_geology[column].to_numpy()
    
    def _is_inside(self,xy):
        """Check whether inside bounding box, not used
        """
        inside = np.zeros(xy.shape[0],dtype=bool)
        inside = np.logical_and(inside,xy[:,0]>self.bounding_box['minx'])
        inside = np.logical_and(inside,xy[:,0]<self.bounding_box['maxx'])
        inside = np.logical_and(inside,xy[:,1]>self.bounding_box['miny'])
        inside = np.logical_and(inside,xy[:,1]<self.bounding_box['maxy'])
        return inside