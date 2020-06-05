from maps.config import *
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon

from maps.roi import ROI


class Project(object):
    def __init__(self,
                 workflow={'model_engine': 'geomodeller'},
                 # region of interest coordinates in metre-based system (or non-degree system)
                 # (minx, miny, maxx, maxy, bottom, top)
                 bbox=(500057, 7455348, 603028, 7567953, -8200, 1200)
                 ):
        self.update_workflow(workflow)
        self.update_roi(bbox)

    # TODO: Make workflow dictionary more editable

    def update_workflow(self, workflow):
        if(workflow['model_engine'] == 'geomodeller'):
            workflow.update({'seismic_section': False,
                             'cover_map': False,
                             'near_fault_interpolations': True,
                             'fold_axial_traces': False,
                             'stereonets': False,
                             'formation_thickness': True,
                             'polarity': False,
                             'strat_offset': False,
                             'contact_dips': True})
        elif(workflow['model_engine'] == 'loopstructural'):
            workflow.update({'seismic_section': False,
                             'cover_map': False,
                             'near_fault_interpolations': False,
                             'fold_axial_traces': True,
                             'stereonets': False,
                             'formation_thickness': True,
                             'polarity': False,
                             'strat_offset': False,
                             'contact_dips': False})
        elif(workflow['model_engine'] == 'gempy'):
            workflow.update({'seismic_section': False,
                             'cover_map': False,
                             'near_fault_interpolations': False,
                             'fold_axial_traces': True,
                             'stereonets': False,
                             'formation_thickness': False,
                             'polarity': False,
                             'strat_offset': False,
                             'contact_dips': False})
        elif(workflow['model_engine'] == 'noddy'):
            workflow.update({'seismic_section': False,
                             'cover_map': False,
                             'near_fault_interpolations': False,
                             'fold_axial_traces': False,
                             'stereonets': False,
                             'formation_thickness': False,
                             'polarity': False,
                             'strat_offset': False,
                             'contact_dips': False})
        else:
            workflow.update({'seismic_section': False,
                             'cover_map': False,
                             'near_fault_interpolations': False,
                             'fold_axial_traces': False,
                             'stereonets': True,
                             'formation_thickness': True,
                             'polarity': False,
                             'strat_offset': True,
                             'contact_dips': False})
        self.workflow = workflow

    def update_roi(self, bbox):
        # TODO: Make crs defaults and specifiable not from config
        self.bbox = bbox
        minx = bbox[0]
        maxx = bbox[1]
        miny = bbox[2]
        maxy = bbox[3]
        lat_point_list = [miny, miny, maxy, maxy, maxy]
        lon_point_list = [minx, maxx, maxx, minx, minx]
        bbox_geom = Polygon(zip(lon_point_list, lat_point_list))
        polygon = gpd.GeoDataFrame(
            index=[0], crs=dst_crs, geometry=[bbox_geom])
        self.roi = ROI(polygon, self.bbox[:4])


# def main():
#     proj = Project()
#     # print(proj.roi)
#     # print(proj.workflow)
#     print(proj)


# if __name__ == "__main__":
#     main()
