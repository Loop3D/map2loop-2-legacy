from maps.cfg import *
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon

from maps.config import Config


class Project(object):
    def __init__(self,
                 workflow={'model_engine': 'geomodeller'},
                 # region of interest coordinates in metre-based system (or non-degree system)
                 # (minx, miny, maxx, maxy, bottom, top)
                 # TODO: Remove harcoding if local files ship with examples
                 geology_file="/map2loop-2/maps/data/hams-test.shp",
                 fault_file="/map2loop-2/maps/data/GEOS_GEOLOGY_LINEARSTRUCTURE_500K_GSD.shp",
                 structure_file="/map2loop-2/maps/data/hams2_structure.shp",
                 mindep_file="/map2loop-2/maps/data/mindeps_2018.shp",
                 src_crs={'init': 'EPSG:4326'},
                 dst_crs={'init': 'EPSG:28350'}
                 ):
        # TODO: Create ways to get and set local files
        self.update_workflow(workflow)
        self.geology_file = geology_file
        self.fault_file = fault_file
        self.structure_file = structure_file
        self.mindep_file = mindep_file
        self.src_crs = src_crs
        self.dst_crs = dst_crs

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

    def update_config(self,
                      bbox_3d={
                          "minx": 500057,
                          "maxx": 7455348,
                          "miny": 603028,
                          "maxy": 7567953,
                          "base": -8200,
                          "top": 1200,
                      },
                      step_out=0.1,
                      ):
        # TODO: Make crs defaults and specifiable not from config
        minx, miny, maxx, maxy = list(bbox_3d.values())[:4]
        lat_point_list = [miny, miny, maxy, maxy, maxy]
        lon_point_list = [minx, maxx, maxx, minx, minx]
        bbox_geom = Polygon(zip(lon_point_list, lat_point_list))
        polygon = gpd.GeoDataFrame(
            index=[0], crs=self.dst_crs, geometry=[bbox_geom])
        self.config = Config(self.geology_file, self.fault_file,
                             self.structure_file, self.mindep_file, bbox_3d, polygon, step_out, c_l)

    def plot(self):
        self.config.preprocess("plot")
