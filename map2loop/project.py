from map2loop.cfg import *
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon

from map2loop.config import Config


class Project(object):
    def __init__(self,
                 # region of interest coordinates in metre-based system (or non-degree system)
                 # (minx, miny, maxx, maxy, bottom, top)
                 # TODO: Remove harcoding if local files ship with examples
                 geology_file,
                 fault_file,
                 structure_file,
                 mindep_file,
                 workflow={'model_engine': 'geomodeller'},
                 ):
        # TODO: Create ways to get and set local files
        self.update_workflow(workflow)
        self.geology_file = geology_file
        self.fault_file = fault_file
        self.structure_file = structure_file
        self.mindep_file = mindep_file

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
                      src_crs={'init': 'EPSG:4326'},
                      dst_crs={'init': 'EPSG:28350'}
                      ):
        # TODO: Make crs defaults and specifiable not from config
        minx, miny, maxx, maxy = tuple([bbox_3d["minx"], bbox_3d["miny"],
                                        bbox_3d["maxx"], bbox_3d["maxy"]])
        lat_point_list = [miny, miny, maxy, maxy, maxy]
        lon_point_list = [minx, maxx, maxx, minx, minx]
        bbox_geom = Polygon(zip(lon_point_list, lat_point_list))
        polygon = gpd.GeoDataFrame(
            index=[0], crs=dst_crs, geometry=[bbox_geom])
        self.config = Config(
            self.geology_file, self.fault_file,
            self.structure_file, self.mindep_file,
            bbox_3d, polygon, step_out, src_crs, dst_crs, c_l
        )

        # Store important data frames and display
        self.config.preprocess("plot")

    def run(self):
        print("Generating topology analyser input...")
        self.config.export_csv()
        self.config.runMap2Model()
        self.config.load_dtm()
        self.config.join_features()

        if(workflow['cover_map']):
            # Depth to basement calculation
            pass
