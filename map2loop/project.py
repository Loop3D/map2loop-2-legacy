import os
import sys
import logging
import urllib.request
import warnings
import hjson
from tqdm import tqdm


import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon

from map2loop.config import Config
from map2loop.m2l_utils import display, enable_quiet_mode, disable_quiet_mode, print


class Project(object):
    """A high level object implementation of the map2loop workflow.
    """

    def __init__(self,
                 # region of interest coordinates in metre-based system (or non-degree system)
                 # (minx, miny, maxx, maxy, bottom, top)
                 # TODO: Remove harcoding if local files ship with examples
                 state,
                 remote=False,
                 geology_file=None,
                 fault_file=None,
                 fold_file=None,
                 structure_file=None,
                 mindep_file=None,
                 # TODO: Abstract metadata away from user when using loop remote files
                 metadata=None,
                 workflow={'model_engine': 'loopstructural'},
                 ):

        warnings.filterwarnings('ignore')

        if state in ['WA', 'NSW', 'VIC', 'SA', 'QLD', 'ACT', 'TAS']:
            self.state = state
        else:
            print(
                "Please provide the a string of the state code.")
            sys.exit(1)

        # If local data files are expected and not given throw an error
        self.local = not remote
        if remote == False:
            if any(source is None for source in [geology_file, fault_file, fold_file, structure_file, mindep_file]):
                sys.exit(
                    "Please pass local file paths as input params if you do not wish to use remote sources.")

        # If remote, these will be set to null for now, otherwise set to local paths
        self.geology_file = geology_file
        self.fault_file = fault_file
        self.fold_file = fold_file
        self.structure_file = structure_file
        self.mindep_file = mindep_file

        meta_error = "When using your own local files or remote data, pass the path or url to a valid metadata file (.json or .hjson) that describes your input column names.\n"
        meta_error += "You can find an example config here https://gist.github.com/yohanderose/a127c29cb88529f049a5bafc881bb1a0"
        # Load in dictionary that describes the column names
        # TODO: Ensure hjson works with regular json files too
        if (metadata is None) and (geology_file is not None):
            # Check metadata is provided if non loop sources are given
            sys.exit(meta_error)
        elif remote:
            # Pass if using loop sources
            pass
        else:
            # Try to read metadata if given
            self.read_metadata(metadata)

        self.set_proj_defaults()
        self.update_workflow(workflow)

    def set_proj_defaults(self):
        # If local, set the maximum bounding box to the bounds of the input files
        # If remote, set bounds to zeros, hopefully this gives the whole map...
        if self.geology_file is not None:
            print("not remote")
            geology = gpd.read_file(self.geology_file)
            self.proj_bounds = geology.total_bounds
            self.proj_crs = geology.crs
        else:
            # TODO: doesn't like zeros
            self.proj_bounds = (0, 0, 0, 0)
            # TODO: If remote, make epsg required
            self.proj_crs = {'init': 'EPSG:28350'}

        self.step_out = 0.1

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
                             'strat_offset': True,
                             'contact_dips': True})
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
                      out_dir,
                      overwrite=False,
                      bbox_3d={
                          "minx": 0,
                          "maxx": 0,
                          "maxx": 0,
                          "maxy": 0,
                          "base": -10000,
                          "top": 1200,
                      },
                      dtm_crs={'init': 'EPSG:4326'},
                      proj_crs=None,
                      step_out=None,
                      quiet=False,
                      **kwargs
                      ):

        if bbox_3d["minx"] == 0 and bbox_3d["maxx"] == 0:
            bbox_3d.update({
                "minx": self.proj_bounds[0],
                "minx": self.proj_bounds[1],
                "minx": self.proj_bounds[2],
                "minx": self.proj_bounds[3],
            })

        if proj_crs is None:
            proj_crs = self.proj_crs

        if step_out is None:
            step_out = self.step_out

        self.quiet = quiet

        bbox = tuple([bbox_3d["minx"], bbox_3d["miny"],
                      bbox_3d["maxx"], bbox_3d["maxy"]])
        minx, miny, maxx, maxy = bbox
        lat_point_list = [miny, miny, maxy, maxy, maxy]
        lon_point_list = [minx, maxx, maxx, minx, minx]
        bbox_geom = Polygon(zip(lon_point_list, lat_point_list))
        polygon = gpd.GeoDataFrame(
            index=[0], crs=proj_crs, geometry=[bbox_geom])

        # Define the url queries if remote flag is set
        if self.geology_file is None:
            self.fetch_sources(bbox)

        self.config = Config(out_dir, overwrite,
                             self.geology_file, self.fault_file, self.fold_file,
                             self.structure_file, self.mindep_file,
                             bbox_3d, polygon, step_out, dtm_crs, proj_crs, self.local, self.quiet, self.c_l,
                             **kwargs
                             )

        self.config.preprocess()

    def read_metadata(self, filename):
        try:
            if filename.startswith("http"):
                with urllib.request.urlopen(filename) as raw_data:
                    self.c_l = hjson.load(raw_data)
            else:
                with open(filename) as raw_data:
                    self.c_l = hjson.load(raw_data)
        except Exception as e:
            print(e)
            sys.exit("ERROR: Unable to read provided filename file")

    def fetch_sources(self, bbox):
        if self.state == "WA":
            bbox_str = "{},{},{},{}".format(bbox[0], bbox[1], bbox[2], bbox[3])
            self.geology_file = 'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=loop:geol_500k&bbox={}&srs=EPSG:28350'.format(
                bbox_str)
            self.fault_file = 'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=linear_500k&bbox={}&srs=EPSG:28350'.format(
                bbox_str)
            self.fold_file = 'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=linear_500k&bbox={}&srs=EPSG:28350'.format(
                bbox_str)
            self.structure_file = 'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=waroxi_wa_28350_bed&bbox={}&srs=EPSG:28350'.format(
                bbox_str)
            self.mindep_file = 'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=loop:mindeps_2018_28350&bbox={}&srs=EPSG:28350'.format(
                bbox_str)
            self.metadata = 'https://gist.githubusercontent.com/yohanderose/8f843de0dde531f009a3973cbdadcc9f/raw/918f412ae488ce1a6bca188306f7730061ecf551/meta_remote.hjson'
            self.read_metadata(self.metadata)
        if self.state == "QLD":
            bbox_str = "{},{},{},{}".format(bbox[0], bbox[1], bbox[2], bbox[3])
            self.structure_file = 'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=outcrops_28355&bbox={}&srsName=EPSG:28355'.format(
                bbox_str)
            self.geology_file = 'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=qld_geol_dissolved_join_fix_mii_clip_wgs84&bbox={}&srsName=EPSG:28355'.format(
                bbox_str)
            self.mindep_file = 'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=qld_mindeps_28355&bbox={}&srsName=EPSG:28355'.format(
                bbox_str)
            self.fault_file = 'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=qld_faults_folds_28355&bbox={}&srsName=EPSG:28355'.format(
                bbox_str)
            self.fold_file = 'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=qld_faults_folds_28355&bbox={}&srsName=EPSG:28355'.format(
                bbox_str)
            self.metadata = "https://gist.githubusercontent.com/yohanderose/7ad2ae1b36e4e0b3f27dff17eeae0cc2/raw/82c2120a9a554d470ab28199d9ee9d6bb41f9d0f/QLD_meta.hjson"
            self.read_metadata(self.metadata)

    # TODO: Create notebooks for lower level use
    # TODO: Run flags that affect the config object functions
    def run(self,
            aus=True,
            deposits="Fe,Cu,Au,NONE",
            dtb='',
            orientation_decimate=0,
            contact_decimate=5,
            intrusion_mode=0,
            interpolation_spacing=500,
            misorientation=30,
            interpolation_scheme='scipy_rbf',
            fault_decimate=5,
            min_fault_length=5000,
            fault_dip=90,
            pluton_dip=45,
            pluton_form='domes',
            dist_buffer=10,
            contact_dip=-999,
            contact_orientation_decimate=5,
            null_scheme='null',
            thickness_buffer=5000,
            max_thickness_allowed=10000,
            fold_decimate=5,
            fat_step=750,
            close_dip=-999,
            use_interpolations=True,
            use_fat=True
            ):

        if self.quiet:
            print("enabling quiet mode")
            enable_quiet_mode()

        with tqdm(total=100, position=0) as pbar:
            pbar.update(0)

            print("Generating topology analyser input...")
            self.config.export_csv()
            self.config.run_map2model(deposits, aus)
            pbar.update(10)

            self.config.load_dtm()
            pbar.update(10)

            self.config.join_features()
            pbar.update(10)

            self.config.calc_depth_grid(dtb)
            pbar.update(10)

            self.config.export_orientations(orientation_decimate)
            pbar.update(10)
            self.config.export_contacts(contact_decimate, intrusion_mode)
            pbar.update(10)
            self.config.test_interpolation(
                interpolation_spacing, misorientation, interpolation_scheme)
            pbar.update(10)

            self.config.export_faults(
                fault_decimate, min_fault_length, fault_dip)
            self.config.process_plutons(
                pluton_dip, pluton_form, dist_buffer, contact_decimate)
            pbar.update(20)

            # Seismic section is in the hamersely model area
            # TODO: Implement this option better and test with Turner
            if (self.workflow['seismic_section']):
                self.config.extract_section_features(seismic_line_file="",
                                                     seismic_bbox_file="",
                                                     seismic_interp_file="")

            if(self.workflow['contact_dips']):
                self.config.propagate_contact_dips(
                    contact_dip, contact_orientation_decimate)

            if(self.workflow['formation_thickness']):
                self.config.calc_thickness(
                    contact_decimate, null_scheme, thickness_buffer, max_thickness_allowed)

            if(self.workflow['fold_axial_traces']):
                self.config.create_fold_axial_trace_points(
                    fold_decimate, fat_step, close_dip)

            # Prepocess model inputs
            inputs = ('')
            if(self.workflow['model_engine'] == 'geomodeller'):
                inputs = ('invented_orientations', 'intrusive_orientations',
                          'fat_orientations', 'fault_tip_contacts', 'contact_orientations')
            elif(self.workflow['model_engine'] == 'loopstructural'):
                inputs = ('invented_orientations',
                          'fat_orientations', 'contact_orientations')
            elif(self.workflow['model_engine'] == 'gempy'):
                inputs = ('invented_orientations', 'interpolated_orientations',
                          'fat_orientations', 'contact_orientations')
            elif(self.workflow['model_engine'] == 'noddy'):
                inputs = ('')

            self.config.postprocess(
                inputs, self.workflow, use_interpolations, use_fat)
            pbar.update(10)
            self.config.create_map_cmap()

        # self.config.update_projectfile()
        # self.config.export_png()

        disable_quiet_mode()
