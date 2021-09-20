import time
import os
import sys
import re
import logging
import urllib.request
import warnings
import hjson
from tqdm import tqdm

import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
from shapely.geometry import Polygon
from matplotlib import pyplot as plt
from .topology import Topology
from .map2graph import Map2Graph
from . import (
    geology_loopdata,
    structure_loopdata,
    fault_loopdata,
    fold_loopdata,
    mindep_loopdata,
    metafiles,
    clut_paths,
)

from .config import Config
from .m2l_utils import display, enable_quiet_mode, disable_quiet_mode, print
from .m2l_geometry import combine_point_data, update_fault_layer

class Project(object):
    """A high level object implementation of the map2loop workflow."""

    def __init__(
        self,
        loopdata_state=None,
        geology_file=None,
        fault_file=None,
        fold_file=None,
        structure_file=None,
        dtm_file=None,
        mindep_file=None,
        section_files=None,
        drillhole_file=None,
        dtb_grid_file=None,
        cover_map_file=None,
        metadata=None,
    ):
        """Creates project that defines the shared source data.

        :param loopdata_state: Indicates use of loop remote sources and which Australian state to use, defaults to None
        :type loopdata_state: string, optional
        :param dtm_file: Local path or URL to digital terrain model source data, defaults to None
        :type dtm_file: string, optional
        :param geology_file: Local path or URL to stratigraphic source data, defaults to None
        :type geology_file: string, optional
        :param fault_file:  Local path or URL to fault source data[description], defaults to None
        :type fault_file: string, optional
        :param fold_file:  Local path or URL to fold source data[description], defaults to None
        :type fold_file: string, optional
        :param structure_file:  Local path or URL to orientation source data, defaults to None
        :type structure_file: string, optional
        :param mindep_file: Local path or URL to mineral deposit source data, defaults to None
        :type mindep_file: string, optional
        :param section_files: Local path or URL to one or more sets of three geological section source data files, defaults to None
        :type section_files: list of strings, optional
        :param drillhole_file: Local path or URL to drillhole source data, defaults to None
        :type drillhole_file: string, optional
        :param dtb_grid_file: Local path or URL to depth to basement data, defaults to None
        :type dtb_grid_file: string, optional
        :param cover_map_file: Local path or URL to map of limit of cover source data, defaults to None
        :type cover_map_file: string, optional
        :param metadata: File that describes the attributes (column names) in given local or remote sources, defaults to None
        :type metadata: string, optional
        """

        warnings.filterwarnings("ignore")

        if loopdata_state is None:
            self.local = True
            self.state = None
            if any(
                source is None
                for source in [
                    geology_file,
                    fault_file,
                    fold_file,
                    structure_file,
                    metadata,
                ]
            ):
                sys.exit(
                    "Please pass local file paths or urls and a metadata file as input params if you do not wish to use Loop's remote sources."
                )
        else:
            self.local = False
            if loopdata_state in ["WA", "NSW", "VIC", "SA", "QLD", "ACT", "TAS"]:
                self.state = loopdata_state
            else:
                sys.exit(
                    "Valid state code not found, expected 'WA', 'NSW', 'VIC', 'SA', 'QLD', 'ACT' or 'TAS'"
                )

        # If remote, these will be set to null for now, otherwise set to local paths
        self.dtm_file = dtm_file
        self.geology_file = geology_file
        self.fault_file = fault_file
        self.fold_file = fold_file
        self.structure_file = structure_file
        self.mindep_file = mindep_file
        self.section_files = section_files
        self.drillhole_file = drillhole_file
        self.dtb_grid_file = dtb_grid_file
        self.cover_map_file = cover_map_file

        meta_error = "When using your own local files or remote data, pass the path or url to a valid metadata file (.json or .hjson) that describes your input column names.\n"
        meta_error += "You can find an example config here https://gist.github.com/yohanderose/a127c29cb88529f049a5bafc881bb1a0"
        # Load in dictionary that describes the column names
        if (metadata is None) and (geology_file is not None):
            # Check metadata is provided if non loop sources are given
            sys.exit(meta_error)
        elif not self.local:
            # Pass if using loop sources that will be set further down when bounding box is set
            pass
        else:
            # Try to read metadata if given
            self.read_metadata(metadata)

        self.set_proj_defaults()

    def set_proj_defaults(self):
        """Set the bounds and projection of the input data to use if not explicitly provided."""
        # If local, set the maximum bounding box to the bounds of the input files
        # If remote, set bounds to zeros, hopefully this gives the whole map...
        if self.geology_file is not None:
            # print("not remote")
            geology = gpd.read_file(self.geology_file)
            self.proj_bounds = geology.total_bounds
            self.proj_crs = geology.crs
        else:
            # TODO: doesn't like zeros
            self.proj_bounds = (0, 0, 0, 0)
            # TODO: If remote, make epsg required
            self.proj_crs = {"init": "EPSG:28350"}

        self.step_out = 0.1

        # Make matplotlib comply with interface/cmd line window managers
        import matplotlib

        gui_env = ["Qt4Agg", "TkAgg", "GTK3Agg", "WXAgg"]
        all_backends = list(set([*gui_env, *matplotlib.rcsetup.all_backends]))

        for gui in all_backends:
            try:
                matplotlib.use(gui, warn=False, force=True)
                from matplotlib import pyplot as plt

                break
            except:
                continue
        # print("Using:", matplotlib.get_backend())

    def update_workflow(self, model_engine):
        """Set unique run flags depending on model engine to tailor outputs.

        :param workflow: Dict containing desired engine to be updated with flags, defaults to {'model_engine': 'loopstructural'}
        :type workflow: dict

        """

        workflow = {'model_engine': model_engine}

        if (workflow['model_engine'] == 'geomodeller'):
            workflow.update({
                'seismic_section': False,
                'cover_map': False,
                'near_fault_interpolations': True,
                'fold_axial_traces': False,
                'stereonets': False,
                'formation_thickness': True,
                'polarity': False,
                'strat_offset': False,
                'contact_dips': True,
                'drillholes': False
            })
        elif (workflow['model_engine'] == 'loopstructural'):
            workflow.update({
                'seismic_section': False,
                'cover_map': False,
                'near_fault_interpolations': False,
                'fold_axial_traces': False,
                'stereonets': False,
                'formation_thickness': True,
                'polarity': False,
                'strat_offset': True,
                'contact_dips': True,
                'drillholes': False,
                'cover_contacts':True,
                'cover_orientations':True
            })
        elif (workflow['model_engine'] == 'gempy'):
            workflow.update({
                'seismic_section': False,
                'cover_map': False,
                'near_fault_interpolations': False,
                'fold_axial_traces': True,
                'stereonets': False,
                'formation_thickness': False,
                'polarity': False,
                'strat_offset': False,
                'contact_dips': True,
                'drillholes': False
           })
        elif (workflow['model_engine'] == 'noddy'):
            workflow.update({
                'seismic_section': False,
                'cover_map': False,
                'near_fault_interpolations': False,
                'fold_axial_traces': False,
                'stereonets': False,
                'formation_thickness': False,
                'polarity': False,
                'strat_offset': False,
                'contact_dips': True,
                'drillholes': False
            })
        else:
            workflow.update(
                {
                    "seismic_section": False,
                    "cover_map": False,
                    "near_fault_interpolations": False,
                    "fold_axial_traces": False,
                    "stereonets": True,
                    "formation_thickness": True,
                    "polarity": False,
                    "strat_offset": True,
                    'contact_dips': True,
                    'drillholes': False
                }
            )

        self.workflow = workflow

    def update_config(self,
                      out_dir,
                      overwrite='false',
                      loopFilename=None,
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
                      quiet='None',
                      clut_path='',
                      model_engine='loopstructural',
                      run_flags=None,
                      **kwargs):
        """Creates a 'sub-project' Config object and pre-processes input data for some area.

        :param out_dir: Path to write output files to.
        :type out_dir: string
        :param overwrite: Allow overwriting the given out_dir if it exists, false, true or in-place, defaults to false
        :type overwrite: string, optional
        :param bbox_3d: 3D bounding box of coordinates and base/top values defining the area, defaults to { "minx": 0, "maxx": 0, "maxx": 0, "maxy": 0, "base": -10000, "top": 1200, }
        :type bbox_3d: dict, optional
        :param dtm_crs: Set the projection of the dtm, defaults to {'init': 'EPSG:4326'}
        :type dtm_crs: dict, optional
        :param proj_crs: Set the projection of the input data, defaults to None
        :type proj_crs: dict, optional
        :param step_out: How far to consider outside the re-projected dtm, defaults to None
        :type step_out: int, optional
        :param quiet: Allow or block print statements and matplotlib figures, 'None' to quiet nothing, 'all' to quiet everything, 'no-figures' to disable plots and allow text output. Defaults to 'None'
        :type quiet: string, optional
        :param clut_path: Path to custom map colours file
        :type clut_path: string, optional
        :param model_engine: Which modelling engine to use and set associated flags for, defaults to loopstructural
        :type model_engine: string, optional
        :param run_flags: Global dictionary that defines custom params such as decimates and fault length see https://github.com/Loop3D/map2loop-2/issues/56  
        :type run_flags: dict, optional
        :param **kwargs:
        """

        self.update_workflow(model_engine)

        self.loopFilename = loopFilename
        if self.loopFilename is not None:
            if not os.path.exists(loopFilename):
                sys.exit("That project file path does not exist.")

        if bbox_3d["minx"] == 0 and bbox_3d["maxx"] == 0:
            bbox_3d.update(
                {
                    "minx": self.proj_bounds[0],
                    "minx": self.proj_bounds[1],
                    "minx": self.proj_bounds[2],
                    "minx": self.proj_bounds[3],
                }
            )

        self.clut_path = ""

        if proj_crs is None:
            proj_crs = self.proj_crs

        if step_out is None:
            step_out = self.step_out

        self.quiet = quiet

        bbox = tuple(
            [bbox_3d["minx"], bbox_3d["miny"], bbox_3d["maxx"], bbox_3d["maxy"]]
        )
        minx, miny, maxx, maxy = bbox
        lat_point_list = [miny, miny, maxy, maxy, maxy]
        lon_point_list = [minx, maxx, maxx, minx, minx]
        bbox_geom = Polygon(zip(lon_point_list, lat_point_list))
        polygon = gpd.GeoDataFrame(
            index=[0], crs=proj_crs, geometry=[bbox_geom])

        # Define the url queries if remote flag is set
        if self.geology_file is None:
            self.fetch_sources(bbox)

        # TODO: Make run flags global vars that can be updated here instead of in run
        if clut_path != "":
            self.clut_path = clut_path

        try:
            # Check if (perhaps editted) run_flags already exist
            self.run_flags = self.run_flags
        except Exception:
            # Otherwise set them up
            self.run_flags = {
                'aus': True,
                'deposits': "Fe,Cu,Au,NONE",
                'dtb': '',
                'dtb_null':-2147483648,
                'cover_dip':10,
                'cover_spacing':5000,
                'orientation_decimate': 0,
                'contact_decimate': 5,
                'intrusion_mode': 0,
                'interpolation_spacing': 500,
                'misorientation': 30,
                'interpolation_scheme': 'scipy_rbf',
                'fault_decimate': 5,
                'min_fault_length': 5000,
                'fault_dip': 90,
                'fault_orientation_clusters':2,
                'fault_length_clusters':2,
                'pluton_dip': 45,
                'pluton_form': 'domes',
                'dist_buffer': 10,
                'contact_dip': -999,
                'contact_orientation_decimate': 5,
                'null_scheme': 'null',
                'thickness_buffer': 5000,
                'max_thickness_allowed': 10000,
                'fold_decimate': 5,
                'fat_step': 750,
                'close_dip': -999,
                'use_interpolations': True,
                'use_fat': True,
                'use_roi_clip': False,
                'roi_clip_path':'',
                'drift_prefix':['None'],
                'fault_fault_weight':3,
                'fault_weight':1,
                'formation_weight':7,
                'formation_formation_weight':9,
                'fault_formation_weight':5
            }

        # And copy in any new settings from the user
        if run_flags is not None:
            try:
                for key in self.run_flags.keys():
                    try:
                        self.run_flags[key] = run_flags[key]
                    except Exception:
                        pass
            except Exception:
                print('run_flags must be a dictionary, setting defaults.')

        kwargs = {'clut_path': self.clut_path,
                  'run_flags': self.run_flags}
        self.config = Config(out_dir, overwrite, self.geology_file,
                             self.fault_file, self.fold_file,
                             self.structure_file, self.mindep_file, self.section_files,self.drillhole_file, 
                             self.dtb_grid_file,self.cover_map_file,bbox_3d,
                             polygon, step_out, dtm_crs, proj_crs, self.local,
                             self.quiet, self.loopFilename, self.c_l, **kwargs)

        self.config.preprocess()

    def read_metadata(self, filename):
        """Helper function that turns json and hjson files into usable configuration dictionaries.

        :param filename: Path or url to metadata file.
        :type filename: string

        """
        try:
            if filename.startswith("http"):
                with urllib.request.urlopen(filename) as raw_data:
                    self.c_l = hjson.load(raw_data)
                    if self.state is not None:
                        # Check for if remote sources are given as local files, state won't exist
                        if self.state == "SA":
                            for key in self.c_l.keys():
                                try:
                                    self.c_l[key] = self.c_l[key].lower()
                                except Exception as e:
                                    pass

            else:
                with open(filename) as raw_data:
                    self.c_l = hjson.load(raw_data)

            # self.cols = []
            # geol = gpd.read_file(self.geology_file)
            # self.cols += list(geol.columns)
            # struct = gpd.read_file(self.structure_file)
            # self.cols += list(struct.columns)
            # faults = gpd.read_file(self.fault_file)
            # self.cols += list(faults.columns)
            # for key in self.c_l.keys():
            # if self.c_l[key] not in self.cols:
            # try:
            # print(self.c_l[key].lower())
            # except Exception as e:
            # print(self.c_l[key])
            # sys.exit('done')
        except Exception as e:
            sys.exit(e)

    def fetch_sources(self, bbox):
        """Fetch remote loop geospatial data and metadata.

        :param bbox: 2D list of ints or floats describing the area of interest.
        :type bbox: list

        """
        bbox_str = "{},{},{},{}".format(bbox[0], bbox[1], bbox[2], bbox[3])
        self.geology_file = geology_loopdata[self.state].replace(
            "bbox=", "bbox=" + bbox_str
        )
        self.structure_file = structure_loopdata[self.state].replace(
            "bbox=", "bbox=" + bbox_str
        )
        self.fault_file = fault_loopdata[self.state].replace(
            "bbox=", "bbox=" + bbox_str
        )
        self.fold_file = fold_loopdata[self.state].replace(
            "bbox=", "bbox=" + bbox_str)
        self.mindep_file = mindep_loopdata[self.state].replace(
            "bbox=", "bbox=" + bbox_str
        )
        self.metadata = metafiles[self.state]
        self.clut_path = clut_paths[self.state]

        if self.metadata is not "":
            self.read_metadata(self.metadata)

    def run(self):

        if self.quiet == "all":
            enable_quiet_mode()

        with tqdm(total=100, position=0) as pbar:
            pbar.update(0)

            print("Generating topology analyser input...")
            self.config.export_csv()
            self.config.run_map2model(
                self.run_flags['deposits'], self.run_flags['aus'])
            #Map2Graph.map2graph('./test_m2g',self.config.geology_file,self.config.fault_file,self.config.mindep_file,self.c_l,self.run_flags['deposits'])
            pbar.update(10)

            self.config.load_dtm(self.dtm_file if self.local else self.state)
            pbar.update(10)

            self.config.join_features()
            pbar.update(10)

            self.config.calc_depth_grid(self.workflow['cover_map'])
            pbar.update(10)

            self.config.export_orientations(
                self.run_flags['orientation_decimate'],self.workflow['cover_map'])
            pbar.update(10)
            self.config.export_contacts(
                self.run_flags['contact_decimate'], self.run_flags['intrusion_mode'],self.workflow['cover_map'])
            pbar.update(10)
            self.config.test_interpolation(self.run_flags['interpolation_spacing'],
                                           self.run_flags['misorientation'],
                                           self.run_flags['interpolation_scheme'],self.workflow['cover_map'])
            pbar.update(10)

            # TODO: make all these internal, the config class already has the run_flags dictionary
            self.config.export_faults(self.run_flags['fault_decimate'], self.run_flags['min_fault_length'],
                                      self.run_flags['fault_dip'],self.workflow['cover_map'])
            self.config.process_plutons(self.run_flags['pluton_dip'], self.run_flags['pluton_form'], self.run_flags['dist_buffer'],
                                        self.run_flags['contact_decimate'],self.workflow['cover_map'])
            pbar.update(20)

            # Seismic section 
            if (self.workflow['seismic_section']):
                self.config.extract_section_features(self.section_files)

            if self.workflow["contact_dips"]:
                self.config.propagate_contact_dips(
                    self.run_flags['contact_dip'], self.run_flags['contact_orientation_decimate'],self.workflow['cover_map'])

            if (self.workflow['formation_thickness']):
                self.config.calc_thickness(self.run_flags['contact_decimate'], self.run_flags['null_scheme'],
                                           self.run_flags['thickness_buffer'],
                                           self.run_flags['max_thickness_allowed'], self.c_l,self.workflow['cover_map'])

            if self.workflow["fold_axial_traces"]:
                self.config.create_fold_axial_trace_points(
                    self.run_flags['fold_decimate'], self.run_flags['fat_step'], self.run_flags['close_dip'],self.workflow['cover_map'])

            if (self.workflow['drillholes']):
                self.config.extract_drillholes(self,self.drillhole_file)

            # Preprocess model inputs
            inputs = ('')
            if (self.workflow['model_engine'] == 'geomodeller'):
                inputs = ('invented_orientations', 'intrusive_orientations',
                          'fat_orientations', 'fault_tip_contacts',
                          'contact_orientations')
            elif (self.workflow['model_engine'] == 'loopstructural'):
                inputs = ('invented_orientations', 'fat_orientations','intrusive_orientations',
                          'contact_orientations','cover_orientations','cover_contacts')
            elif (self.workflow['model_engine'] == 'gempy'):
                inputs = ('invented_orientations', 'interpolated_orientations',
                          'fat_orientations', 'contact_orientations')
            elif (self.workflow['model_engine'] == 'noddy'):
                inputs = ('')

            self.config.postprocess(inputs, self.workflow, self.run_flags['use_interpolations'],
                                    self.run_flags['use_fat'])
            pbar.update(10)

            self.config.save_cmap(self.workflow['cover_map'])

            # Combine multiple outputs into single graph
            config_out={
                'project_path':os.path.normcase(self.config.project_path),
                'geology_file':os.path.normcase(self.config.geology_file),
                'fault_file':os.path.normcase(self.config.fault_file),
                'fold_file':os.path.normcase(self.config.fold_file),
                'structure_file':os.path.normcase(self.config.structure_file),
                'mindep_file':os.path.normcase(self.config.mindep_file),
                'section_files':self.config.section_files,
                'drillhole_file':self.config.drillhole_file,
                'dtb_grid_file':self.config.dtb_grid_file,
                'cover_map_file':self.config.cover_map_file,
                'bbox_3d':self.config.bbox_3d,
                'step_out':self.config.step_out,
                'dtm_crs':self.config.dtm_crs,
                'proj_crs':self.config.proj_crs,
                'local':self.config.local
            }
            try:
                point_data=combine_point_data(self.config.output_path,self.config.tmp_path)
            except:
                print("combine_point_data failed")          
            
            try:
                Gloop=Topology.make_Loop_graph(self.config.tmp_path,self.config.output_path,self.run_flags['fault_orientation_clusters'],
                                            self.run_flags['fault_length_clusters'],point_data,os.path.normcase(self.config.dtm_reproj_file),self.config.proj_crs,
                                            self.config.c_l,self.config.run_flags,config_out,self.config.bbox_3d)
                nx.write_gml(Gloop, os.path.join(self.config.output_path,'loop.gml'))
                Topology.colour_Loop_graph(self.config.output_path,'loop')
            except:
                print("Topology.make_Loop_graph failed")          
            
            try:
                Map2Graph.map2graph(self.config.output_path,self.config.geology_file,self.config.fault_file,
                    self.config.mindep_file,self.config.c_l,self.config.run_flags['deposits'],
                    self.config.run_flags['fault_orientation_clusters'],self.config.run_flags['fault_length_clusters'],
                    self.config.run_flags['fault_fault_weight'],
                    self.config.run_flags['fault_weight'],
                    self.config.run_flags['formation_weight'],
                    self.config.run_flags['formation_formation_weight'],
                    self.config.run_flags['fault_formation_weight'])
            except:
                print("Topology.map2graph failed")          

            try:
                Map2Graph.granular_map2graph(self.config.output_path,self.config.geology_file,self.config.fault_file,self.config.mindep_file,self.config.c_l,self.config.run_flags['deposits'],
                    self.config.run_flags['fault_fault_weight'],
                    self.config.run_flags['fault_weight'],
                    self.config.run_flags['formation_weight'],
                    self.config.run_flags['formation_formation_weight'],
                    self.config.run_flags['fault_formation_weight'])
            except:
                print("Topology.granular_map2graph failed")          

            try:
                update_fault_layer(self.config.tmp_path,self.config.output_path,self.c_l)
            except:
                print("update_fault_layer failed")          

            try:              
                self.config.update_projectfile()
            except:
                print("self.config.update_projectfile failed")          

            try:              
                self.config.export_png()
            except:
                print("self.config.export_png failed")          

        disable_quiet_mode()
