from logging.handlers import WatchedFileHandler
import shutil
import os
import warnings
from tqdm import tqdm

import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from shapely.geometry import Polygon
from .topology import Topology
import map2model
from . import (m2l_interpolation, m2l_utils, m2l_geometry)
from .map2graph import Map2Graph
from . import (geology_loopdata, structure_loopdata, fault_loopdata, fold_loopdata,
               mindep_loopdata, metafiles, clut_paths, m2l_export)

from .mapdata import MapData
from .config import Config
from .m2l_enums import Datastate, Datatype, ErrorState, VerboseLevel

import beartype

def warning_without_codeline(message, category, filename, lineno, line=''):
    return str(message)
warnings.formatwarning = warning_without_codeline

from map2loop import mapdata

class Project(object):
    """ A high level object implementation of the map2loop workflow.
    """

    @m2l_utils.timer_decorator
    @beartype.beartype
    def __init__(
        self,
        loopdata_state:str="",
        verbose_level:VerboseLevel=VerboseLevel.ALL,
        project_path:str=".\\m2l_data",
        overwrite:str='check',
        model_engine:str='loopstructural',
        working_projection=None,
        geology_filename:str="",
        fault_filename:str="",
        fold_filename:str="",
        structure_filename:str="",
        dtm_filename:str="",
        mindep_filename:str="",
        section_filenames:str="",
        drillhole_filename:str="",
        dtb_grid_filename:str="",
        cover_map_filename:str="",
        metadata_filename:str="",
        loop_project_filename:str="",
        **kwargs):
        """ Creates project that defines the shared source data.

        Parameters
        ----------
        loopdata_state: string, optional
            Indicates use of loop remote sources and which Australian state to use, defaults to blank string
        quiet: string, optional, deprecated
            A filter for turning off matplotlib graphs and log output to console, defaults to "None" (options "None", "no-figures", "All")
        verbose_level: m2l_enums.VerboseLevel, optional
             An enum to set the filter for turning off matplotlib graphs and log output to console, defaults to VerboseLevel.ALL
        project_path: string, optional
            The path to store temporary project data, defaults to ".\\m2l_data"
        overwrite: string, optional
            A string setting whether the project should re-create, abort if present or append to the project directory if it exists, defaults to "check" (options "true", "false", "in-place", "check")
        dtm_filename: string, optional
            Local path or URL to digital terrain model source data, defaults to blank string
        geology_filename: string, optional
            Local path or URL to stratigraphic source data, defaults to blank string
        fault_filename: string, optional
            Local path or URL to fault source data[description], defaults to blank string
        fold_filename: string, optional
            Local path or URL to fold source data[description], defaults to blank string
        structure_filename: string, optional
            Local path or URL to orientation source data, defaults to blank string
        mindep_filename: string, optional
            Local path or URL to mineral deposit source data, defaults to blank string
        section_filenames: list of strings, optional
            Local path or URL to one or more sets of three geological section source data files, defaults to blank string
        drillhole_filename: string, optional
            Local path or URL to drillhole source data, defaults to blank string
        dtb_grid_filename: string, optional
            Local path or URL to depth to basement data, defaults to blank string
        cover_map_filename: string, optional
            Local path or URL to map of limit of cover source data, defaults to blank string
        metadata_filename: string, optional
            Local path or UEL to file that describes the attributes (column names) in given local or remote sources, defaults to blank string
        loop_project_filename: string, optional
            The filename of the loop project file to export to, defaults to blank string

        Returns
        -------
        """
        # Deal with deprecated parameters
        if "quiet" in kwargs:
            warnings.warn("quiet is deprecated and will be removed in later versions of map2loop.  Please use the parameter verbose_level instead", DeprecationWarning)
            if kwargs['quiet'] == "None":
                verbose_level = VerboseLevel.ALL
            elif kwargs['quiet'] == "no-figures":
                verbose_level = VerboseLevel.TEXTONLY
            elif kwargs['quiet'] == "all":
                verbose_level = VerboseLevel.NONE

        self.errorState = ErrorState.NONE
        self.errorStateMsg = ""
        self.map_data = MapData()

        # Sanity check on working projection parameter
        if type(working_projection) == str:
            self.map_data.set_working_projection(working_projection)
        elif type(working_projection) == dict:
            warnings.warn("NOTE: the 'working_projection' parameter will soon no longer accept projection in {'init':'EPSG:28350'} format, please use the string 'EPSG:28350' format instead")
            if 'init' in working_projection:
                self.map_data.set_working_projection(working_projection['init'])
            else:
                warnings.warn("working_projection dict not formatted correctly, please set projection in 'EPSG:28350' format")
        else:
            warnings.warn("No working projection set, will attempt to use the projection of the geology map")

        self.config = Config(self.map_data,verbose_level)
        self.map_data.set_config(self.config)

        self.state = loopdata_state
        if self.state in ["WA", "NSW", "VIC", "SA", "QLD", "ACT", "TAS"]:
            self.use_remote = True
            self.map_data.set_filename(Datatype.GEOLOGY,geology_loopdata[self.state])
            self.map_data.set_filename(Datatype.STRUCTURE,structure_loopdata[self.state])
            self.map_data.set_filename(Datatype.FAULT,fault_loopdata[self.state])
            self.map_data.set_filename(Datatype.FOLD,fold_loopdata[self.state])
            self.map_data.set_filename(Datatype.MINERAL_DEPOSIT,mindep_loopdata[self.state])
            self.map_data.set_filename(Datatype.METADATA,metafiles[self.state])
            self.map_data.set_filename(Datatype.DTM,"AU")

        # Set or replace filenames with incoming filenames
        if geology_filename != "":
            self.map_data.set_filename(Datatype.GEOLOGY,geology_filename)
        if structure_filename != "":
            self.map_data.set_filename(Datatype.STRUCTURE,structure_filename)
        if fault_filename != "":
            self.map_data.set_filename(Datatype.FAULT,fault_filename)
        if fold_filename != "":
            self.map_data.set_filename(Datatype.FOLD,fold_filename)
        if mindep_filename != "":
            self.map_data.set_filename(Datatype.MINERAL_DEPOSIT,mindep_filename)
        if metadata_filename != "":
            self.map_data.set_filename(Datatype.METADATA,metadata_filename)
        if dtm_filename != "":
            self.map_data.set_filename(Datatype.DTM,dtm_filename)
        self.map_data.set_filename(Datatype.SECTION,section_filenames)
        self.map_data.set_filename(Datatype.DRILLHOLE,drillhole_filename)
        self.map_data.set_filename(Datatype.DTB_GRID,dtb_grid_filename)
        self.map_data.set_filename(Datatype.COVER_MAP,cover_map_filename)

        self.topology_graph = None
        self.topology_graph_labels = []
        self.overwrite = overwrite
        self.project_path = project_path
        if type(self.overwrite) == bool:
            print("WARNING: Overwrite should be a string value {\"true\", \"in-place\", \"false\"} ...")
            if self.overwrite == True:
                self.overwrite = "true"
            else:
                self.overwrite = "false"
                self.check_overwrite()

        # Create proj root dir if it doesn't exist
        if (not os.path.exists(project_path)):
            os.mkdir(project_path)
        elif overwrite == "in-place":
            # Pass if proj root exists and complete overwrite not wanted
            pass
        elif overwrite == "true":
            # Remove if exists and accept user's direction
            self.remove_directory_and_wait(project_path)
            os.mkdir(project_path)
        elif overwrite == "check":
            self.check_overwrite()

        # Create each sub directory as needed
        for subdir in ["tmp","data","output","dtm","vtk","graph"]:
            if not os.path.isdir(os.path.join(project_path,subdir)):
                os.mkdir(os.path.join(project_path,subdir))
        self.project_path = project_path

        self.local = True
        if self.state in ["WA", "NSW", "VIC", "SA", "QLD", "ACT", "TAS"]:
            self.local = False
        self.loop_project_filename = loop_project_filename
        self.update_workflow(model_engine)
        self.setup_matplotlib()

        # Check that sufficient files exist to proceed with map2loop process
        # If not set error state so that further processing stops with error message
        if self.map_data.get_filename(Datatype.GEOLOGY) == "" or \
           self.map_data.get_filename(Datatype.STRUCTURE) == "" or \
           self.map_data.get_filename(Datatype.FAULT) == "" or \
           self.map_data.get_filename(Datatype.METADATA) == "":
        #    self.map_data.get_filename(Datatype.FOLD) == "" or \
        #    self.map_data.get_filename(Datatype.MINERAL_DEPOSIT) == "" or \
            self.errorState = ErrorState.URLERROR
            self.errorStateMsg = "Insufficient data files specified on initialisation"
            print("ERROR: Cannot continue processing map as not all datatypes present.  Please check the data specified")
            return

    def check_overwrite(self):
        """ Function to confirm user wants to overwrite previous calculations
        """
        allow = input("Directory \"{}\" exists, overwrite? (y/[n])".format(self.project_path))
        if allow == "y":
            self.remove_directory_and_wait(self.project_path)
            os.mkdir(self.project_path)
            self.overwrite = "true"
        else:
            self.errorState = ErrorState.CONFIGERROR
            self.errorStateMsg = 'Either set overwrite to true or specify a different output_path.'

    @beartype.beartype
    def remove_directory_and_wait(self,directory:str=""):
        """ Function to remove a specified directory and wait until it is completely gone.
        Note potential hang state if os cannot release directory

        Parameters
        ----------
        directory: str, optional
            The path to remove and confirm it is removed. Defaults to "".
        """
        if directory != "":
            shutil.rmtree(directory)
            # TODO: Check for hang state adding sleep and counter, raise if wait too long
            while os.path.exists(directory):
                pass

    def setup_matplotlib(self):
        """ Sets the backend of matplotlib by preference of permissive to restrictive licences
        """
        # Make matplotlib comply with interface/cmd line window managers
        import matplotlib

        # Put Qt4Agg last as it includes GPL code through pyqt
        # and is not included in Loop distributions
        gui_env = ["PS", "TkAgg", "GTK3Agg", "WXAgg", "Qt4Agg"]
        all_backends = list(set([*gui_env, *matplotlib.rcsetup.all_backends]))

        for gui in all_backends:
            try:
                matplotlib.use(gui, warn=False, force=True)
                from matplotlib import pyplot as plt
                break
            except:
                continue

    @beartype.beartype
    def update_workflow(self, model_engine:str):
        """ Set unique run flags depending on model engine to tailor outputs.

        Parameters
        ----------
        model_engine: string
            The program or output type for the project to aim for
        """

        workflow = {'model_engine': model_engine,
                    "seismic_section": False,
                    "cover_map": False,
                    "near_fault_interpolations": False,
                    "fold_axial_traces": False,
                    "stereonets": True,
                    "formation_thickness": True,
                    "polarity": False,
                    "strat_offset": True,
                    'contact_dips': True,
                    'drillholes': False,
                    'cover_contacts':True,
                    'cover_orientations':True
        }

        if (workflow['model_engine'] == 'geomodeller'):
            workflow.update({
                'near_fault_interpolations': True,
                'stereonets': False,
                'strat_offset': False,
            })
        elif (workflow['model_engine'] == 'loopstructural'):
            workflow.update({
                'stereonets': False,
            })
        elif (workflow['model_engine'] == 'gempy'):
            workflow.update({
                'fold_axial_traces': True,
                'stereonets': False,
                'formation_thickness': False,
                'strat_offset': False,
           })
        elif (workflow['model_engine'] == 'noddy'):
            workflow.update({
                'stereonets': False,
                'formation_thickness': False,
                'strat_offset': False,
            })
        self.workflow = workflow

    @m2l_utils.timer_decorator
    @beartype.beartype
    def update_config(self,
                      out_dir:str,
                      overwrite=None,
                      loop_project_filename:str="",
                      bbox_3d=None,
                      dtm_crs:str='EPSG:4326',
                      project_crs:str=None,
                      step_out:float=0.1,
                      clut_path='',
                      model_engine="",
                      run_flags={},
                      **kwargs):
        """ Updates the project variables in the config structure needed for the map2loop process to run

        Parameters
        ----------
        out_dir: string
            The name of the sub-directory to store working and output files
        overwrite: string/None, deprecated
            Old flag to indicate whether to overwrite the output directory
        loop_project_filename: string, optional
            The filename of the loop project file to export to
        bbox_3d: tuple/None, optional
            The bounding box of the region to apply map2loop to. Defaults to limits of geology_filename data if not specified
        dtm_crs: string, optional
            The projection of the digital terrain map. Defaults to 'EPSG:4326' (WGS84 Lat/Long)
        project_crs: string
            The projection use internally for map2loop. Defaults to 'EPSG:28350' (MGA Zone 50 ~Western Australia)
        step_out: float, optional
            A buffer area in degrees to expand the digital terrain map to ensure coverage of region of interest. Default is 0.1
        quiet: string, optional, deprecated
            A string to indicate the level of display and printed output ('None' to quiet nothing, 'all' to quiet everything, 'no-figures' to disable plots and allow text output. Defaults to 'None')
        clut_path: string, optional
            The path to a custom map colouring lookup file (crs formatted with columns for formation and colour string)
        model_engine: string, optional
            String to re-define the output format and map2loop process. Overrides the Project init function setting
        run_flags: dict, optional
            Additional dictionary options to use such as decimates and fault length limits see https://github.com/Loop3D/map2loop-2/issues/56
        """
        # Deal with deprecated parameters
        if "quiet" in kwargs:
            warnings.warn("quiet is deprecated and will be removed in later versions of map2loop.  Please use the parameter verbose_level instead", DeprecationWarning)
            if kwargs['quiet'] == "None":
                self.config.verbose_level = VerboseLevel.ALL
            elif kwargs['quiet'] == "no-figures":
                self.config.verbose_level = VerboseLevel.TEXTONLY
            elif kwargs['quiet'] == "all":
                self.config.verbose_level = VerboseLevel.NONE

        if overwrite != None:
            warnings.warn('NOTE: the "overwrite" parameter has moved from the update_config function to the project initialisation call')

        if self.errorState != ErrorState.NONE:
            print("ERROR: Error state -", self.errorStateMsg)
            return

        if model_engine != "":
            self.update_workflow(model_engine)

        if loop_project_filename != "":
            self.loop_project_filename = loop_project_filename

        bbox_valid = False
        if isinstance(bbox_3d,dict) and 'minx' in bbox_3d and 'maxx' in bbox_3d \
                                    and 'miny' in bbox_3d and 'maxy' in bbox_3d:
            if bbox_3d['minx'] < bbox_3d['maxx'] and \
               bbox_3d['miny'] < bbox_3d['maxy'] and \
               bbox_3d['base'] < bbox_3d['top']:
                bbox_valid = True
        if bbox_valid == False:
            print("Invalid bounding box specified, attempting to get one from the map files")
            bbox_3d, project_crs = self.map_data.calculate_bounding_box_and_projection()

        if bbox_3d is None:
            bbox_3d = {'minx':0,'maxx':1000,',miny':0,'maxy':1000,'base':-10000,'top':1200}
        if project_crs is None:
            project_crs = 'EPSG:28350'
            if "proj_crs" in kwargs:
                warnings.warn("proj_crs is deprecated and will be removed in later versions of map2loop.  Please use the parameter project_crs instead", DeprecationWarning)
                project_crs = kwargs['proj_crs']

        bbox_str = "{},{},{},{},{}".format(bbox_3d['minx'], bbox_3d['miny'], bbox_3d['maxx'], bbox_3d['maxy'],project_crs)
        self.map_data.update_filenames_with_bounding_box(bbox_str)

        bbox = tuple(
            [bbox_3d["minx"], bbox_3d["miny"], bbox_3d["maxx"], bbox_3d["maxy"]]
        )
        minx, miny, maxx, maxy = bbox
        lat_point_list = [miny, miny, maxy, maxy, maxy]
        lon_point_list = [minx, maxx, maxx, minx, minx]
        bbox_geom = Polygon(zip(lon_point_list, lat_point_list))
        polygon = gpd.GeoDataFrame(
            index=[0], crs=project_crs, geometry=[bbox_geom])

        if clut_path == "":
            if self.state in ["WA", "NSW", "VIC", "SA", "QLD", "ACT", "TAS"]:
                clut_path = clut_paths[self.state]

        self.config.update(
                            out_dir,
                            bbox_3d,
                            polygon,
                            step_out,
                            dtm_crs,
                            project_crs,
                            self.local,
                            self.map_data.get_filename(Datatype.METADATA),
                            clut_path=clut_path,
                            run_flags=run_flags
                        )

        self.map_data.load_all_map_data(self.config)

    @m2l_utils.timer_decorator
    def run(self):
        """ The main map2loop process for reading, converting, reprojecting and interpolating map data into event, observation
        and relationship data suitable for further processing
        """

        if self.errorState != ErrorState.NONE:
            print("ERROR: Error state set at ", self.errorStateMsg)
            return

        # set orientation/structure data to recreate with dirty flag and unloaded state
        # because half way through run process (merge_structure_with geology)
        # the geopandas is fundamentally changed which breaks export_wkt_format_files call
        # TODO: Fix so that both version of orientations are possible without reloading the file
        self.map_data.data_states[Datatype.STRUCTURE] = Datastate.UNLOADED
        self.map_data.dirtyflags[Datatype.STRUCTURE] = True
        with tqdm(total=100, position=0) as pbar:
            pbar.update(0)

            if self.config.verbose_level != VerboseLevel.NONE:
                print("Generating topology analyser input...")
            self.map_data.export_wkt_format_files()
            self.__run_map2model()
            #Map2Graph.map2graph('./test_m2g',self.config.geology_file,self.config.fault_file,self.config.mindep_file,self.config.c_l,self.config.run_flags['deposits'])
            pbar.update(10) #10%

            # self.map_data.load_dtm(self.config.dtm_filename if self.local else self.state)
            pbar.update(10) #20%

            # self.config.join_features(self.map_data.get_map_data(Datatype.GEOLOGY))
            Topology.save_group(self.topology_graph,
                                self.topology_graph_labels,
                                self.config,
                                self.map_data)
            pbar.update(5) #25%

            self.map_data.calc_depth_grid(self.workflow)
            pbar.update(5) #30%

            self.map_data.merge_structure_with_geology(self.config.c_l)
            self.map_data.export_orientations(self.workflow)
            pbar.update(10) #40%

            self.map_data.export_contacts(self.workflow)
            pbar.update(10) #50%

            self.__test_interpolation()
            pbar.update(10) #60%

            self.map_data.export_faults(self.workflow, self.map_data.dip_grid, self.map_data.dip_dir_grid)
            pbar.update(10) #70%

            self.__process_plutons()
            pbar.update(10) #80%

            # Seismic section
            if (self.workflow['seismic_section']):
                self.__extract_section_features()

            if self.workflow["contact_dips"]:
                self.__propagate_contact_dips()

            if (self.workflow['formation_thickness']):
                self.__calc_thickness()

            if self.workflow["fold_axial_traces"]:
                self.__create_fold_axial_trace_points()

            if (self.workflow['drillholes']):
                self.__extract_drillholes()

            self.__postprocess()

            self.config.save_cmap(self.workflow)

            point_data=m2l_geometry.combine_point_data(self.config.output_path,self.config.tmp_path)

            Gloop = Topology.make_Loop_graph(self.config, self.map_data, point_data)
            nx.write_gml(Gloop, os.path.join(self.config.output_path,'loop.gml'))
            Topology.colour_Loop_graph(self.config.output_path,'loop')

            # TODO: Fix map2graph and make it a proper object
            if(self.config.run_flags['map2graph']):
                try:
                    Map2Graph.map2graph(
                        self.config.output_path,
                        self.config.geology_filename,
                        self.config.fault_filename,
                        self.config.mindep_filename,
                        self.config.c_l,
                        self.config.run_flags['deposits'],
                        self.config.run_flags['fault_orientation_clusters'],
                        self.config.run_flags['fault_length_clusters'],
                        self.config.run_flags['fault_fault_weight'],
                        self.config.run_flags['fault_weight'],
                        self.config.run_flags['formation_weight'],
                        self.config.run_flags['formation_formation_weight'],
                        self.config.run_flags['fault_formation_weight'])
                except:
                    print("Topology.map2graph failed")

            # TODO: Fix map2graph and make it a proper object
            if(self.config.run_flags['granular_map2graph']):
                try:
                    Map2Graph.granular_map2graph(
                        self.config.output_path,
                        self.config.geology_file,
                        self.config.fault_file,
                        self.config.mindep_file,
                        self.config.c_l,
                        self.config.run_flags['deposits'],
                        self.config.run_flags['fault_fault_weight'],
                        self.config.run_flags['fault_weight'],
                        self.config.run_flags['formation_weight'],
                        self.config.run_flags['formation_formation_weight'],
                        self.config.run_flags['fault_formation_weight'])
                except:
                    print("Topology.granular_map2graph failed")

            m2l_geometry.update_fault_layer(self.config, self.map_data)

            self.update_loop_project_file()

            self.__export_png()

            pbar.update(20) #100%

    def update_loop_project_file(self):
        """ A function to convert multiple csv and map2loop output files into a single loop project file
        """
        m2l_export.export_to_projectfile(self.loop_project_filename, self.config)

    def __run_map2model(self):
        mindep_filename = ""
        if self.map_data.get_map_data(Datatype.MINERAL_DEPOSIT) is not None:
            mindep_filename = self.config.mindep_filename_wkt

        # TODO: Simplify. Note: this is external so have to match fix to map2model module
        run_log = map2model.run(self.config.graph_path,
                                self.config.geology_filename_wkt,
                                self.config.fault_filename_wkt,
                                mindep_filename,
                                self.config.bbox_3d,
                                self.config.c_l,
                                self.config.verbose_level != VerboseLevel.NONE,
                                self.config.run_flags['deposits']
                                )
        if self.config.verbose_level != VerboseLevel.NONE:
            print(run_log)

        if self.config.run_flags['aus']:
            if self.config.verbose_level != VerboseLevel.NONE:
                print("Resolving ambiguities using ASUD...", end='\toutput_dir:')
            Topology.use_asud(self.config.strat_graph_filename, self.config.graph_path)
            self.config.strat_graph_filename = os.path.join(self.config.graph_path,'ASUD_strat.gml')
            if self.config.verbose_level != VerboseLevel.NONE:
                print("Done.")

        if self.config.verbose_level != VerboseLevel.NONE:
            print("Generating topology graph display and unit groups...")

        self.topology_graph = nx.read_gml(self.config.strat_graph_filename, label='id')
        selected_nodes = [n for n, v in self.topology_graph.nodes(data=True) if n >= 0]

        if self.config.verbose_level != VerboseLevel.NONE:
            nx.draw_networkx(self.topology_graph,
                             pos=nx.kamada_kawai_layout(self.topology_graph),
                             arrows=True,
                             nodelist=selected_nodes)
            plt.savefig(os.path.join(self.config.tmp_path,"topology-fig.png"))
            print("Topology figure saved to", os.path.join(self.config.tmp_path, "topology-fig.png"))

        if self.config.verbose_level != VerboseLevel.NONE:
            nlist = list(self.topology_graph.nodes.data('LabelGraphics'))
            nlist.sort()
            for node in nlist:
                if node[0] >= 0:
                    elem = str(node[1]).replace("{'text':", "").replace(
                        ", 'fontSize': 14}", "")
                    # second = elem.split(":").replace("'", "")
                    print(node[0], " ", elem)

        # Save groups of stratigraphic units
        groups, self.topology_graph_labels, topology_graph = Topology.get_series(self.config.strat_graph_filename, 'id')
        Topology.save_units(topology_graph, self.topology_graph_labels, self.config)

    # Interpolates a regular grid of orientations from an shapefile of
    # arbitrarily-located points and saves out four csv files of l, m & n
    # direction cosines and dip dip direction data
    def __test_interpolation(self):
        # TODO: Add basal contacts to storage structure
        basal_contacts_filename = os.path.join(self.config.tmp_path, 'basal_contacts.shp')

        group_girdle = m2l_utils.plot_bedding_stereonets(self.config, self.map_data)
        super_groups, self.use_gcode3 = Topology.super_groups_and_groups(group_girdle, self.config, self.map_data, self.workflow)

        orientation_interp, contact_interp, combo_interp = m2l_interpolation.interpolation_grids(
            self.config, self.map_data, basal_contacts_filename, super_groups)

        spacing = self.config.run_flags['interpolation_spacing']  # grid spacing in meters
        x = self.config.width
        y = self.config.height
        bbox = self.config.bbox
        dip_grid = np.ones((y, x))
        dip_grid = dip_grid * -999
        dip_dir_grid = np.ones((y, x))
        dip_dir_grid = dip_dir_grid * -999
        contact_grid = np.ones((y, x))
        contact_grid = contact_grid * -999
        polarity_grid = np.ones((y, x))
        polarity_grid = polarity_grid * -999
        for i, row in combo_interp.iterrows():
            r = int((row[1] - bbox[1]) / spacing)
            c = int((row[0] - bbox[0]) / spacing)
            dip_grid[r, c] = float(row[5])
            dip_dir_grid[r, c] = float(row[6])
            polarity_grid[r,c]=float(row[4])
        for i, row in contact_interp.iterrows():
            r = int((row[1] - bbox[1]) / spacing)
            c = int((row[0] - bbox[0]) / spacing)
            contact_grid[r, c] = float(row[4])

        self.map_data.dip_grid = dip_grid
        self.map_data.dip_dir_grid = dip_dir_grid
        self.map_data.polarity_grid = polarity_grid

        if self.config.verbose_level == VerboseLevel.ALL:
            plt.imshow(dip_grid,
                       cmap="hsv",
                       origin='lower',
                       vmin=-90,
                       vmax=90)
            plt.title('Interpolated Dips')
            plt.show()

            plt.imshow(dip_dir_grid,
                       cmap="hsv",
                       origin='lower',
                       vmin=0,
                       vmax=360)
            plt.title('Interpolated Dip Directions')
            plt.show()

            plt.imshow(contact_grid,
                       cmap="hsv",
                       origin='lower',
                       vmin=-360,
                       vmax=360)
            plt.title('Interpolated Contacts')
            plt.show()

    def __process_plutons(self):
        m2l_geometry.process_plutons(self.config, self.map_data, self.workflow)

    def __extract_section_features(self):
        # Extract faults and basal contacts of groups from seismic sections
        # input geology file (if local)

        for section in self.map_data.get_map_data(Datatype.SECTION):
            fig, ax = plt.subplot(1,1)
            seismic_line_filename = section[0]
            seismic_line = gpd.read_file(seismic_line_filename)  # import map
            seismic_line.plot(figsize=(10, 10),
                              edgecolor='#000000',
                              title="Seismic Line",
                              ax=ax,
                              linewidth=0.2)  # display map
            m2l_utils.display(seismic_line)

            # input geology file (if local)
            seismic_bbox_filename = section[1]
            seismic_bbox = gpd.read_file(seismic_bbox_filename)  # import map
            seismic_bbox.set_index('POSITION', inplace=True)

            # input geology file (if local)
            seismic_interp_filename = section[2]
            seismic_interp = gpd.read_file(seismic_interp_filename)  # import map
            if self.config.verbose_level == VerboseLevel.ALL:
                seismic_interp.plot(column='FEATURE',
                                    figsize=(10, 10),
                                    edgecolor='#000000',
                                    title="Seismic Interpolation",
                                    linewidth=0.5)  # display map
                m2l_utils.display(seismic_interp)

            m2l_geometry.extract_section( self.config,
                                         self.map_data,
                                         seismic_line,
                                         seismic_bbox,
                                         seismic_interp,
                                         surface_cut=2000)

            contacts = pd.read_csv(os.path.join(self.config.output_path, 'contacts4.csv'),sep=",")
            seismic_contacts = pd.read_csv(os.path.join(self.config.output_path, 'seismic_base.csv'), sep=",")
            all_contacts = pd.concat([contacts, seismic_contacts], sort=False)
            all_contacts.to_csv(os.path.join(self.config.output_path, 'contacts4.csv'),
                                index=None,
                                header=True)

            faults = pd.read_csv(os.path.join(self.config.output_path, 'faults.csv'),sep=",")
            seismic_faults = pd.read_csv(os.path.join(self.config.output_path, 'seismic_faults.csv'), sep=",")
            all_faults = pd.concat([faults, seismic_faults], sort=False)
            all_faults.to_csv(os.path.join(self.config.output_path, 'faults.csv'),
                              index=None,
                              header=True)

    def __propagate_contact_dips(self):
        if self.config.verbose_level != VerboseLevel.NONE:
            print("Propagating dips along contacts...")
        m2l_geometry.save_basal_contacts_orientations_csv(self.config, self.map_data, self.workflow)

    def __calc_thickness(self):
        # Estimate formation thickness and normalised formation thickness
        # basal_contacts_filename = os.path.join(self.config.tmp_path, 'basal_contacts.shp')
        m2l_interpolation.save_contact_vectors(self.config, self.map_data, self.workflow)

        # TODO: multi thread / numba jit
        m2l_geometry.calc_thickness_with_grid(self.config, self.map_data)
        m2l_geometry.calc_min_thickness_with_grid(self.config, self.map_data)

        m2l_geometry.normalise_thickness(self.config.output_path)

        if self.config.verbose_level == VerboseLevel.ALL:
            m2l_utils.plot_points(
                os.path.join(self.config.output_path, 'formation_thicknesses_norm.csv'),
                self.map_data.get_map_data(Datatype.GEOLOGY),
                'norm_th', 'x', 'y', True, 'numeric', "Normalised Thicknesses"
            )

    def __create_fold_axial_trace_points(self):
        folds_clip = self.map_data.get_map_data(Datatype.FOLD).copy()
        # folds_clip = gpd.read_file(self.fold_filename)
        if (len(folds_clip) > 0):
            m2l_geometry.save_fold_axial_traces(self.config, self.map_data, self.workflow)

            # TODO : better approximation / multithread / numba
            m2l_geometry.save_fold_axial_traces_orientations(self.config, self.map_data, self.workflow)

    def __postprocess(self):

        if (self.workflow['model_engine'] == 'geomodeller'):
            inputs = ('invented_orientations',
                      'intrusive_orientations',
                      'fat_orientations',
                      'fault_tip_contacts',
                      'contact_orientations')
        elif (self.workflow['model_engine'] == 'loopstructural'):
            inputs = ('invented_orientations',
                      'fat_orientations',
                      'intrusive_orientations',
                      'contact_orientations',
                      'cover_orientations',
                      'cover_contacts')
        elif (self.workflow['model_engine'] == 'gempy'):
            inputs = ('invented_orientations',
                      'interpolated_orientations',
                      'fat_orientations',
                      'contact_orientations')
        elif (self.workflow['model_engine'] == 'noddy'):
            inputs = ('')
        else:
            inputs = ('')

        m2l_geometry.tidy_data(self.config, self.map_data, self.use_gcode3, inputs)
        if(self.workflow['cover_map']):
            self.map_data.dtb.close()

        # Calculate polarity of original bedding orientation data
        if (self.workflow['polarity']):
            if self.config.verbose_level != VerboseLevel.NONE:
                print("Processing polarity")
            m2l_geometry.save_orientations_with_polarity(self.config, self.map_data)

            if self.config.verbose_level == VerboseLevel.ALL:
                m2l_utils.plot_points(
                    os.path.join(self.config.output_path,'orientations_polarity.csv'),
                    self.map_data.get_map_data(Datatype.GEOLOGY),
                    'polarity', 'X', 'Y', True, 'alpha',"Orientation Polarities")

        # Calculate minimum fault offset from stratigraphy and stratigraphic fault offset
        if (self.workflow['strat_offset']):
            if self.config.verbose_level != VerboseLevel.NONE:
                print("Processing strat offsets")
            fault_test = pd.read_csv(os.path.join(self.config.output_path, 'fault_dimensions.csv'), sep=',')
            if (len(fault_test) > 0):

                m2l_geometry.fault_strat_offset(self.config, self.map_data)

                if self.config.verbose_level == VerboseLevel.ALL:
                    m2l_utils.plot_points(
                        os.path.join(self.config.output_path,'fault_strat_offset3.csv'),
                        self.map_data.get_map_data(Datatype.GEOLOGY), 'min_offset', 'X', 'Y', True,
                        'numeric',"Fault Min Offsets")
                    m2l_utils.plot_points(
                        os.path.join(self.config.output_path,'fault_strat_offset3.csv'),
                        self.map_data.get_map_data(Datatype.GEOLOGY), 'strat_offset', 'X', 'Y', True,
                        'numeric',"Fault Strat Offsets")

        # Analyse fault topologies
        if self.config.verbose_level != VerboseLevel.NONE:
            print("Parsing fault relationships")
        Topology.parse_fault_relationships(self.config)
        m2l_geometry.save_interpolation_parameters(self.config)

    def __export_png(self):
        filename = self.loop_project_filename
        if filename == "":
            filename = os.path.join(self.project_path, '{}'.format(os.path.basename(self.project_path)))
        if self.config.verbose_level != VerboseLevel.NONE:
            print("Exporting graphical map...")
        try:
            fig, ax = plt.subplots()
            plt.tight_layout()
            ax.ticklabel_format(axis='both', useOffset=False, style='plain')
            ax.margins(0.0)
            fig.set_facecolor("#ffffff00")
            geology_figure = self.map_data.get_map_data(Datatype.GEOLOGY).copy().plot(
                column=self.config.c_l['c'],
                ax=ax,
                figsize=(10, 10),
                edgecolor='#000000',
                linewidth=0.2,
                cmap=self.config.cmap).get_figure()
            geology_figure.savefig("{}.png".format(filename))
            if self.config.verbose_level != VerboseLevel.NONE:
                print("Geology graphic exported to: {}.png".format(filename))
            if self.config.verbose_level == VerboseLevel.ALL:
                plt.title('Geology Map')
                plt.show()
            else:
                plt.close(geology_figure)
        except Exception as e:
            warnings.warn(e)
            warnings.warn("WARNING: Could not save geology graphic")

    def __extract_drillholes(self):
        dhdb_filename = self.map_data.get_map_data(Datatype.DRILLHOLE) # input drill hole file (if local)
        dh = gpd.read_file(dhdb_filename, bbox=self.bbox)
        if(len(dh) > 0):
            dh = dh[['X', 'Y', 'Z', 'formation']]
            contacts = pd.read_csv(os.path.join(self.config.tmp_path, 'contacts4.csv'))
            all_contacts = pd.concat([contacts, dh], sort=False)
            all_contacts.reset_index(inplace=True)
            all_contacts.drop(labels='index', axis=1, inplace=True)
            all_contacts.to_csv(os.path.join(self.config.tmp_path, 'contacts4.csv'), index=None, header=True)
            if self.config.verbose_level != VerboseLevel.NONE:
                print(len(dh), 'drillhole contacts added')
        else:
            if self.config.verbose_level != VerboseLevel.NONE:
                print('No drillhole data for merging.')
