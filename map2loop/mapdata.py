import geopandas
import rasterio
import fiona
from .config import Config
import beartype
import os
from . import m2l_map_checker
from .m2l_enums import Datatype, Datastate, VerboseLevel
import warnings
from .topology import Topology
from . import (m2l_utils, m2l_geometry, m2l_interpolation)
import time
import matplotlib.pyplot as plt
import numpy
import sys

class MapData:
    """
    A data structure containing all the map data loaded from map files

    Attributes
    ----------
    data: list of geopandas.DataFrames
        A list containing the geopanda data frames or rasterio image of all the different map data
    filenames: list of data source filenames
        A list of the filenames/urls of where to find the data to be loaded
    dirtyflags: list of booleans
        A list of flags indicating whether the data has been loaded or dirtied
    data_states: intEnum
        Enums representing the state of each of the data sets
    working_projection: str
        A string containing the projection e.g. "EPSG:28350"
    config: Config
        A link to the config structure which is defined in config.py
    """
    def __init__(self):
        self.data = [None] * len(Datatype)
        self.filenames = [None] * len(Datatype)
        self.output_filenames = [None] * len(Datatype)
        self.dirtyflags = [True] * len(Datatype)
        self.data_states = [Datastate.UNNAMED] * len(Datatype)
        self.working_projection = None
        self.config = None

        # TODO: Add proper getter/setters to these dataframes
        self.merged_structure = None
        self.basal_contacts = None
        self.basal_contacts_no_faults = None
        self.dtb = None
        self.dtb_null = 0
        self.dip_grid = None
        self.dip_dir_grid = None
        self.polarity_grid = None

    def set_working_projection(self, projection):
        """
        Setter function for the working projection for the map data

        Args:
            projection (int or string): The projection to use for map reprojection
        """
        if type(projection) == int:
            projection = "EPSG:" + str(projection)
        if type(projection) == str:
            self.working_projection = projection
        else:
            warnings.warn(f"Warning: Unknown projection set {projection}. Leaving all map data in original projection\n")

    @beartype.beartype
    def set_config(self, config:Config):
        """
        Setter function for config structure

        Args:
            config (Config) : The configuration structure for the map data
        """
        self.config = config

    @beartype.beartype
    def set_filename(self,datatype:Datatype,filename:str):
        """
        Setter function for the filename for a specific datatype 

        Args:
            datatype (Datatype): The datatype for the filename specified
            filename (str): The filename to store
        """
        if self.filenames[datatype] != filename:
            self.filenames[datatype] = filename
            self.data_states[datatype] = Datastate.UNLOADED
            self.dirtyflags[datatype] = True

    @beartype.beartype
    def get_filename(self,datatype:Datatype):
        """
        Getter function for a filename of a specified datatype

        Args:
            datatype (Datatype): The datatype of the filename wanted

        Returns:
            string: The filename of the datatype specified
        """
        if self.data_states != Datastate.UNNAMED:
            return self.filenames[datatype]
        else:
            warnings.warn(f"Requested filename for {str(type(datatype))} is not set\n")
            return None
    
    @beartype.beartype
    def load_all_map_data(self, config:Config=None):
        """
        Function to load all the map data for each datatype.  Cycles through each type and loads it

        Args:
            config (Config, optional): The config structure to use for map data if not already specified. Defaults to None.
        """
        if config == None:
            config = self.config
        else:
            self.config = config
        # print("Loading data frames")
        for i in [Datatype.GEOLOGY,Datatype.STRUCTURE,Datatype.FAULT,Datatype.FOLD,Datatype.MINERAL_DEPOSIT]:
            self.load_map_data(i)
        for i in [Datatype.DTB_GRID]:
            self.load_rasterio_map_data(i)
        self.load_dtm()

    @beartype.beartype
    def load_map_data(self, datatype:Datatype):
        """
        Function to load map data from file, reproject and clip it and then check data is valid

        Args:
            datatype (Datatype): The datatype to load
        """
        if self.filenames[datatype] == None or self.data_states[datatype] == Datastate.UNNAMED:
            warnings.warn(f"Datatype {datatype.name} is not set and so cannot be loaded\n")
        elif self.dirtyflags[datatype] == True:
            if self.data_states[datatype] == Datastate.UNLOADED:
                # Load data from file
                try:
                    self.data[datatype] = geopandas.read_file(self.filenames[datatype])
                    self.data_states[datatype] = Datastate.LOADED
                except:
                    sys.stdout.flush()
                    sys.stderr.flush()
                    warnings.warn(f"Failed to open {datatype.name} file called '{self.filenames[datatype]}'\n")
                    sys.stderr.flush()
                    warnings.warn(f"Cannot continue as {datatype.name} was not loaded\n")
                    sys.stderr.flush()
                    return
            if self.data_states[datatype] == Datastate.LOADED:
                # Reproject geopanda to required CRS
                self.set_working_projection_on_map_data(datatype)
                self.data_states[datatype] = Datastate.REPROJECTED
            if self.data_states[datatype] == Datastate.REPROJECTED:
                # Clip geopanda to bounding polygon
                self.data[datatype] = geopandas.clip(self.data[datatype],self.config.polygon)
                self.data_states[datatype] = Datastate.CLIPPED
            if self.data_states[datatype] == Datastate.CLIPPED:
                # Convert column names using codes_and_labels dictionary
                self.check_map(datatype)
                self.data_states[datatype] = Datastate.CONVERTED
            if self.data_states[datatype] == Datastate.CONVERTED:
                self.data_states[datatype] = Datastate.COMPLETE
            self.dirtyflags[datatype] = False
    
    @beartype.beartype
    def load_rasterio_map_data(self, datatype:Datatype):
        """
        Function to load rasterio map data from file, reproject and clip it and then check data is valid
        Note: currently no reprojection, clipping or data checking is available (just loading)

        Args:
            datatype (Datatype): The rasterio datatype to load
        """
        if self.filenames[datatype] == None or self.data_states[datatype] == Datastate.UNNAMED:
            warnings.warn(f"Datatype {datatype.name} is not set and so cannot be loaded\n")
        elif self.dirtyflags[datatype] == True:
            if self.data_states[datatype] == Datastate.UNLOADED:
                # Load data from file
                try:
                    self.data[datatype] = rasterio.open(self.filenames[datatype])
                    self.data_states[datatype] = Datastate.LOADED
                except:
                    sys.stdout.flush()
                    sys.stderr.flush()
                    warnings.warn(f"Failed to open {datatype.name} file called '{self.filenames[datatype]}'\n")
                    sys.stderr.flush()
            if self.data_states[datatype] == Datastate.LOADED:
                # Reproject geopanda to required CRS
                # self.set_working_projection_on_map_data(datatype)
                self.data_states[datatype] = Datastate.REPROJECTED
            if self.data_states[datatype] == Datastate.REPROJECTED:
                # Clip geopanda to bounding polygon
                # self.data[datatype] = geopandas.clip(self.data[datatype],self.config.polygon)
                self.data_states[datatype] = Datastate.CLIPPED
            if self.data_states[datatype] == Datastate.CLIPPED:
                # Convert column names using codes_and_labels dictionary
                # self.check_map(datatype)
                self.data_states[datatype] = Datastate.COMPLETE
            self.dirtyflags[datatype] = False

    @beartype.beartype
    def check_map(self,datatype:Datatype):
        """
        Function to check the validity of a map data from file

        Args:
            datatype (Datatype): The rasterio datatype to check
        """
        _warnings = []
        _errors = []
        if datatype == Datatype.GEOLOGY:
            self.data[Datatype.GEOLOGY] = m2l_map_checker.check_geology_map(self.data[Datatype.GEOLOGY],self.config.c_l,self.config.run_flags['drift_prefix'],_warnings,_errors,self.config.verbose_level)
        if datatype == Datatype.STRUCTURE:
            self.data[Datatype.STRUCTURE] = m2l_map_checker.check_structure_map(self.data[Datatype.STRUCTURE],self.config.c_l,_warnings,_errors,self.config.verbose_level)
        if datatype == Datatype.FAULT:
            self.data[Datatype.FAULT] = m2l_map_checker.check_fault_map(self.data[Datatype.FAULT],self.config.c_l,_warnings,_errors,self.config.verbose_level)
        if datatype == Datatype.FOLD:
            self.data[Datatype.FOLD] = m2l_map_checker.check_fold_map(self.data[Datatype.FOLD],self.config.c_l,_warnings,_errors,self.config.verbose_level)
        if datatype == Datatype.MINERAL_DEPOSIT:
            self.data[Datatype.MINERAL_DEPOSIT] = m2l_map_checker.check_mindep_map(self.data[Datatype.MINERAL_DEPOSIT],self.config.c_l,_warnings,_errors,self.config.verbose_level)
        if(len(_warnings) > 0):
            print("\nWarnings:")
            for w in _warnings:
                print("   -> " + str(w))
            sys.stdout.flush()
            warnings.warn(f'The warnings listed above were issued checking the {datatype.name} map\n')
            sys.stderr.flush()
        if(len(_errors) > 0):
            print("\nErrors:")
            for e in _errors:
                print("   -> " +  str(e))
            sys.stdout.flush()
            warnings.warn(f'The errors listed above on {datatype.name} must be fixed prior to rerunning map2loop\n')
            sys.stderr.flush()
            sep = '\12'
            raise NameError(sep.join(_errors) + '\n map2loop error: Fix errors before running again')
        else:
            self.data_states[datatype] = Datastate.CONVERTED

    @beartype.beartype
    def set_working_projection_on_map_data(self,datatype:Datatype):
        if self.working_projection == None:
            print("No working projection set leaving map data in original projection")
        else:
            if self.data_states[datatype] >= Datastate.LOADED:
                if self.data[datatype].crs == None:
                    print("No projection on original map data, assigning working_projection",self.working_projection)
                    self.data[datatype].crs = self.working_projection
                else:
                    self.data[datatype].to_crs(crs=self.working_projection,inplace=True)

    @beartype.beartype
    def save_all_map_data(self,project_dir:str):
        for i in Datatype:
            self.save_map_data(project_dir,i) 

    @beartype.beartype
    def save_map_data(self, project_dir:str, datatype:Datatype):
        if self.data_states[datatype] == Datastate.CONVERTED:
            try:
                filename = os.path.join(project_dir,"tmp",datatype.name,".csv")
                self.data[datatype].write_csv(filename)
            except:
                warnings.warn(f"Failed to save {datatype.name} to file named {filename}\n")

    @beartype.beartype
    def get_map_data(self,datatype:Datatype):
        if self.data_states[datatype] != Datastate.COMPLETE:
            self.load_map_data(datatype)
        return self.data[datatype]

    @beartype.beartype
    def update_filenames_with_bounding_box(self,boundingbox_str:str):
        for i in Datatype:
            self.filenames[i] = self.filenames[i].replace("{BBOX_STR}",boundingbox_str)

    @beartype.beartype
    def update_filenames_with_projection(self,projection_str:str):
        for i in Datatype:
            self.filenames[i] = self.filenames[i].replace("{PROJ_STR}",projection_str)

    def calculate_bounding_box_and_projection(self):
        if self.filenames[Datatype.GEOLOGY] == None:
            return None, None
        temp_geology_filename = self.filenames[Datatype.GEOLOGY].replace("{BBOX_STR}","")
        temp_geology_filename = temp_geology_filename.replace("{PROJ_STR}","")
        try:
            geology_data = geopandas.read_file(temp_geology_filename)
            XYbounds = geology_data.total_bounds
            return {"minx":XYbounds[0],"maxx":XYbounds[1],"miny":XYbounds[2],"maxy":XYbounds[3]}, geology_data.crs
        except Exception:
            print("Could not open geology file", temp_geology_filename, "no bounding box or projection found")
            return None, None

    def merge_structure_with_geology(self,c_l):
        if self.data[Datatype.GEOLOGY] is not None and self.data[Datatype.STRUCTURE] is not None:
            columns = ['geometry', c_l['d'], c_l['dd'], c_l['sf'],c_l['bo']]
            if 'sl' in c_l:
                columns.append(c_l['sl'])
            columns = list(set(columns))
            structure_code = geopandas.sjoin(self.data[Datatype.STRUCTURE][columns], self.data[Datatype.GEOLOGY], how="left", predicate="within")
            if "index_left" in structure_code.columns:
                structure_code.drop(columns=["index_left"],inplace=True)
            if "index_right" in structure_code.columns:
                structure_code.drop(columns=["index_right"],inplace=True)
            self.data[Datatype.STRUCTURE] = structure_code[~structure_code[c_l['o']].isnull()]

    @beartype.beartype
    def export_wkt_format_files(self):
        # TODO: - Move away from tab seperators entirely (topology and map2model)

        # Check geology data status and export to a WKT format file
        if self.dirtyflags[Datatype.GEOLOGY] == True:
            self.load_map_data(Datatype.GEOLOGY)
        if self.data_states[Datatype.GEOLOGY] == Datastate.COMPLETE:
            geology = self.get_map_data(Datatype.GEOLOGY)
            hint_flag = False  # use GSWA strat database to provide topology hints
            columns = ['geometry', self.config.c_l['o'], self.config.c_l['c'], self.config.c_l['g'],
                self.config.c_l['u'], self.config.c_l['min'], self.config.c_l['max'], self.config.c_l['ds'],
                self.config.c_l['r1'], self.config.c_l['r2']]
            sub_geol = geology[columns]
            sub_geol.reset_index(inplace=True,drop=True)
            # print(sub_geol, self.config.geology_filename_wkt, self.config.c_l, hint_flag)
            Topology.save_geol_wkt(sub_geol, self.config.geology_filename_wkt, self.config.c_l, hint_flag, self.config.verbose_level)
        else:
            print(f"Cannot export geology data as it only loaded to {self.data_states[Datatype.GEOLOGY]} status")

        # Check mineral deposits data status and export to a WKT format file
        if self.dirtyflags[Datatype.MINERAL_DEPOSIT] == True:
            self.load_map_data(Datatype.MINERAL_DEPOSIT)
        if self.data_states[Datatype.MINERAL_DEPOSIT] == Datastate.COMPLETE:
            mindeps = self.get_map_data(Datatype.MINERAL_DEPOSIT)
            if mindeps is not None:
                columns = ['geometry', self.config.c_l['msc'], self.config.c_l['msn'], self.config.c_l['mst'],
                    self.config.c_l['mtc'], self.config.c_l['mscm'], self.config.c_l['mcom']]
                sub_mindep = mindeps[columns]
                Topology.save_mindep_wkt(sub_mindep, self.config.mindep_filename_wkt, self.config.c_l)
        elif self.get_filename(Datatype.MINERAL_DEPOSIT) == "":
            pass
        else:
            print(f"Cannot export mineral deposit data as it only loaded to {self.data_states[Datatype.MINERAL_DEPOSIT]} status")


        # Check structure data status and export to a WKT format file
        if self.dirtyflags[Datatype.STRUCTURE] == True:
            self.load_map_data(Datatype.STRUCTURE)
        if self.data_states[Datatype.STRUCTURE] == Datastate.COMPLETE:
            orientations = self.get_map_data(Datatype.STRUCTURE)
            columns = ['geometry', self.config.c_l['gi'], self.config.c_l['d'], self.config.c_l['dd']]
            sub_pts = orientations[columns]
            Topology.save_structure_wkt(sub_pts, self.config.structure_filename_wkt, self.config.c_l)
        else:
            print(f"Cannot export structure data as it only loaded to {self.data_states[Datatype.STRUCTURE]} status")

        # Check faults data status and export to a WKT format file
        if self.dirtyflags[Datatype.FAULT] == True:
            self.load_map_data(Datatype.FAULT)
        if self.data_states[Datatype.FAULT] == Datastate.COMPLETE:
            faults = self.get_map_data(Datatype.FAULT)
            columns = ['geometry', self.config.c_l['o'], self.config.c_l['f']]
            sub_lines = faults[columns]
            Topology.save_faults_wkt(sub_lines, self.config.fault_filename_wkt, self.config.c_l)
        else:
            print(f"Cannot export fault data as it only loaded to {self.data_states[Datatype.FAULT]} status")

    def load_dtm(self):
        source = self.filenames[Datatype.DTM]
        # group all Australian states codes under the global country code (ISO 3166 ALPHA-2)
        polygon_ll = self.config.polygon.to_crs(self.config.dtm_crs)
        minlong = polygon_ll.total_bounds[0] - self.config.step_out
        maxlong = polygon_ll.total_bounds[2] + self.config.step_out
        minlat = polygon_ll.total_bounds[1] - self.config.step_out
        maxlat = polygon_ll.total_bounds[3] + self.config.step_out
        if self.config.verbose_level != VerboseLevel.NONE:
        # if self.config.quiet != "All":
            print("Fetching DTM... ", end=" bbox:")
            print(minlong, maxlong, minlat, maxlat)
            print('source=', source)

        if source in ("WA", "NSW", "VIC", "SA", "QLD", "ACT", "TAS"):
            source = 'AU'

        success = False
        num_attempts = 10
        if source.startswith('http') or source == 'AU':
            i, done = 0, False
            while not done:
                if i >= num_attempts:
                    break
                try:
                    dtm = m2l_utils.load_and_reproject_dtm(self.config.polygon, self.working_projection, url=source)
                    done = True
                    success = True
                except Exception as e:
                    print(e)
                    time.sleep(1)
                    i += 1
        else:
            dtm = m2l_utils.load_and_reproject_dtm(self.config.polygon, self.working_projection, url=source)
            success = True
        if success == False:
            raise NameError(
                f'map2loop error: Could not access DTM server after {num_attempts} attempts'
            )

        self.data[Datatype.DTM] = dtm
        self.dirtyflags[Datatype.DTM] = False
        self.data_states[Datatype.DTM] = Datastate.COMPLETE
        if self.config.verbose_level == VerboseLevel.ALL:
            with dtm.open() as dtmtmp:
                plt.imshow(dtmtmp.read(1), cmap='terrain', vmin=numpy.percentile(dtmtmp.read(1), 5), vmax=numpy.percentile(dtmtmp.read(1), 95))
                plt.title('DTM Reprojected')
                plt.show()

    @beartype.beartype
    def calc_depth_grid(self, workflow:dict):
        dtm = self.get_map_data(Datatype.DTM).open()
        if (self.get_map_data(Datatype.DTB_GRID) == None):
            self.dtb = 0
            self.dtb_null = 0

            if self.config.verbose_level != VerboseLevel.NONE:
                print("dtb and dtb_null set to 0")
            return
        
        dtb_null=self.config.run_flags['dtb_null'] 
        cover_dip=self.config.run_flags['cover_dip'] 
        cover_spacing=self.config.run_flags['cover_spacing']
        # print(workflow['cover_map'],dtb_null,cover_dip,cover_spacing,self.get_filename(Datatype.COVER_MAP))
        # TODO: DTB need to be defined, every function call bellow here that has a False boolean is referencing to the workflow['cover_map'] flag
        # dtb grid clipped to cover map
        dtb_clip_filename = os.path.join(self.config.output_path, 'young_cover_grid_clip.tif')

        with fiona.open(self.get_filename(Datatype.COVER_MAP), "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]

        with rasterio.open(self.dtb_grid_file) as src:
            out_image, out_transform = rasterio.mask.mask(src,
                                                          shapes,
                                                          crop=True)
            out_meta = src.meta.copy()

        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })

        with rasterio.open(dtb_clip_filename, "w", **out_meta) as dest:
            dest.write(out_image)

        dtb = rasterio.open(dtb_clip_filename)
        self.dtb=dtb
        self.dtb_null=dtb_null

        m2l_geometry.process_cover(self.config, self, workflow)

    @beartype.beartype
    def export_orientations(self, workflow:dict):
        m2l_geometry.save_orientations(self.config, self, workflow)

        if self.config.verbose_level == VerboseLevel.ALL:
            m2l_utils.plot_points(
                os.path.join(self.config.output_path, 'orientations.csv'),
                self.get_map_data(Datatype.GEOLOGY), 'formation', 'X', 'Y', False, 'alpha', 'Orientations')

        # Create arbitrary points for series without orientation data
        m2l_geometry.create_orientations(self.config, self, workflow)

    @beartype.beartype
    def export_contacts(self, workflow:dict):

        self.basal_contacts = m2l_geometry.save_basal_contacts(self.config, self, workflow)

        # Remove basal contacts defined by faults, no decimation
        self.basal_contacts_no_faults = m2l_geometry.save_basal_no_faults(self.config, self)

        # TODO: Get project run to use self.basal_contacts_no_fault instead of via basal_contacts.shp file
        # Remove faults from decimated basal contacts then save
        m2l_geometry.save_basal_contacts_csv(self.basal_contacts_no_faults, self.config, self, workflow)
        
        if self.config.verbose_level == VerboseLevel.ALL:
            m2l_utils.plot_points(
                os.path.join(self.config.output_path, 'contacts4.csv'),
                self.get_map_data(Datatype.GEOLOGY), 'formation', 'X', 'Y', False, 'alpha',"Contacts")

    @beartype.beartype
    def export_faults(self, workflow:dict, dip_grid, dip_dir_grid):

        m2l_geometry.save_faults(self.config, self, workflow)

        faults = self.get_map_data(Datatype.FAULT)
        if type(faults) != None and len(faults) > 0:
            m2l_interpolation.process_fault_throw_and_near_faults_from_grid(self.config, self, workflow, dip_grid, dip_dir_grid)