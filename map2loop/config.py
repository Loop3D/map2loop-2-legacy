import os
import random
import urllib.request
import hjson
import pandas as pd
import matplotlib.colors as colors
import beartype
import warnings
from .m2l_enums import Datatype, VerboseLevel

class Config(object):
    """Object that represents a sub-project. It is defined by some source data, 
    a region of interest (bounding box or polygon) and some execution flags.
    """

    def __init__(self,map_data,verbose_level=VerboseLevel.ALL):
        self.project_path = None
        self.bbox_3d = {"minx":0,
                        "maxx":1000,
                        "miny":0,
                        "maxy":1000,
                        "base":-10000,
                        "top" :1200}
        self.polygon = None
        self.step_out = 0.1
        self.dtm_crs = "EPSG:4326"
        self.project_crs = None
        
        # TODO: Need to add check for map_data and verbose_level types
        self.map_data = map_data
        self.verbose_level = verbose_level

        # Set default run flags
        self.run_flags = {
            'aus': True,
            'deposits': "Fe,Cu,Au,NONE",
            'dtb': '',
            'dtb_null':-2147483648,
            'cover_dip':10,
            'cover_spacing':5000,
            'orientation_decimate': 1,
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
            'min_pluton_area':10000000,
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
            'fault_formation_weight':5,
            'map2graph':False,
            'granular_map2graph':False
        }

    @beartype.beartype
    def update(self,
                 project_path:str,
                 bbox_3d,
                 polygon,
                 step_out,
                 dtm_crs,
                 project_crs,
                 local,
                 metadata_filename:str="",
                 clut_path:str="",
                 run_flags=dict(),
                 ):

        self.project_path = project_path
        self.graph_path = os.path.join(project_path, 'graph')
        self.tmp_path = os.path.join(project_path, 'tmp')
        self.output_path = os.path.join(project_path, 'output')

        # filenames for WKT format (map2graph requirement)
        self.fault_filename_wkt = os.path.join(project_path,'tmp',"faults.csv")
        # self.fault_output_file_csv = os.path.join(project_path,'output',"faults.csv")
        self.structure_filename_wkt = os.path.join(project_path,'tmp', "structure.csv")
        self.geology_filename_wkt = os.path.join(project_path,'tmp', "geology.csv")
        self.mindep_filename_wkt = os.path.join(project_path,'tmp', "mindep.csv")

        self.strat_graph_filename = os.path.join(project_path,"graph", "graph_strat_NONE.gml")
        # self.dtm_filename = os.path.join(project_path,"dtm", 'dtm.tif')
        # self.dtm_reproj_filename = os.path.join(project_path,"dtm", 'dtm_rp.tif')
        self.clut_path = clut_path

        if isinstance(run_flags, dict):
            for key in run_flags:
                if not key in self.run_flags:
                    warnings.warn(f"config run_flags key {key} is not a valid option")
            self.run_flags.update(run_flags)
            # TODO: Add sanity checks for run_flags updated from user
            if self.run_flags['interpolation_spacing'] < 0:
                self.run_flags['interpolation_spacing'] = - self.run_flags['interpolation_spacing']
            self.run_flags['orientation_decimate'] = max(self.run_flags['orientation_decimate'],1)
            self.run_flags['contact_decimate'] = max(self.run_flags['contact_decimate'],1)
            self.run_flags['fault_decimate'] = max(self.run_flags['fault_decimate'],1)
            self.run_flags['contact_orientation_decimate'] = max(self.run_flags['contact_orientation_decimate'],1)
            self.run_flags['fold_decimate'] = max(self.run_flags['fold_decimate'],1)
        else:
            print('run_flags must be a dictionary, setting config run flags to the defaults.')

        self.bbox_3d = bbox_3d
        self.bbox = tuple([bbox_3d["minx"], bbox_3d["miny"], bbox_3d["maxx"], bbox_3d["maxy"]])

        # Recalculate height and width from bbox and interpolation_spacing
        self.width = int((self.bbox[2] - self.bbox[0]) / self.run_flags['interpolation_spacing']) + 1
        self.height = int((self.bbox[3] - self.bbox[1]) / self.run_flags['interpolation_spacing']) + 1
        self.polygon = polygon
        self.step_out = step_out

        self.dtm_crs = dtm_crs
        self.project_crs = project_crs
        self.local = local

        self.read_metadata(metadata_filename)
        self.create_cmap()

    @beartype.beartype
    def read_metadata(self, filename:str):
        """
        Helper function that turns json and hjson files into usable configuration dictionaries.

        Args:
            filename (string) : Path or url to metadata file
        """
        if (filename == "" and self.map_data.get_map_data(Datatype.GEOLOGY) is not None):
            metadata_error = "When using your own local files or remote data, pass the path or url to a valid metadata file (.json or .hjson) that describes your input column names.\n"
            metadata_error += "You can find an example config here https://gist.github.com/yohanderose/a127c29cb88529f049a5bafc881bb1a0"
            raise ValueError(metadata_error)
        
        # Load in dictionary that describes the column names
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
        except Exception as e:
            self.errorState = str(e)

    # TODO: Make it more obvious when formations are not found and random colours are assigned
    def create_cmap(self):
        # Make colours consistent from map to model
        formations = sorted([
            formation.replace(" ", "_").replace('-', '_') for formation in
            list(set(self.map_data.get_map_data(Datatype.GEOLOGY)['UNIT_NAME'].to_numpy()))
        ])
        temp_colours = [""] * len(formations)
        self.colour_dict = dict(zip(formations, temp_colours))
        try:
            # Try to retrieve the clut reference
            colour_ref = pd.read_csv(self.clut_path)
            for formation in formations:
                key = formation.replace("-","_")
                colour = None
                try:
                    colour = colour_ref[colour_ref['code'] == key]['colour'].to_numpy()[0]
                except Exception as e:
                    try:
                        colour = colour_ref[colour_ref['UNITNAME'] == key]['colour'].to_numpy()[0]
                    except Exception as e:
                        print('formation ' + formation + ' does not have a colour reference')
                        colour = ('#%02X%02X%02X' %
                                  (random.randint(0, 255), random.randint(
                                      0, 255), random.randint(0, 255)))

                self.colour_dict[formation] = colour
                #print(key, colour)

        except Exception as e:
            # Otherwise, just append a random set
            self.clut_path = ""
            random_colours = [
                '#%02X%02X%02X' % (random.randint(
                    0, 255), random.randint(0, 255), random.randint(0, 255))
                for i in range(len(formations))
            ]
            i = 0
            for key in self.colour_dict.keys():
                self.colour_dict[key] = random_colours[i]
        
        self.cmap = colors.ListedColormap(self.colour_dict.values(),
                                          name='geol_key')

    @beartype.beartype
    def save_cmap(self,workflow:dict):
        """Create a colourmap for the model using the colour code
        """
        all_sorts = pd.read_csv(os.path.join(self.tmp_path, 'all_sorts_clean.csv'))

        if(workflow['cover_map']):
            self.colour_dict['cover']='#FFFFC0'
            self.colour_dict['cover_up']='#FFFFC0'
        colours = []
        newKeys = [key.replace("-","_") for key in self.colour_dict.keys()]
        newColourDict = dict(zip(newKeys,list(self.colour_dict.values())))
        for code in all_sorts['code']:
            # colours.append([self.colour_dict[code]])
            code = code.replace("-","_")
            colours.append([newColourDict[code]])
        data = colours
        expected_extra_cols = pd.DataFrame(columns=['colour'], data=data)
        all_sorts = pd.concat([all_sorts, expected_extra_cols], axis=1)
        all_sorts.to_csv(os.path.join(self.tmp_path, 'all_sorts_clean.csv'), ",", index=None)
