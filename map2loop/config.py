import os
import sys
import time
import shutil
import re
import csv

import numpy as np
import pandas as pd
import geopandas as gpd
from map2loop import _clut_path
from map2loop.topology import Topology
from map2loop import m2l_utils
from map2loop import m2l_geometry
from map2loop import m2l_interpolation
from map2loop import m2l_map_checker
from map2loop.m2l_utils import display, enable_quiet_mode, disable_quiet_mode, print
from map2loop.m2l_export import export_to_projectfile
import map2model

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import rasterio
import shapely


class Config(object):
    """Object that represents a sub-project. It is defined by some source data, 
    a region of interest (bounding box or polygon) and some execution flags.
    """

    def __init__(self,
                 project_path,
                 overwrite,
                 geology_file,
                 fault_file,
                 fold_file,
                 structure_file,
                 mindep_file,
                 bbox_3d,
                 polygon,
                 step_out,
                 dtm_crs,
                 proj_crs,
                 local,
                 quiet,
                 c_l={},
                 **kwargs
                 ):
        """Creates the config object.

        :param geology_file: Path to local file or remote database containing stratigraphic information.
        :type geology_file: string
        :param fault_file: Path to local file or remote database containing fault information.
        :type fault_file: string
        :param fold_file: Path to local file or remote database containing fold information.
        :type fold_file: string
        :param structure_file: Path to local file or remote database containing structure information.
        :type structure_file: string
        :param mindep_file: Path to local file or remote database containing mindep information.
        :type mindep_file: string
        :param bbox_3d: Bounding box region of interest with x, y and z bounds.
        :type bbox_3d: dict
        :param polygon: Region of interest dataframe for creating figures and preprocessing.
        :type polygon: geopandas.GeoDataFrame
        :param step_out: Buffer to ensure DTM covers entire expected area on conversion.
        :type step_out: float
        :param dtm_crs: Projection the DTM is in.
        :type dtm_crs: dict
        :param proj_crs: Projection the project sources are in.
        :type proj_crs: dict
        :param local: Flag to indicate if remote or local data is being used.
        :type local: bool
        :param c_l: Map of configuration flags that dictate the flow of execution, defaults to {}
        :type c_l: dict, optional
        """
        self.clut_path = "https://gist.githubusercontent.com/yohanderose/8f7e2d57db9086fbe1a7c651b9e25996/raw/ac5062e68d251c21bbc24b811ee5b17cc2f98ce3/500kibg_colours.csv"
        # Create directory structure if possible
        if (not os.path.exists(project_path)):
            os.mkdir(project_path)
        else:
            if overwrite:
                shutil.rmtree(project_path)
                time.sleep(1)
                os.mkdir(project_path)
            else:
                allow = input(
                    "Directory \"{}\" exists, overwrite? (y/[n])".format(project_path))
                if allow == "y":
                    shutil.rmtree(project_path)
                    os.mkdir(project_path)
                else:
                    sys.exit(0)

        self.project_path = project_path

        self.graph_path = self.project_path+'/graph/'
        self.tmp_path = self.project_path+'/tmp/'
        self.data_path = self.project_path+'/data/'
        self.dtm_path = self.project_path+'/dtm/'
        self.output_path = self.project_path+'/output/'
        self.vtk_path = self.project_path+'/vtk/'

        self.fault_file_csv = self.tmp_path + "faults.csv"
        self.structure_file_csv = self.tmp_path + "structure.csv"
        self.geology_file_csv = self.tmp_path + "geology.csv"
        self.mindep_file_csv = self.tmp_path + "mindep.csv"

        self.strat_graph_file = self.graph_path + "graph_strat_NONE.gml"
        self.dtm_file = self.dtm_path+'dtm.tif'
        self.dtm_reproj_file = self.dtm_path+'dtm_rp.tif'

        if(not os.path.isdir(self.tmp_path)):
            os.mkdir(self.tmp_path)
        if(not os.path.isdir(self.data_path)):
            os.mkdir(self.data_path)
        if(not os.path.isdir(self.output_path)):
            os.mkdir(self.output_path)
        if(not os.path.isdir(self.dtm_path)):
            os.mkdir(self.dtm_path)
        if(not os.path.isdir(self.vtk_path)):
            os.mkdir(self.vtk_path)
        if(not os.path.isdir(self.graph_path)):
            os.mkdir(self.graph_path)

        if quiet:
            enable_quiet_mode()

        self.bbox_3d = bbox_3d
        self.bbox = tuple([bbox_3d["minx"], bbox_3d["miny"],
                           bbox_3d["maxx"], bbox_3d["maxy"]])
        self.polygon = polygon
        self.step_out = step_out
        if quiet:
            plt.ioff()
        self.quiet = quiet
        self.c_l = c_l

        self.dtm_crs = dtm_crs
        self.proj_crs = proj_crs

        # Check input maps for missing values
        drift_prefix = kwargs.get('drift_prefix', ['None'])
        self.local = local
        #       - Check if fold file is always the same as fault or needs to be seperated
        # TODO: Allow for input as a polygon, not just a bounding box.
        structure_file, geology_file, fault_file, mindep_file, fold_file, c_l = m2l_map_checker.check_map(
            structure_file, geology_file, fault_file, mindep_file, fold_file, self.tmp_path, self.bbox, c_l, proj_crs, self.local, drift_prefix)

        # Process and store workflow params
        self.geology_file = geology_file
        self.structure_file = structure_file
        self.fault_file = fault_file
        self.fold_file = fold_file
        self.mindep_file = mindep_file

        disable_quiet_mode()

    def preprocess(self):
        """[summary]

        :param command: [description], defaults to ""
        :type command: str, optional
        """

        if self.quiet:
            enable_quiet_mode()

        geology = gpd.read_file(self.geology_file, bbox=self.bbox)
        geology[self.c_l['g']].fillna(geology[self.c_l['g2']], inplace=True)
        geology[self.c_l['g']].fillna(geology[self.c_l['c']], inplace=True)
        faults = gpd.read_file(self.fault_file, bbox=self.bbox)
        folds = gpd.read_file(self.fold_file, bbox=self.bbox)
        structures = gpd.read_file(self.structure_file, bbox=self.bbox)
        mindeps = gpd.read_file(self.mindep_file, bbox=self.bbox)

        # Fix crs to project default and overwrite source
        # TODO: - Maybe do these checks in map_checker
        geology.crs = self.proj_crs
        faults.crs = self.proj_crs
        folds.crs = self.proj_crs
        structures.crs = self.proj_crs
        mindeps.crs = self.proj_crs

        geology.to_file(self.geology_file)
        faults.to_file(self.fault_file)
        folds.to_file(self.fold_file)
        if len(mindeps) > 0:
            mindeps.to_file(self.mindep_file)

        # Check if bedding data uses the strike convention instead of dip direction
        if(self.c_l['otype'] == 'strike'):
            structures['azimuth2'] = structures.apply(
                lambda row: row[self.c_l['dd']]+90.0, axis=1)
            self.c_l['dd'] = 'azimuth2'
            self.c_l['otype'] = 'dip direction'
            structures.to_file(self.structure_file)

        self.geology = geology
        self.faults = faults
        self.structures = structures
        self.mindeps = mindeps

        try:
            self.geology_figure = geology.copy().plot(column=self.c_l['c'], figsize=(
                10, 10), edgecolor='#000000', linewidth=0.2).get_figure()

            base = geology.plot(column=self.c_l['c'], figsize=(
                10, 10), edgecolor='#000000', linewidth=0.2, legend=True)
            leg = base.get_legend()
            leg.set_bbox_to_anchor((1.04, 1))

            structures.plot(ax=base, color='none', edgecolor='black')

            faults.plot(ax=base, cmap='rainbow',
                        column=self.c_l['f'], figsize=(10, 10), linewidth=0.4)
            structures[['geometry', self.c_l['gi'],
                        self.c_l['d'], self.c_l['dd']]].plot(ax=base)

            fig = self.polygon.plot(ax=base, color='none', edgecolor='black').set_title(
                "Input {}".format(self.bbox)).get_figure()
            fig.savefig(self.tmp_path+"/input-data.png")

            print("Input graphic saved to: " +
                  self.tmp_path + "input-fig.png")

            if not self.quiet:
                plt.show()

            return
        except Exception as e:
            print(e)

        disable_quiet_mode()

    def export_csv(self):
        # TODO: - Move away from tab seperators entirely (topology and map2model)

        # Save geology polygons
        hint_flag = False  # use GSWA strat database to provide topology hints
        sub_geol = self.geology[['geometry', self.c_l['o'], self.c_l['c'], self.c_l['g'],
                                 self.c_l['u'], self.c_l['min'], self.c_l['max'], self.c_l['ds'], self.c_l['r1'], self.c_l['r2']]]
        Topology.save_geol_wkt(
            sub_geol, self.geology_file_csv, self.c_l, hint_flag)

        # Save mineral deposits
        sub_mindep = self.mindeps[['geometry', self.c_l['msc'], self.c_l['msn'],
                                   self.c_l['mst'], self.c_l['mtc'], self.c_l['mscm'], self.c_l['mcom']]]
        Topology.save_mindep_wkt(
            sub_mindep, self.mindep_file_csv, self.c_l)

        # Save orientation data
        sub_pts = self.structures[[
            'geometry', self.c_l['gi'], self.c_l['d'], self.c_l['dd']]]
        Topology.save_structure_wkt(
            sub_pts, self.structure_file_csv, self.c_l)

        # Save faults
        sub_lines = self.faults[['geometry', self.c_l['o'], self.c_l['f']]]
        Topology.save_faults_wkt(sub_lines, self.fault_file_csv, self.c_l)

    def update_parfile(self):
        Topology.save_parfile(self, self.c_l, self.output_path, self.geology_file_csv, self.fault_file_csv, self.structure_file_csv,
                              self.mindep_file_csv, self.bbox[0], self.bbox[1], self.bbox[2], self.bbox[3], 500.0, 'Fe,Cu,Au,NONE')

    def run_map2model(self, deposits, aus):
        run_log = map2model.run(self.graph_path, self.geology_file_csv,
                                self.fault_file_csv, self.mindep_file_csv,
                                self.bbox_3d,
                                self.c_l,
                                self.quiet,
                                deposits
                                )

        print(run_log)

        print("Resolving ambiguities using ASUD...", end='\toutput_dir:')
        if aus:
            Topology.use_asud(self.strat_graph_file,  self.graph_path)
            self.strat_graph_file = self.graph_path+'ASUD_strat.gml'
        print("Done.")

        print("Generating topology graph display and unit groups...")
        self.G = nx.read_gml(self.strat_graph_file, label='id')
        selected_nodes = [n for n, v in self.G.nodes(data=True) if n >= 0]

        if not self.quiet:
            nx.draw_networkx(self.G, pos=nx.kamada_kawai_layout(
                self.G), arrows=True, nodelist=selected_nodes)

        nlist = list(self.G.nodes.data('LabelGraphics'))
        nlist.sort()
        for node in nlist:
            if node[0] >= 0:
                elem = str(node[1]).replace(
                    "{'text':", "").replace(", 'fontSize': 14}", "")
                # second = elem.split(":").replace("'", "")
                print(node[0], " ", elem)

        # plt.savefig(self.tmp_path+"topology-fig.png")
        print("Topology figure saved to", self.tmp_path+"topology-fig.png")

        # Save groups of stratigraphic units
        groups, self.glabels, G = Topology.get_series(
            self.strat_graph_file, 'id')
        Topology.save_units(G, self.tmp_path, self.glabels,
                            Australia=True, asud_strat_file="https://gist.githubusercontent.com/yohanderose/3b257dc768fafe5aaf70e64ae55e4c42/raw/8598c7563c1eea5c0cd1080f2c418dc975cc5433/ASUD.csv",
                            quiet=self.quiet)

        print("Done")

    def load_dtm(self):

        polygon_ll = self.polygon.to_crs(self.dtm_crs)

        minlong = polygon_ll.total_bounds[0]-self.step_out
        maxlong = polygon_ll.total_bounds[2]+self.step_out
        minlat = polygon_ll.total_bounds[1]-self.step_out
        maxlat = polygon_ll.total_bounds[3]+self.step_out

        print("Fetching DTM... ", end=" bbox:")
        print(minlong, maxlong, minlat, maxlat)
        downloaded = False
        i = 0
        print('Attempt: 0 ', end='')
        while downloaded == False:
            try:
                m2l_utils.get_dtm(self.dtm_file, minlong,
                                  maxlong, minlat, maxlat)
                downloaded = True
            except:
                time.sleep(10)
                i = i+1
                print(' ', i, end='')
        if(i == 100):
            raise NameError(
                'map2loop error: Could not access DTM server after 100 attempts')
        print('Done.')

        geom_rp = m2l_utils.reproject_dtm(
            self.dtm_file, self.dtm_reproj_file, self.dtm_crs, self.proj_crs)
        self.dtm = rasterio.open(self.dtm_reproj_file)

        if not self.quiet:
            plt.imshow(self.dtm.read(1), cmap='terrain',
                       vmin=0, vmax=1000)

            plt.title('DTM')
            plt.show()

    def join_features(self):
        # Faults
        self.faults_clip = gpd.read_file(self.fault_file, bbox=self.bbox)
        self.faults_clip.crs = self.proj_crs
        self.faults_clip_file = self.tmp_path + "faults_clip.shp"
        self.faults_clip.to_file(self.faults_clip_file)

        # Geology
        self.geol_clip = m2l_utils.explode(self.geology)
        self.geol_clip.crs = self.proj_crs
        self.geol_clip_file = self.tmp_path + "geol_clip.shp"
        self.geol_clip.to_file(self.geol_clip_file)

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)

        # Structures
        list1 = ['geometry', self.c_l['d'],
                 self.c_l['dd'], self.c_l['sf'], self.c_l['bo']]
        list2 = list(set(list1))
        sub_pts = self.structures[list2]
        structure_code = gpd.sjoin(
            sub_pts, self.geol_clip, how="left", op="within")

        minx, miny, maxx, maxy = self.bbox
        y_point_list = [miny, miny, maxy, maxy, miny]
        x_point_list = [minx, maxx, maxx, minx, minx]

        bbox_geom = shapely.geometry.Polygon(zip(x_point_list, y_point_list))

        # TODO: 'polygo' is never used
        polygo = gpd.GeoDataFrame(
            index=[0], crs=self.proj_crs, geometry=[bbox_geom])
        is_bed = structure_code[self.c_l['sf']].str.contains(
            self.c_l['bedding'], regex=False)

        structure_clip = structure_code[is_bed]
        structure_clip.crs = self.proj_crs

        if(self.c_l['otype'] == 'strike'):
            structure_clip['azimuth2'] = structure_clip.apply(
                lambda row: row[self.c_l['dd']]+90.0, axis=1)
            self.c_l['dd'] = 'azimuth2'
            self.c_l['otype'] = 'dip direction'

        self.structure_clip = structure_clip[~structure_clip[self.c_l['o']].isnull(
        )]
        self.structure_clip_file = self.tmp_path+'structure_clip.shp'
        self.structure_clip.to_file(self.structure_clip_file)

        # Save geology clips
        Topology.save_group(Topology, self.G, self.tmp_path,
                            self.glabels, self.geol_clip, self.c_l, self.quiet)

    def calc_depth_grid(self, dtb):
        dtm = self.dtm

        if dtb == "":
            self.dtb = 0
            self.dtb_null = 0

            print("dtb and dtb_null set to 0")
            return

        # TODO: DTB need to be defined, every function call bellow here that has a False boolean is referencing to the workflow['cover_map'] flag
        # dtb_grid = data_path+'young_cover_grid.tif' #obviously hard-wired for the moment
        # dtb_null = '-2147483648' #obviously hard-wired for the moment
        # cover_map_path = data_path+'Young_Cover_FDS_MGA_clean.shp' #obviously hard-wired for the moment
        # dtb_clip = output_path+'young_cover_grid_clip.tif' #obviously hard-wired for the moment
        # cover_dip = 10 # dip of cover away from contact
        # cover_spacing = 5000 # of contact grid in metres

        dtb_raw = rasterio.open(dtb_grid)

        cover = gpd.read_file(cover_map_path)

        with fiona.open(cover_map_path, "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]

        with rasterio.open(dtb_grid) as src:
            out_image, out_transform = rasterio.mask.mask(
                src, shapes, crop=True)
            out_meta = src.meta.copy()

        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})

        with rasterio.open(dtb_clip, "w", **out_meta) as dest:
            dest.write(out_image)

        dtb = rasterio.open(dtb_clip)

        m2l_geometry.process_cover(output_path, dtm, dtb, dtb_null, cover,
                                   workflow['cover_map'], cover_dip, bbox, proj_crs, cover_spacing, contact_decimate=3, use_vector=True, use_grid=True)

    def export_orientations(self, orientation_decimate):
        m2l_geometry.save_orientations(
            self.structure_clip, self.output_path, self.c_l, orientation_decimate, self.dtm, self.dtb, self.dtb_null, False)

        if not self.quiet:
            m2l_utils.plot_points(self.output_path+'orientations.csv',
                                  self.geol_clip, 'formation', 'X', 'Y', False, 'alpha')

        # Create arbitrary points for series without orientation data
        m2l_geometry.create_orientations(
            self.tmp_path, self.output_path, self.dtm, self.dtb, self.dtb_null, False, self.geol_clip, self.structure_clip, self.c_l)

    def export_contacts(self, contact_decimate, intrusion_mode):

        ls_dict, ls_dict_decimate = m2l_geometry.save_basal_contacts(
            self.tmp_path, self.dtm, self.dtb, self.dtb_null, False, self.geol_clip, contact_decimate, self.c_l, intrusion_mode)

        # Remove basal contacts defined by faults, no decimation
        m2l_geometry.save_basal_no_faults(
            self.tmp_path+'basal_contacts.shp', self.tmp_path+'faults_clip.shp', ls_dict, 10, self.c_l, self.proj_crs)

        # Remove faults from decimated basal contacts then save
        contacts = gpd.read_file(self.tmp_path+'basal_contacts.shp')
        m2l_geometry.save_basal_contacts_csv(
            contacts, self.output_path, self.dtm, self.dtb, self.dtb_null, False, contact_decimate, self.c_l)
        # False in this call was already false and isn't the cover flag
        if not self.quiet:
            m2l_utils.plot_points(self.output_path+'contacts4.csv',
                                  self.geol_clip, 'formation', 'X', 'Y', False, 'alpha')

    # Interpolates a regular grid of orientations from an shapefile of
    # arbitrarily-located points and saves out four csv files of l, m & n
    # direction cosines and dip dip direction data
    def test_interpolation(self, interpolation_spacing, misorientation, interpolation_scheme):

        geology_file = self.geol_clip_file
        structure_file = self.structure_clip_file
        basal_contacts = self.tmp_path+'basal_contacts.shp'
        self.spacing = interpolation_spacing  # grid spacing in meters
        # misorientation = misorientation
        self.scheme = interpolation_scheme
        orientations = self.structures

        group_girdle = m2l_utils.plot_bedding_stereonets(
            orientations, self.geology, self.c_l, self.quiet)
        super_groups, self.use_gcode3 = Topology.super_groups_and_groups(
            group_girdle, self.tmp_path, misorientation)
        # print(super_groups)
        # print(self.geology['GROUP_'].unique())
        bbox = self.bbox

        orientation_interp, contact_interp, combo_interp = m2l_interpolation.interpolation_grids(
            geology_file, structure_file, basal_contacts, bbox, self.spacing, self.proj_crs, self.scheme, super_groups, self.c_l)

        with open(self.tmp_path+'interpolated_orientations.csv', 'w') as f:
            f.write('X, Y, l, m, n, dip, dip_dir\n')
            for row in orientation_interp:
                ostr = '{}, {}, {}, {}, {}, {}, {}\n'.format(
                    row[0], row[1], row[2], row[3], row[4], row[5], row[6])
                f.write(ostr)
        with open(self.tmp_path+'interpolated_contacts.csv', 'w') as f:
            f.write('X, Y, l, m, angle\n')
            for row in contact_interp:
                ostr = '{}, {}, {}, {}, {}\n'.format(
                    row[0], row[1], row[2], row[3], row[4])
                f.write(ostr)
        with open(self.tmp_path+'interpolated_combined.csv', 'w') as f:
            f.write('X, Y, l, m, n, dip, dip_dir\n')
            for row in combo_interp:
                ostr = '{}, {}, {}, {}, {}, {}, {}\n'.format(
                    row[0], row[1], row[2], row[3], row[4], row[5], row[6])
                f.write(ostr)

        if(self.spacing < 0):
            self.spacing = -(bbox[2]-bbox[0])/spacing
        self.x = int((bbox[2]-bbox[0])/self.spacing)+1
        self.y = int((bbox[3]-bbox[1])/self.spacing)+1
        x = self.x
        y = self.y
        print(x, y)
        dip_grid = np.ones((y, x))
        dip_grid = dip_grid*-999
        dip_dir_grid = np.ones((y, x))
        dip_dir_grid = dip_dir_grid*-999
        contact_grid = np.ones((y, x))
        contact_grid = dip_dir_grid*-999
        for row in combo_interp:
            r = int((row[1]-bbox[1])/self.spacing)
            c = int((row[0]-bbox[0])/self.spacing)
            dip_grid[r, c] = float(row[5])
            dip_dir_grid[r, c] = float(row[6])

        for row in contact_interp:
            r = int((row[1]-bbox[1])/self.spacing)
            c = int((row[0]-bbox[0])/self.spacing)
            contact_grid[r, c] = float(row[4])

        self.dip_grid = dip_grid
        self.dip_dir_grid = dip_dir_grid

        if not self.quiet:
            print('interpolated dips')
            plt.imshow(self.dip_grid, cmap="hsv",
                       origin='lower', vmin=-90, vmax=90)
            plt.show()

            print('interpolated dip directions')
            plt.imshow(self.dip_dir_grid, cmap="hsv",
                       origin='lower', vmin=0, vmax=360)
            plt.show()

            print('interpolated contacts')
            plt.imshow(contact_grid, cmap="hsv",
                       origin='lower', vmin=-360, vmax=360)
            plt.show()

    def create_map_cmap(self):
        """Create a colourmap for the model using the colour code
        """
        asc = pd.read_csv(self.tmp_path+'all_sorts_clean.csv', ",")
        # display(asc)
        colours = pd.read_csv(self.clut_path, ",")
        if self.c_l['c'] == 'CODE':
            code = self.c_l['c'].lower()
        else:
            code = self.c_l['c']

        colours = []  # container for the discrete colours we are using
        i = 0
        # initialise a colour index attribute column
        self.geol_clip['colour_index'] = np.nan
        for ind, strat in asc.iterrows():
            self.geol_clip[self.c_l['c']].str.replace(" ", "_")
            self.geol_clip.loc[self.geol_clip[self.c_l['c']] ==
                               strat['code'].replace("_", " "), 'colour_index'] = i
            colours.append(strat['colour'])
            i = i+1
        self.cmap = colors.ListedColormap(colours)

    def export_faults(self, fault_decimate, min_fault_length, fault_dip):
        # fault_decimate = 5
        # min_fault_length = 5000
        # fault_dip = 90

        m2l_geometry.save_faults(
            self.tmp_path+'faults_clip.shp', self.output_path, self.dtm, self.dtb, self.dtb_null, False, self.c_l, fault_decimate, min_fault_length, fault_dip)

        faults = pd.read_csv(self.fault_file_csv, sep='\t')
        faults_len = len(faults)
        if(faults_len > 0):
            m2l_interpolation.process_fault_throw_and_near_faults_from_grid(self.tmp_path, self.output_path, self.dtm_reproj_file, self.dtb, self.dtb_null, False, self.c_l, self.proj_crs, self.bbox,
                                                                            self.scheme, self.dip_grid, self.dip_dir_grid, self.x, self.y, self.spacing)

    def process_plutons(self, pluton_dip, pluton_form, dist_buffer, contact_decimate):
        # pluton_dip = 45
        pluton_dip = str(pluton_dip)
        self.pluton_form = pluton_form  # 'domes'

        # dist_buffer = 10
        # contact_decimate = 5  # store every nth contact point (in object order)

        m2l_geometry.process_plutons(self.tmp_path, self.output_path, self.geol_clip, self.local,
                                     self.dtm, self.dtb, self.dtb_null, False, self.pluton_form, pluton_dip, contact_decimate, self.c_l)

    def extract_section_features(self, seismic_line_file, seismic_bbox_file, seismic_interp_file):
        # Extract faults and basal contacts of groups from seismic section
        # input geology file (if local)

        seismic_line_file = seismic_line_file
        seismic_line = gpd.read_file(seismic_line_file)  # import map
        seismic_line.plot(figsize=(10, 10), edgecolor='#000000',
                          linewidth=0.2)  # display map
        display(seismic_line)

        # input geology file (if local)
        seismic_bbox_file = seismic_bbox_file
        seismic_bbox = gpd.read_file(seismic_bbox_file)  # import map
        seismic_bbox.set_index('POSITION', inplace=True)

        # input geology file (if local)
        seismic_interp_file = seismic_interp_file
        seismic_interp = gpd.read_file(seismic_interp_file)  # import map
        seismic_interp.plot(column='FEATURE', figsize=(
            10, 10), edgecolor='#000000', linewidth=0.5)  # display map
        display(seismic_interp)

        surface_cut = 2000

        m2l_geometry.extract_section(self.tmp_path, self.output_path, seismic_line, seismic_bbox,
                                     seismic_interp, self.dtm, self.dtb, self.dtb_null, False, surface_cut)

        contacts = pd.read_csv(self.output_path+'contacts4.csv', ", ")
        seismic_contacts = pd.read_csv(
            self.output_path+'seismic_base.csv', ", ")
        all_contacts = pd.concat([contacts, seismic_contacts], sort=False)
        all_contacts.to_csv(self.output_path+'contacts4.csv',
                            index=None, header=True)

        faults = pd.read_csv(self.output_path+'faults.csv', ", ")
        seismic_faults = pd.read_csv(
            self.output_path+'seismic_faults.csv', ", ")
        all_faults = pd.concat([faults, seismic_faults], sort=False)
        all_faults.to_csv(self.output_path+'faults.csv',
                          index=None, header=True)

    def propagate_contact_dips(self, contact_dip, contact_orientation_decimate):
        print("Propagating dips along contacts...")
        orientations = pd.read_csv(self.output_path+'orientations.csv', ", ")
        # This is supposed to be a csv but my csv doesn't have a geometry part
        contacts = gpd.read_file(self.tmp_path + 'basal_contacts.shp')
        # contact_dip = -999
        # contact_orientation_decimate = 5
        m2l_geometry.save_basal_contacts_orientations_csv(contacts, orientations, self.geol_clip, self.tmp_path, self.output_path, self.dtm, self.dtb,
                                                          self.dtb_null, False, contact_orientation_decimate, self.c_l, contact_dip, self.dip_grid, self.spacing, self.bbox)

    def calc_thickness(self, contact_decimate, null_scheme, thickness_buffer, max_thickness_allowed):
        # Estimate formation thickness and normalised formation thickness
        geology_file = self.tmp_path+'basal_contacts.shp'
        # contact_decimate = 5
        # null_scheme = 'null'
        m2l_interpolation.save_contact_vectors(
            geology_file, self.tmp_path, self.dtm, self.dtb, self.dtb_null, False, self.bbox, self.c_l, null_scheme, contact_decimate)

        # buffer = 5000
        # max_thickness_allowed = 10000

        m2l_geometry.calc_thickness_with_grid(self.tmp_path, self.output_path, thickness_buffer, max_thickness_allowed,
                                              self.c_l, self.bbox, self.dip_grid, self.dip_dir_grid, self.x, self.y, self.spacing)

        m2l_geometry.calc_min_thickness_with_grid(self.tmp_path, self.output_path, thickness_buffer, max_thickness_allowed,
                                                  self.c_l, self.bbox, self.dip_grid, self.dip_dir_grid, self.x, self.y, self.spacing)

        m2l_geometry.normalise_thickness(self.output_path)

        if not self.quiet:
            m2l_utils.plot_points(self.output_path+'formation_thicknesses_norm.csv',
                                  self.geol_clip, 'norm_th', 'x', 'y', True, 'numeric')

    # TODO: This needs a shorter name
    def create_fold_axial_trace_points(self, fold_decimate, fat_step, close_dip):
        # fold_decimate = 5
        folds_clip = gpd.read_file(self.fold_file)
        if(len(folds_clip) > 0):

            m2l_geometry.save_fold_axial_traces(
                self.fold_file, self.output_path, self.dtm, self.dtb, self.dtb_null, False, self.c_l, fold_decimate)

            m2l_geometry.save_fold_axial_traces_orientations(self.fold_file, self.output_path, self.tmp_path, self.dtm, self.dtb, self.dtb_null, False, self.c_l, self.proj_crs,
                                                             fold_decimate, fat_step, close_dip, self.scheme, self.bbox, self.spacing, self.dip_grid, self.dip_dir_grid)

    def postprocess(self, inputs, workflow, use_interpolations, use_fat):
        clut_path = self.clut_path
        # use_interpolations = True
        # use_fat = True

        m2l_geometry.tidy_data(self.output_path, self.tmp_path, clut_path, self.use_gcode3,
                               use_interpolations, use_fat, self.pluton_form, inputs, workflow, self.c_l)
        model_top = round(np.amax(self.dtm.read(1)), -2)

        # self.dtm.close()
        # if(workflow['cover_map']):
        #     dtb.close()

        # Calculate polarity of original bedding orientation data
        if(workflow['polarity']):
            m2l_geometry.save_orientations_with_polarity(
                self.output_path+'orientations.csv', self.output_path, self.c_l, self.tmp_path+'basal_contacts.shp', self.tmp_path+'all_sorts.csv', )

            m2l_utils.plot_points(self.output_path+'orientations_polarity.csv',
                                  self.geol_clip, 'polarity', 'X', 'Y', True, 'alpha')

        # Calculate minimum fault offset from stratigraphy and stratigraphic fault offset
        if(workflow['strat_offset']):
            fault_test = pd.read_csv(
                self.output_path+'fault_dimensions.csv', ', ')
            if(len(fault_test) > 0):

                m2l_geometry.fault_strat_offset(self.output_path, self.c_l, self.proj_crs, self.output_path+'formation_summary_thicknesses.csv', self.tmp_path +
                                                'all_sorts.csv', self.tmp_path+'faults_clip.shp', self.tmp_path+'geol_clip.shp', self.output_path+'fault_dimensions.csv')

                m2l_utils.plot_points(self.output_path+'fault_strat_offset3.csv',
                                      self.geol_clip, 'min_offset', 'X', 'Y', True, 'numeric')
                m2l_utils.plot_points(self.output_path+'fault_strat_offset3.csv',
                                      self.geol_clip, 'strat_offset', 'X', 'Y', True, 'numeric')

        # Analyse fault topologies
        Topology.parse_fault_relationships(
            self.graph_path, self.tmp_path, self.output_path)

        # TODO: Figures sometimes look a bit squashed in notebooks

    def update_projectfile(self, loopFilename):
        self.loop_projectfile = export_to_projectfile(
            loopFilename, self.output_path, self.bbox_3d, self.proj_crs)

        print("PROJECTFILE FOUND AT", self.loop_projectfile)

    def export_png(self):
        self.geology_figure.savefig("{}.png".format(self.loop_projectfile))
        print("Geology graphic exported to: " + self.tmp_path+"geology.png")
