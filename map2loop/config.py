import os
import time

import numpy as np
import pandas as pd
import geopandas as gpd
from map2loop.topology import Topology
from map2loop import m2l_utils
import map2model

import networkx as nx
import matplotlib.pyplot as plt
import rasterio
import shapely


class Config(object):
    def __init__(self, geology_file, fault_file, structure_file, mindep_file, bbox_3d, polygon, step_out, src_crs, dst_crs, c_l={}):
        self.bbox_3d = bbox_3d
        self.bbox = tuple([bbox_3d["minx"], bbox_3d["miny"],
                           bbox_3d["maxx"], bbox_3d["maxy"]])
        self.polygon = polygon
        self.step_out = step_out
        self.c_l = c_l

        self.src_crs = src_crs
        self.dst_crs = dst_crs

        # Create file structure for model instance
        self.geology_file = geology_file
        self.structure_file = structure_file
        self.fault_file = fault_file
        self.mindep_file = mindep_file

        self.project_path = 'model-test'

        self.graph_path = self.project_path+'/graph'
        self.tmp_path = self.project_path+'/tmp'
        self.data_path = self.project_path+'/data'
        self.dtm_path = self.project_path+'/dtm'
        self.output_path = self.project_path+'/output'
        self.vtk_path = self.project_path+'/vtk'

        self.fault_file_csv = self.tmp_path + "/faults.csv"
        self.structure_file_csv = self.tmp_path + "/structure.csv"
        self.geology_file_csv = self.tmp_path + "/geology.csv"
        self.mindep_file_csv = self.tmp_path + "/mindep.csv"

        self.strat_graph_file = self.graph_path + "/graph_strat_NONE.gml"
        self.dtm_file = self.dtm_path+'/dtm.tif'
        self.dtm_reproj_file = self.dtm_path+'/dtm_rp.tif'

        if(not os.path.isdir(self.project_path)):
            os.mkdir(self.project_path)
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

    def preprocess(self, command=""):
        geology = gpd.read_file(self.geology_file, bbox=self.bbox)
        faults = gpd.read_file(self.fault_file, bbox=self.bbox)
        structures = gpd.read_file(self.structure_file, bbox=self.bbox)
        mindeps = gpd.read_file(self.mindep_file, bbox=self.bbox)

        self.geology = geology
        self.faults = faults
        self.structures = structures
        self.mindeps = mindeps

        if command == "plot":
            try:
                base = geology.plot(column=self.c_l['c'], figsize=(
                    10, 10), edgecolor='#000000', linewidth=0.2, legend=True)
                leg = base.get_legend()
                leg.set_bbox_to_anchor((1.04, 1))
                self.geology_figure = base.get_figure()

                structures.plot(ax=base, color='none', edgecolor='black')

                faults.plot(ax=base, cmap='rainbow',
                            column=self.c_l['f'], figsize=(10, 10), linewidth=0.4)
                structures[['geometry', self.c_l['gi'],
                            self.c_l['d'], self.c_l['dd']]].plot(ax=base)

                fig = self.polygon.plot(ax=base, color='none', edgecolor='black').set_title(
                    "Input {}".format(self.bbox)).get_figure()
                fig.savefig(self.tmp_path+"/input-data.png")

                print("Input graphic saved to: " +
                      self.tmp_path + "/input-fig.png")

                self.export_png()
                plt.show()

                return
            except Exception as e:
                print(e)

    def export_png(self):
        self.geology_figure.savefig(self.tmp_path+"/geology.png")
        print("Geology graphic exported to: " + self.tmp_path+"/geology.png")

    def export_csv(self):
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
        sub_faults = self.faults[['geometry', self.c_l['o'], self.c_l['f']]]
        Topology.save_faults_wkt(sub_faults, self.fault_file_csv, self.c_l)

    def update_parfile(self):
        Topology.save_parfile(self, self.c_l, self.output_path, self.geology_file_csv, self.fault_file_csv, self.structure_file_csv,
                              self.mindep_file_csv, self.bbox[0], self.bbox[1], self.bbox[2], self.bbox[3], 500.0, 'Fe,Cu,Au,NONE')

    def runMap2Model(self):
        print(map2model.run(self.graph_path, self.geology_file_csv,
                            self.fault_file_csv, self.mindep_file_csv,
                            self.bbox_3d,
                            self.c_l,
                            "Fe,Cu,Au,NONE"))

        print("Resolving ambiguities using ASUD...", end='\toutput_dir:')
        aus = True
        if aus:
            Topology.use_asud(self.strat_graph_file,  self.graph_path)
            self.strat_graph_file = self.graph_path+'/ASUD_strat.gml'
        print("Done.")

        print("Generating topology graph display and unit groups...")
        self.G = nx.read_gml(self.strat_graph_file, label='id')
        selected_nodes = [n for n, v in self.G.nodes(data=True) if n >= 0]

        plt.figure(1, figsize=(19, 12))
        nx.draw_networkx(self.G, pos=nx.kamada_kawai_layout(
            self.G), arrows=True, nodelist=selected_nodes)
        nlist = list(self.G.nodes.data('LabelGraphics'))
        nlist.sort()
        for node in nlist:
            if node[0] >= 0:
                elem = str(node[1]).replace(
                    "{'text':", "").replace(", 'fontSize': 14}", "")
                # second=elem.split(":").replace("'","")
                print(node[0], " ", elem)

        plt.axis("off")
        plt.savefig(self.tmp_path+"/topology-fig.png")
        print("Topology figure saved to", self.tmp_path+"/topology-fig.png")

        # Save groups of stratigraphic units
        groups, self.glabels, G = Topology.get_series(
            self.strat_graph_file, 'id')
        Topology.save_units(G, self.tmp_path, self.glabels,
                            Australia=True, asud_strat_file="https://gist.githubusercontent.com/yohanderose/3b257dc768fafe5aaf70e64ae55e4c42/raw/8598c7563c1eea5c0cd1080f2c418dc975cc5433/ASUD.csv")

        print("Done")

    def load_dtm(self):
        polygon_ll = self.polygon.to_crs(self.src_crs)

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
            self.dtm_file, self.dtm_reproj_file, self.src_crs, self.dst_crs)
        self.dtm = rasterio.open(self.dtm_reproj_file)
        fig = plt.imshow(self.dtm.read(1), cmap='terrain',
                         vmin=0, vmax=1000).get_figure()
        plt.tight_layout(5)
        plt.show()

    def join_features(self):
        # Geology
        geol_clip = m2l_utils.explode(self.geology)
        geol_clip.crs = self.dst_crs
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)

        # Structures
        list1 = ['geometry', self.c_l['d'],
                 self.c_l['dd'], self.c_l['sf'], self.c_l['bo']]
        list2 = list(set(list1))
        sub_pts = self.structures[list2]
        structure_code = gpd.sjoin(sub_pts, geol_clip, how="left", op="within")

        minx, miny, maxx, maxy = self.bbox
        y_point_list = [miny, miny, maxy, maxy, miny]
        x_point_list = [minx, maxx, maxx, minx, minx]

        bbox_geom = shapely.geometry.Polygon(zip(x_point_list, y_point_list))

        # TODO: 'polygo' is never used
        polygo = gpd.GeoDataFrame(
            index=[0], crs=self.dst_crs, geometry=[bbox_geom])
        is_bed = structure_code[self.c_l['sf']].str.contains(
            self.c_l['bedding'], regex=False)

        structure_clip = structure_code[is_bed]
        structure_clip.crs = self.dst_crs

        if(self.c_l['otype'] == 'strike'):
            structure_clip['azimuth2'] = structure_clip.apply(
                lambda row: row[self.c_l['dd']]+90.0, axis=1)
            self.c_l['dd'] = 'azimuth2'
            self.c_l['otype'] = 'dip direction'

        structure_clip = structure_clip[~structure_clip[self.c_l['o']].isnull(
        )]
        structure_clip.to_file(self.tmp_path+'/structure_clip.shp')

        # Save geology clips
        Topology.save_group(self.G, self.tmp_path,
                            self.glabels, geol_clip, self.c_l)
