import geopandas as gpd
import os
from maps.topology import Topology


class Config(object):
    def __init__(self, geology_file, fault_file, structure_file, mindep_file, bbox_3d, polygon, step_out, c_l={}):
        self.bbox_3d = bbox_3d
        self.bbox = tuple(list(bbox_3d.values())[:4])
        self.polygon = polygon
        self.step_out = step_out
        self.c_l = c_l

        # Create file structure for model instance
        self.geology_file = geology_file
        self.structure_file = structure_file
        self.fault_file = fault_file
        self.mindep_file = mindep_file

        self.project_path = './model-test'

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

    def preprocess(self, command="plot"):
        geology = gpd.read_file(self.geology_file, bbox=self.bbox)
        lines = gpd.read_file(self.fault_file, bbox=self.bbox)
        structures = gpd.read_file(self.structure_file, bbox=self.bbox)
        mindep = gpd.read_file(self.mindep_file, bbox=self.bbox)

        if command == "plot":
            try:
                base = geology.plot(column=self.c_l['c'], figsize=(
                    10, 10), edgecolor='#000000', linewidth=0.2)
                structures.plot(ax=base, color='none', edgecolor='black')
                lines.plot(ax=base, cmap='rainbow',
                           column=self.c_l['f'], figsize=(10, 10), linewidth=0.4)
                structures[['geometry', self.c_l['gi'],
                            self.c_l['d'], self.c_l['dd']]].plot(ax=base)
                self.polygon.plot(ax=base, color='none', edgecolor='black')

            except Exception as e:
                print(e)

        elif command == "export_csv":
            # Save geology polygons
            hint_flag = False  # use GSWA strat database to provide relative age hints
            sub_geol = geology[['geometry', self.c_l['o'], self.c_l['c'], self.c_l['g'],
                                self.c_l['u'], self.c_l['min'], self.c_l['max'], self.c_l['ds'], self.c_l['r1'], self.c_l['r2']]]
            Topology.save_geol_wkt(
                sub_geol, self.geology_file_csv, self.c_l, hint_flag)

            # Save mineral deposits
            sub_mindep = mindep[['geometry', self.c_l['msc'], self.c_l['msn'],
                                 self.c_l['mst'], self.c_l['mtc'], self.c_l['mscm'], self.c_l['mcom']]]
            Topology.save_mindep_wkt(
                sub_mindep, self.mindep_file_csv, self.c_l)

            # Save orientation data
            sub_pts = structures[[
                'geometry', self.c_l['gi'], self.c_l['d'], self.c_l['dd']]]
            Topology.save_structure_wkt(
                sub_pts, self.structure_file_csv, self.c_l)

            # Save faults
            sub_lines = lines[['geometry', self.c_l['o'], self.c_l['f']]]
            Topology.save_faults_wkt(sub_lines, self.fault_file_csv, self.c_l)

            # Return those sub frames for debugging
            # return sub_geol, sub_lines, sub_pts, sub_mindep

    def update_parfile(self):
        Topology.save_parfile(self, self.c_l, self.output_path, self.geology_file_csv, self.fault_file_csv, self.structure_file_csv,
                              self.mindep_file_csv, self.bbox[0], self.bbox[1], self.bbox[2], self.bbox[3], 500.0, 'Fe,Cu,Au,NONE')
