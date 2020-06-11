import geopandas as gpd
import os


class Model(object):
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

        self.test_data_path = './model-test/'

        self.graph_path = self.test_data_path+'graph/'
        self.tmp_path = self.test_data_path+'tmp/'
        self.data_path = self.test_data_path+'data/'
        self.dtm_path = self.test_data_path+'dtm/'
        self.output_path = self.test_data_path+'output/'
        self.vtk_path = self.test_data_path+'vtk/'

        self.fault_file_csv = self.tmp_path + "faults.csv"
        self.fault_file_csv = self.tmp_path + "structure.csv"
        self.geology_file_csv = self.tmp_path + "geology.csv"
        self.mindep_file_csv = self.tmp_path + "mindep.csv"

        self.strat_graph_file = self.graph_path + "graph_strat_NONE.gml"
        self.dtm_file = self.dtm_path+'dtm.tif'
        self.dtm_reproj_file = self.dtm_path+'dtm_rp.tif'

        if(not os.path.isdir(self.test_data_path)):
            os.mkdir(self.test_data_path)
        if(not os.path.isdir(self.tmp_path)):
            os.mkdir(self.tmp_path)
        if(not os.path.isdir(self.output_path)):
            os.mkdir(self.output_path)
        if(not os.path.isdir(self.dtm_path)):
            os.mkdir(self.dtm_path)
        if(not os.path.isdir(self.vtk_path)):
            os.mkdir(self.vtk_path)
        if(not os.path.isdir(self.graph_path)):
            os.mkdir(self.graph_path)

    def plot(self):
        geology = gpd.read_file(self.geology_file, bbox=self.bbox)
        lines = gpd.read_file(self.fault_file, bbox=self.bbox)
        structures = gpd.read_file(self.structure_file, bbox=self.bbox)

        try:
            base = geology.plot(column=self.c_l['c'], figsize=(
                10, 10), edgecolor='#000000', linewidth=0.2)
            structures.plot(ax=base, color='none', edgecolor='black')
            lines.plot(ax=base, cmap='rainbow',
                       column=self.c_l['f'], figsize=(10, 10), linewidth=0.4)
            self.polygon.plot(ax=base, color='none', edgecolor='black')

        except Exception as e:
            print(e)
