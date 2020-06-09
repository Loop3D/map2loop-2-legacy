import geopandas as gpd


class ROI(object):
    def __init__(self, geology_file, fault_file, structure_file, mindep_file, bbox, polygon, c_l={}):
        self.bbox_dict = bbox
        self.bbox = tuple(list(bbox.values())[:4])
        self.polygon = polygon

        self.geology = gpd.read_file(geology_file, bbox=self.bbox)
        self.lines = gpd.read_file(fault_file, bbox=self.bbox)
        self.structures = gpd.read_file(structure_file, bbox=self.bbox)

        self.c_l = c_l

    def plot(self):
        try:
            base = self.geology.plot(column=self.c_l['c'], figsize=(
                10, 10), edgecolor='#000000', linewidth=0.2)
            self.structures.plot(ax=base, color='none', edgecolor='black')
            self.lines.plot(ax=base, cmap='rainbow',
                            column=self.c_l['f'], figsize=(10, 10), linewidth=0.4)
            # TODO: - polygon roi not plotting correctly
            #       - make a plot function in project that plots the above stuff in here
            # polygon.plot(ax=base, color='none', edgecolor='black')

        except Exception as e:
            print(e)
