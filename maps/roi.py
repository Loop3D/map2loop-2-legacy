import geopandas as gpd


class ROI(object):
    def __init__(self, geology_file, fault_file, structure_file, mindep_file, bbox, polygon, c_l={}):
        self.bbox = bbox
        self.polygon = polygon

        try:
            geology = gpd.read_file(geology_file, bbox=bbox)
            lines = gpd.read_file(fault_file, bbox=bbox)
            structures = gpd.read_file(structure_file, bbox=bbox)

            base = geology.plot(column=c_l['c'], figsize=(
                10, 10), edgecolor='#000000', linewidth=0.2)
            structures.plot(ax=base, color='none', edgecolor='black')
            lines.plot(ax=base, cmap='rainbow',
                       column=c_l['f'], figsize=(10, 10), linewidth=0.4)
            # TODO: - polygon roi not plotting correctly
            #       - make a plot function in project that plots the above stuff in here
            # polygon.plot(ax=base, color='none', edgecolor='black')

        except Exception as e:
            print(e)
