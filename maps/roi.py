import geopandas as gpd


class ROI(object):
    def __init__(self, geology_file, fault_file, structure_file, mindep_file, bbox, c_l={}):

        try:
            geology = gpd.read_file(geology_file, bbox=bbox)
            lines = gpd.read_file(fault_file, bbox=bbox)
            structures = gpd.read_file(structure_file, bbox=bbox)

            base = geology.plot(column=c_l['c'], figsize=(
                10, 10), edgecolor='#000000', linewidth=0.2)

        except Exception as e:
            print(e)
