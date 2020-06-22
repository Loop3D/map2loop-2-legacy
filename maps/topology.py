import numpy as np
import pandas as pd


class Topology(object):
    def __init__(self):
        pass

    def save_geol_wkt(sub_geol, geology_file_csv, c_l, hint_flag):
        # hint_flag=False
        if(hint_flag == True):
            print("Using ENS age hints")
            code_hint_file = '../test_data3/data/code_age_hint.csv'  # input code_hint file
            code_hints = pd.read_csv(code_hint_file, sep=',')
            code_hints.drop_duplicates(inplace=True)
            code_hints.set_index('code',  inplace=True)

            hint_list = []
            for indx, row in code_hints.iterrows():
                if(not indx in hint_list):
                    hint_list.append(indx)

        # TODO: Understand and apply hints
        sub_geol.columns.values[0] = "WKT"
        # map2model needs '' as unassigned cells, make nones NaNs
        sub_geol.replace(to_replace=[None], value=np.nan, inplace=True)
        sub_geol.fillna(value=np.nan, inplace=True)
        # then make NaNs empties
        sub_geol.replace(np.nan, '', regex=True)
        sub_geol.replace('None', '', inplace=True)
        sub_geol.to_csv(geology_file_csv, index=False)

    def save_mindep_wkt(sub_mindep, mindep_file_csv, c_l):
        sub_mindep.columns.values[0] = "WKT"
        sub_mindep.to_csv(mindep_file_csv, index=False)

    def save_structure_wkt(sub_pts, structure_file_csv, c_l):
        sub_pts.columns.values[0] = "WKT"
        sub_pts.to_csv(structure_file_csv, index=False)

    def save_faults_wkt(sub_lines, fault_file_csv, c_l):
        sub_lines.columns.values[0] = "WKT"
        faults_only = sub_lines[sub_lines["FEATURE"].str.contains(
            'fault', case=False)]
        faults_only.to_csv(fault_file_csv, index=False)

    def save_parfile(self, c_l, data_path, geology_file_csv, fault_file_csv, structure_file_csv, mindep_file_csv, minx, maxx, miny, maxy, deposit_dist, commodities):
        with open('Parfile', 'w') as f:
            f.write(
                '--- COLUMN NAMES IN CSV DATA FILES: -------------------------------------------------------------\n')
            f.write('OBJECT COORDINATES              =WKT\n')
            f.write('FAULT: ID                       ='+c_l['o']+'\n')
            f.write('FAULT: FEATURE                  ='+c_l['f']+'\n')
            f.write('POLYGON: ID                     ='+c_l['o']+'\n')
            f.write('POLYGON: LEVEL1 NAME            ='+c_l['u']+'\n')
            f.write('POLYGON: LEVEL2 NAME            ='+c_l['g']+'\n')
            f.write('POLYGON: MIN AGE                ='+c_l['min']+'\n')
            f.write('POLYGON: MAX AGE                ='+c_l['max']+'\n')
            f.write('POLYGON: CODE                   ='+c_l['c']+'\n')
            f.write('POLYGON: DESCRIPTION            ='+c_l['ds']+'\n')
            f.write('POLYGON: ROCKTYPE1              ='+c_l['r1']+'\n')
            f.write('POLYGON: ROCKTYPE2              ='+c_l['r2']+'\n')
            f.write('DEPOSIT: SITE CODE              ='+c_l['msc']+'\n')
            f.write('DEPOSIT: SITE TYPE              ='+c_l['mst']+'\n')
            f.write('DEPOSIT: SITE COMMODITY         ='+c_l['mscm']+'\n')

            f.write(
                '--- SOME CONSTANTS: ----------------------------------------------------------------------------\n')
            f.write('FAULT AXIAL FEATURE NAME        ='+c_l['fold']+'\n')
            f.write('SILL UNIT DESCRIPTION CONTAINS  ='+c_l['sill']+'\n')
            f.write('IGNEOUS ROCKTYPE CONTAINS                           =' +
                    c_l['intrusive']+'\n')
            f.write('VOLCANIC ROCKTYPE CONTAINS                          =' +
                    c_l['volcanic']+'\n')
            f.write(
                'IGNORE DEPOSITS WITH SITE TYPE                      =Infrastructure\n')
            f.write('Intersect Contact With Fault: angle epsilon (deg)   =1.0\n')
            f.write('Intersect Contact With Fault: distance epsilon (m)  =15.0\n')
            f.write('Distance buffer (fault stops on another fault) (m)  =20.0\n')
            f.write('Distance buffer (point on contact) (m)              =' +
                    str(deposit_dist)+'\n')
            f.write(
                '------------------------------------------------------------------------------------------------\n')
            f.write(
                'Path to the output data folder                      ='+data_path+'\n')
            f.write('Path to geology data file                           =' +
                    geology_file_csv+'\n')
            f.write('Path to faults data file                            =' +
                    fault_file_csv+'\n')
            f.write('Path to mineral deposits data file                  =' +
                    mindep_file_csv+'\n')
            f.write(
                '------------------------------------------------------------------------------------------------\n')
            f.write('Clipping window X1 Y1 X2 Y2 (zeros for infinite)    =' +
                    str(minx)+' '+str(miny)+' '+str(maxx)+' '+str(maxy)+'\n')
            f.write('Min length fraction for strat/fault graphs          =0.0\n')
            f.write(
                'Graph edge width categories (three doubles)         =2000. 20000. 200000.\n')
            f.write('Graph edge direction (0-min age, 1-max age, 2-avg)  =2\n')
            f.write(
                'Deposit names for adding info on the graph          ='+commodities+'\n')
            f.write('Partial graph polygon ID                            =32\n')
            f.write('Partial graph depth                                 =4\n')
            f.write(
                'Map subregion size dx, dy [m] (zeros for full map)  =0. 0.\n')
            f.write(
                '------------------------------------------------------------------------------------------------\n')
