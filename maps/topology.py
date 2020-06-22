import numpy as np
import pandas as pd

import functools
import operator


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

        f = open(geology_file_csv, "w+")
        f.write('WKT\t'+c_l['o'].replace("\n", "")+'\t'+c_l['u'].replace("\n", "")+'\t'+c_l['g'].replace("\n", "")+'\t'+c_l['min'].replace("\n", "")+'\t'+c_l['max'].replace(
            "\n", "")+'\t'+c_l['c'].replace("\n", "")+'\t'+c_l['r1'].replace("\n", "")+'\t'+c_l['r2'].replace("\n", "")+'\t'+c_l['ds'].replace("\n", "")+'\n')
        # display(sub_geol)
        print(len(sub_geol), " polygons")
        # print(sub_geol)
        for i in range(0, len(sub_geol)):
            # print('**',sub_geol.loc[i][[c_l['o']]],'++')
            f.write("\""+str(sub_geol.loc[i].geometry)+"\"\t")
            f.write("\""+str(sub_geol.loc[i][c_l['o']])+"\"\t")
            f.write("\""+str(sub_geol.loc[i][c_l['c']])+"\"\t")
            # since map2model is looking for "" not "None"
            f.write(
                "\""+str(sub_geol.loc[i][c_l['g']]).replace("None", "")+"\"\t")

            if(hint_flag == True and sub_geol.loc[i][c_l['c']] in hint_list):
                hint = code_hints.loc[sub_geol.loc[i][c_l['c']]]['hint']
            else:
                hint = 0
            min = float(sub_geol.loc[i][c_l['min']])+float(hint)
            max = float(sub_geol.loc[i][c_l['max']])+float(hint)
            # f.write("\""+str(sub_geol.loc[i][c_l['min']])+"\"\t")
            # f.write("\""+str(sub_geol.loc[i][c_l['max']])+"\"\t")
            f.write("\""+str(min)+"\"\t")
            f.write("\""+str(max)+"\"\t")
            f.write("\""+str(sub_geol.loc[i][c_l['u']])+"\"\t")
            f.write("\""+str(sub_geol.loc[i][c_l['r1']])+"\"\t")
            f.write("\""+str(sub_geol.loc[i][c_l['r2']])+"\"\t")
            f.write("\""+str(sub_geol.loc[i][c_l['ds']])+"\"\n")
        f.close()

    def save_mindep_wkt(sub_mindep, mindep_file_csv, c_l):
        f = open(mindep_file_csv, "w+")
        f.write('WKT\t'+c_l['msc']+'\t'+c_l['msn']+'\t'+c_l['mst'] +
                '\t'+c_l['mtc']+'\t'+c_l['mscm']+'\t'+c_l['mcom']+'\n')

        print(len(sub_mindep), " points")

        for i in range(0, len(sub_mindep)):
            line = "\""+str(sub_mindep.loc[i].geometry)+"\"\t\""+str(sub_mindep.loc[i][c_l['msc']])+"\"\t\"" +\
                str(sub_mindep.loc[i][c_l['msn']])+"\"\t\""+str(sub_mindep.loc[i][c_l['mst']])+"\"\t\"" +\
                str(sub_mindep.loc[i][c_l['mtc']])+"\"\t\""+str(sub_mindep.loc[i][c_l['mscm']])+"\"\t\"" +\
                str(sub_mindep.loc[i][c_l['mcom']])+"\"\n"
            f.write(functools.reduce(operator.add, (line)))

        f.close()

    def save_parfile(self, c_l, graph_path, geology_file_csv, fault_file_csv, structure_file_csv, mindep_file_csv, minx, maxx, miny, maxy, deposit_dist, commodities):
        f = open('/bin/Parfile', 'w')
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
            'Path to the output data folder                      ='+graph_path+'\n')
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
        f.close()
