import os
from map2loop.topology import Topology
import pandas as pd
import geopandas as gpd
import numpy as np
import networkx as nx
from geopandas import GeoDataFrame
from pandas import DataFrame
from shapely.geometry import  Point
from shapely.ops import snap
from map2loop import m2l_utils
from math import acos, degrees
from . import m2l_utils
from . import project as Project
from . import config as Config
from .m2l_utils import display, print

close_f=1000
close_b=1000
close_i=1000
fault_clip_buffer=.01
snap_buffer=.5


class Map2Graph(object):

    def __init__(self,Gloop=None,c_l=None):
        self.Gloop=nx.DiGraph()
        pass

    def clean_inputs(geology_file,fault_file,mindep_file):
        mindep=gpd.read_file(mindep_file)
        fault=gpd.read_file(fault_file)
        geology=gpd.read_file(geology_file)

        fault_clean = fault.dropna(subset=['geometry'])
        geology_clean_nona = geology.dropna(subset=['geometry'])

        fault_zone = fault_clean.buffer(fault_clip_buffer)
        all_fz = fault_zone.unary_union

        geology_clean=geology_clean_nona.difference(all_fz)
        geology_clean=gpd.GeoDataFrame(geology, crs=geology.crs, geometry=geology_clean.geometry) 
        geology_clean.reset_index(inplace=True)
        geology_clean['idx']=geology_clean.index        
        geology_clean_nona.reset_index(inplace=True)
        geology_clean_nona['idx']=geology_clean_nona.index

        return(mindep,fault_clean,geology_clean,geology_clean_nona)
    
    def strat_graph(geology_clean,output_path,c_l):

        all_contacts=[]
        for ind,g in geology_clean.iterrows():
            for ind2,g2 in geology_clean.iterrows():
                if(not ind >= ind2):
                    if(g.geometry.intersects(g2.geometry)):
                        g1_snapped = snap(g.geometry, g2.geometry,snap_buffer)
                        all_contacts.append([ind,ind2,g1_snapped.buffer(0).intersection(g2.geometry.buffer(0))])

        groups=geology_clean.drop_duplicates(subset=[c_l['g']])
        #groups.reset_index(inplace=True)
        groups.index= pd.RangeIndex(start=0, stop=len(groups), step=1)
        groups['order'] = groups.index
        groups.set_index([c_l['g']],inplace=True)
            
        strats=geology_clean.drop_duplicates(subset=[c_l['c']])
        strats[c_l['c']]=strats[c_l['c']].str.replace(' ','_').replace('-','_')
        strats.set_index([c_l['c']],inplace=True)

        Gloop=nx.DiGraph()

        i=[]
        index=0
        for ind,g in groups.iterrows():
            Gloop.add_node(g.name.replace(" ","_").replace("-","_")+"_gp",id=index, isGroup=1,LabelGraphics ='[ text "'+g.name.replace(" ","_").replace("-","_")+"_gp"+'" anchor "n" fontStyle "bold" fontSize 14 ]')
            Gloop.nodes[g.name.replace(" ","_").replace("-","_")+"_gp"]['ntype']='group'
            index=index+1
            i.append(index)
        groups['gid']=i

        for ind,s in strats.iterrows():
            Gloop.add_node(s.name,id=index, gid=str(groups.loc[s[c_l['g']]]['gid']),LabelGraphics ='[ text "'+s.name.replace(" ","_").replace("-","_")+'" anchor "n" fontStyle "bold" fontSize 14 ]')
            Gloop.nodes[s.name]['ntype']='formation'
            
            index=index+1


        for c in all_contacts:
            ave_age0=int(geology_clean.iloc[c[0]][c_l['min']])+((int(geology_clean.iloc[c[0]][c_l['max']])-int(geology_clean.iloc[c[0]][c_l['min']]))/2)
            ave_age1=int(geology_clean.iloc[c[1]][c_l['min']])+((int(geology_clean.iloc[c[1]][c_l['max']])-int(geology_clean.iloc[c[1]][c_l['min']]))/2)

            if(ave_age0>ave_age1):
                Gloop.add_edge(geology_clean.iloc[c[1]][c_l['c']].replace(" ","_").replace("-","_"),geology_clean.iloc[c[0]][c_l['c']].replace(" ","_").replace("-","_"))
                Gloop[geology_clean.iloc[c[1]][c_l['c']].replace(" ","_").replace("-","_")][geology_clean.iloc[c[0]][c_l['c']].replace(" ","_").replace("-","_")]['etype']='formation_formation'
            else:
                Gloop.add_edge(geology_clean.iloc[c[0]][c_l['c']].replace(" ","_").replace("-","_"),geology_clean.iloc[c[1]][c_l['c']].replace(" ","_").replace("-","_"))
                Gloop[geology_clean.iloc[c[0]][c_l['c']].replace(" ","_").replace("-","_")][geology_clean.iloc[c[1]][c_l['c']].replace(" ","_").replace("-","_")]['etype']='formation_formation'
                
        nx.write_gml(Gloop, os.path.join(output_path,'pre_loop.gml'))
        return(groups,strats,all_contacts,Gloop)

    def fix_Loop_graph(output_path,prefix):


        with open(os.path.join(output_path,prefix+'.gml')) as graph:
            lines = graph.readlines()
            graph.close()
            while(not graph.closed):
                pass
            new_graph=open(os.path.join(output_path,prefix+'.gml'),"w")
            for l in lines:
                l=l.replace('&#34;','"').replace('"[','[').replace(']"',']')
                if('gid "' in l):
                    tmp=l.split('"')
                    new_graph.write('    gid '+str(int(tmp[1])-1)+'\n')                
                else:
                    new_graph.write(l)
                if('ntype "formation"' in l):
                    new_graph.write('    graphics [ type "rectangle" fill "#FFFFFF" ]\n')                

            new_graph.close()

        
    
    def basal_contacts(Gloop,output_path,geology_clean,all_contacts,groups,c_l,formation_weight):
        Gasud=nx.read_gml(os.path.join(output_path,'ASUD_strat.gml'))

        cycles = nx.simple_cycles(Gasud)

        for cy in cycles:
            if(len(cy)==1):
                start=cy[0]
                end=cy[0]
            else:
                start=cy[0]
                end=cy[1]
                
            warning_msg = 'map2loop warning: Stratigraphic relationship: ' + \
                str(Gasud.nodes[start]['LabelGraphics']['text'])+' overlies '\
                +str(Gasud.nodes[end]['LabelGraphics']['text'])+' removed to prevent cycle'
            print(warning_msg)
            if Gasud.has_edge(start,end):
                Gasud.remove_edge(start,end)    

        strat_order=list(nx.topological_sort(Gasud))
        strat_max=len(strat_order)-len(groups)
        strat_order=strat_order[:strat_max]

        index=1
        mini_strat={}
        for s in strat_order:
            Gasud.add_node(s,order=index)
            Gloop.add_node(Gasud.nodes[s]['LabelGraphics']['text'],order=index,formation_weight=formation_weight)
            print('strat,order',Gasud.nodes[s]['LabelGraphics']['text'],index)
            
            mini_strat[index]={"strat":Gasud.nodes[s]['LabelGraphics']['text'],"order":index}
            index=index+1   

        mini_strat_df=pd.DataFrame.from_dict(mini_strat,orient="index")
        mini_strat_df.set_index("strat",inplace=True)

        igneous_contacts={}
        not_igneous_contacts={}
        i=0
        for c in all_contacts:
            if(c[2].geom_type=='MultiLineString' or c[2].geom_type=='LineString' ): 
                if( c_l['intrusive'] in geology_clean.iloc[c[0]][c_l['r1']] ):
                    igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[0]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[1]][c_l['c']], "geometry": c[2]}
                elif(c_l['intrusive'] in geology_clean.iloc[c[1]][c_l['r1']] ):
                    igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[1]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[0]][c_l['c']], "geometry": c[2]}
                else:
                    if(c[0]][c_l['c'].replace(" ","_").replace("-","_") in mini_strat_df.index):
                        if(mini_strat_df.loc[geology_clean.iloc[c[0]][c_l['c']].replace(" ","_").replace("-","_")]['order']>
                        mini_strat_df.loc[geology_clean.iloc[c[1]][c_l['c']].replace(" ","_").replace("-","_")]['order']):
                            not_igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[1]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[0]][c_l['c']], "geometry": c[2]}
                        else:
                            not_igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[0]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[1]][c_l['c']], "geometry": c[2]}
                        
                i=i+1
            elif(c[2].geom_type=='GeometryCollection' ):
                for geom in c[2]:
                    if(geom.geom_type=='MultiLineString' or geom.geom_type=='LineString' ):
                        if( c_l['intrusive'] in geology_clean.iloc[c[0]][c_l['r1']] ):
                            igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[0]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[1]][c_l['c']], "geometry": geom}
                        elif(c_l['intrusive'] in geology_clean.iloc[c[1]][c_l['r1']] ):
                            igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[1]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[0]][c_l['c']], "geometry": geom}
                        else:
                            if(mini_strat_df.loc[geology_clean.iloc[c[0]][c_l['c']].replace(" ","_").replace("-","_")]['order']>
                            mini_strat_df.loc[geology_clean.iloc[c[1]][c_l['c']].replace(" ","_").replace("-","_")]['order']):
                                not_igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[1]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[0]][c_l['c']], "geometry": geom}
                            else:
                                not_igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[0]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[1]][c_l['c']], "geometry": geom}
                        i=i+1
            
        df = DataFrame.from_dict(igneous_contacts, "index")
        if(len(df)>0):
            i_contacts_gdf = GeoDataFrame(df, crs=geology_clean.crs, geometry='geometry')
        else:
            i_contacts_gdf=DataFrame.from_dict({}, "index")

        df = DataFrame.from_dict(not_igneous_contacts, "index")
        if(len(df)>0):
            b_contacts_gdf = GeoDataFrame(df, crs=geology_clean.crs, geometry='geometry')
        else:
            b_contacts_gdf=DataFrame.from_dict({}, "index")

        return(Gloop,i_contacts_gdf,b_contacts_gdf)

    def fault_intersections(Gloop,output_path,fault_clean,c_l,fault_weight,fault_fault_weight):

        Gfault=nx.DiGraph()
        fault_intersections=[]

        for ind,f in fault_clean.iterrows():
            Gloop.add_node('Fault_'+f[c_l['o']])
            Gloop.nodes['Fault_'+f[c_l['o']]]['ntype']='fault'
            Gloop.nodes['Fault_'+f[c_l['o']]]['weight']=fault_weight
            Gfault.add_node('Fault_'+f[c_l['o']])
            Gfault.nodes['Fault_'+f[c_l['o']]]['ntype']='fault'
            for ind2,f2 in fault_clean.iterrows():
                if(not fault_clean.iloc[ind][c_l['o']]>=fault_clean.iloc[ind2][c_l['o']]):
                    if(f.geometry.intersects(f2.geometry)):
                        fault_intersections.append([fault_clean.iloc[ind][c_l['o']],fault_clean.iloc[ind2][c_l['o']],f.geometry.intersection(f2.geometry)])

        fault_name=fault_clean.set_index(c_l['o'])
        eps=.01
        min_ang=30
        for ind,f in fault_name.iterrows():
            for e in fault_intersections:
                if(e[0]==f.name):
                    other=e[1]
                else:
                    other=e[0]
                if(e[2].distance(Point(f.geometry.coords[0]))<eps):
                    ind2=0
                    for pt in fault_name.loc[other].geometry.coords[1:-1]:
                        if(Point(pt).distance(e[2])<eps):
                            l1,m1=m2l_utils.pts2dircos(f.geometry.coords[0][0],f.geometry.coords[0][1],f.geometry.coords[1][0],f.geometry.coords[1][1])
                            l2,m2=m2l_utils.pts2dircos(fault_name.loc[other].geometry.coords[ind2-1][0],fault_name.loc[other].geometry.coords[ind2-1][1],
                                        fault_name.loc[other].geometry.coords[ind2+1][0],fault_name.loc[other].geometry.coords[ind2+1][1])
                            ang=degrees(acos(l1*l2+m1*m2))
                            if(ang>90):
                                ang=180-ang
                            if(ang<min_ang):
                                ang_code='T'
                            else:
                                ang_code='X'

                            Gloop.add_edge('Fault_'+f.name,'Fault_'+other)
                            Gloop['Fault_'+f.name]['Fault_'+other]['etype']='fault_fault'
                            Gloop['Fault_'+f.name]['Fault_'+other]['fault1']='Fault_'+f.name
                            Gloop['Fault_'+f.name]['Fault_'+other]['fault2']='Fault_'+other
                            Gloop['Fault_'+f.name]['Fault_'+other]['angle']=int(ang)
                            Gloop['Fault_'+f.name]['Fault_'+other]['topol']=ang_code
                            Gloop['Fault_'+f.name]['Fault_'+other]['weight']=fault_fault_weight
                            
                            """Gfault.add_edge('Fault_'+f.name,'Fault_'+other)
                            Gfault['Fault_'+f.name]['Fault_'+other]['etype']='fault_fault'
                            Gfault['Fault_'+f.name]['Fault_'+other]['angle']=int(ang)
                            Gfault['Fault_'+f.name]['Fault_'+other]['topol']=ang_code"""
                        ind2=ind2+1
                
                elif(e[2].distance(Point(f.geometry.coords[-1]))<eps):
                    ind2=0
                    for pt in fault_name.loc[other].geometry.coords[1:-1]:
                        if(Point(pt).distance(e[2])<eps):
                            l1,m1=m2l_utils.pts2dircos(f.geometry.coords[-1][0],f.geometry.coords[-1][1],f.geometry.coords[-2][0],f.geometry.coords[-2][1])
                            l2,m2=m2l_utils.pts2dircos(fault_name.loc[other].geometry.coords[ind2-1][0],fault_name.loc[other].geometry.coords[ind2-1][1],
                                        fault_name.loc[other].geometry.coords[ind2+1][0],fault_name.loc[other].geometry.coords[ind2+1][1])
                            ang=degrees(acos(l1*l2+m1*m2))
                            if(ang>90):
                                ang=180-ang
                            if(ang<min_ang):
                                ang_code='T'
                            else:
                                ang_code='X'
                            Gloop.add_edge('Fault_'+f.name,'Fault_'+other)
                            Gloop['Fault_'+f.name]['Fault_'+other]['etype']='fault_fault'
                            Gloop['Fault_'+f.name]['Fault_'+other]['fault1']='Fault_'+f.name
                            Gloop['Fault_'+f.name]['Fault_'+other]['fault2']='Fault_'+other
                            Gloop['Fault_'+f.name]['Fault_'+other]['angle']=int(ang)
                            Gloop['Fault_'+f.name]['Fault_'+other]['topol']=ang_code
                            Gloop['Fault_'+f.name]['Fault_'+other]['weight']=fault_fault_weight
                            
                            """Gfault.add_edge('Fault_'+f.name,'Fault_'+other)
                            Gfault['Fault_'+f.name]['Fault_'+other]['etype']='fault_fault'
                            Gfault['Fault_'+f.name]['Fault_'+other]['angle']=int(ang)
                            Gfault['Fault_'+f.name]['Fault_'+other]['topol']=ang_code"""

                        ind2=ind2+1

            nx.write_gml(Gfault, os.path.join(output_path,'pre_loop_fault_network.gml'))

            return(Gloop)
    
    def fault_formation_intersections(Gloop,output_path,fault_clean,groups,geology_clean_nona,c_l,fault_formation_weight):
        fault_group_array=np.zeros((len(groups),len(fault_clean)), dtype=np.int8)
        fault_names='Fault_'+fault_clean[c_l['o']]
        group_names=list(groups.index)
        group_names=[w.replace(" ","_").replace("-","_") for w in group_names]

        i1=0
        for ind,f in fault_clean.iterrows():
            for ind2,g in geology_clean_nona.iterrows():
                if(f.geometry.intersects(g.geometry)):
                    g1_snapped = snap(g.geometry, f.geometry,snap_buffer)
                    x=g1_snapped.intersection(f.geometry)
                    if(x.geom_type!='Point' and x.geom_type!='MultiPoint'):
                        Gloop.add_edge('Fault_'+f[c_l['o']],g[c_l['g']].replace(" ","_").replace("-","_")+'_gp')
                        Gloop['Fault_'+f[c_l['o']]][g[c_l['g']].replace(" ","_").replace("-","_")+'_gp']['etype']='fault_group'

                        Gloop.add_edge('Fault_'+f[c_l['o']],g[c_l['c']].replace(" ","_").replace("-","_"))
                        Gloop['Fault_'+f[c_l['o']]][g[c_l['c']].replace(" ","_").replace("-","_")]['etype']='fault_formation'
                        fault_group_array[groups.loc[g[c_l['g'].replace(" ","_").replace("-","_")]]['order'],i1]=1
                        Gloop['Fault_'+f[c_l['o']]][g[c_l['c']].replace(" ","_").replace("-","_")]['weight']=fault_formation_weight
                        Gloop['Fault_'+f[c_l['o']]][g[c_l['c']].replace(" ","_").replace("-","_")]['fault']='Fault_'+f[c_l['o']]
                        Gloop['Fault_'+f[c_l['o']]][g[c_l['c']].replace(" ","_").replace("-","_")]['formation']=str(g[c_l['c']].replace(" ","_").replace("-","_"))+'_'+str(g['idx'])
            i1=i1+1
        
        fault_group_df=pd.DataFrame(fault_group_array,index=group_names,columns=fault_names)
        fault_group_df.to_csv(os.path.join(output_path,'group-fault-relationships.csv'),index_label='group') 

        return(Gloop)

    def remove_group_info(Gloop):
        for n in Gloop.nodes():
            if('isGroup' in Gloop.nodes[n]):
                del Gloop.nodes[n]['isGroup']
            elif('gid' in Gloop.nodes[n]):
                del Gloop.nodes[n]['gid']  
        return(Gloop)

    def  group_formation_intersections(Gloop,strats,c_l): 
        for ind,s in strats.iterrows():
            Gloop.add_edge(s[c_l['g']].replace(' ','_').replace('-','_')+'_gp',s.name)
            Gloop[s[c_l['g']].replace(' ','_').replace('-','_')+'_gp'][s.name]['etype']="group_formation"
        return(Gloop)

    def mineralisation_proximity(Gloop,output_path,commodity,geology,fault,b_contacts_gdf,i_contacts_gdf,mindep,c_l):
        commodities=commodity.split(',')
        if(not 'NONE' in commodity):
            commodities.append('NONE')
        fault_tmp=fault.copy()
        b_contacts_gdf_tmp=b_contacts_gdf.copy()
        i_contacts_gdf_tmp=i_contacts_gdf.copy()
        mindep_geology = gpd.sjoin(mindep, geology, how="left", op="within")
        
        first=True
        for com in commodities:
            if(not com == 'NONE'):
                if(first):
                    mindep_sub=mindep[mindep[c_l['mscm']]==com].copy()
                    first=False
                else:
                    mindep_sub=pd.concat([mindep_sub, mindep[mindep[c_l['mscm']]==com]], sort=False)

        for ind,m in mindep_sub.iterrows():
            fault_tmp['d_'+str(ind)]=fault_tmp.distance(Point(m.geometry))

        for ind,m in mindep_sub.iterrows():
            if(len(b_contacts_gdf_tmp)>0):
                b_contacts_gdf_tmp['d_'+str(ind)]=b_contacts_gdf_tmp.distance(Point(m.geometry))

        for ind,m in mindep_sub.iterrows():
            if(len(i_contacts_gdf_tmp)>0):
                i_contacts_gdf_tmp['d_'+str(ind)]=i_contacts_gdf_tmp.distance(Point(m.geometry))

        for com in commodities:
            if(not com == 'NONE'):
                     
                mindep_com=mindep_sub[mindep_sub[c_l['mscm']]==com]
                if(len(mindep)>0):

                    bc=[]
                    for ind,b in b_contacts_gdf_tmp.iterrows():
                        basal_count=0
                        for ind2,m in mindep_com.iterrows():
                            if(b['d_'+str(ind2)]<close_b ):
                                basal_count=basal_count+1
                        bc.append(basal_count)
                    if(len(b_contacts_gdf_tmp)>0):
                        b_contacts_gdf[com+'_min']=bc    

                    fc=[]
                    for ind,f in fault_tmp.iterrows():
                        fault_count=0
                        for ind2,m in mindep_com.iterrows():
                            if(f['d_'+str(ind2)]<close_f):
                                fault_count=fault_count+1
                        fc.append(fault_count)
                    fault[com+'_min']=fc    

                    ic=[]
                    for ind,i in i_contacts_gdf_tmp.iterrows():
                        igneous_count=0
                        for ind2,m in mindep_com.iterrows():
                            if(i['d_'+str(ind2)]<close_i):
                                igneous_count=igneous_count+1
                        ic.append(igneous_count)
                    if(len(i_contacts_gdf_tmp)>0):
                        i_contacts_gdf[com+'_min']=ic 

                    for ind,b in b_contacts_gdf.iterrows():
                        Gloop.nodes[b[c_l['c']].replace(" ","_").replace("-","_")][com+'_min']=b[com+'_min']

                    for ind,f in fault.iterrows():
                        Gloop.nodes['Fault_'+f[c_l['o']]][com+'_min']=f[com+'_min']

                    for ind,i in i_contacts_gdf.iterrows():
                        Gloop.nodes[i[c_l['c']].replace(" ","_").replace("-","_")][com+'_min']=i[com+'_min']
                else:
                    for ind,b in b_contacts_gdf.iterrows():
                        Gloop.nodes[b[c_l['c']].replace(" ","_").replace("-","_")][com+'_min']=-1

                    for ind,f in fault.iterrows():
                        Gloop.nodes['Fault_'+f[c_l['o']]][com+'_min']=-1

                    for ind,i in i_contacts_gdf.iterrows():
                        Gloop.nodes[i[c_l['c']].replace(" ","_").replace("-","_")][com+'_min']=-1

                    if(len(i_contacts_gdf_tmp)>0):
                        i_contacts_gdf[com+'_min']=-1 
                    if(len(b_contacts_gdf_tmp)>0):
                        b_contacts_gdf[com+'_min']=-1 
                    fault[com+'_min']=-1
        
        if(len(i_contacts_gdf)>0):
            i_contacts_gdf.to_file(os.path.join(output_path,'igneous_contacts.shp'))
        if(len(b_contacts_gdf)>0):
            b_contacts_gdf.to_file(os.path.join(output_path,'basal_contacts.shp'))
        if(len(mindep_geology)>0):
            mindep_geology.to_file(os.path.join(output_path,'mindep_geology.shp'))  
        if(len(fault)>0):
            fault.to_file(os.path.join(output_path,'mindep_fault.shp'))
        return(Gloop)
            
               
    
    def feature_geometries(Gloop,fault,b_contacts_gdf,i_contacts_gdf,strats,c_l):

        for ind,s in strats.iterrows():
            bcontacts_strat=b_contacts_gdf[b_contacts_gdf[c_l['c']].str.replace(" ","_").replace("-","_")==s.name]
            length=0
            for ind2,b in bcontacts_strat.iterrows():
                length=length+b.geometry.length
            if(length>0):
                Gloop.nodes[s.name]['contact_length']=int(length)

        for ind,f in fault.iterrows():
            length=length+f.geometry.length
            Gloop.nodes['Fault_'+f[c_l['o']]]['fault_length']=int(length)

            x1=f.geometry.coords[0][0]
            y1=f.geometry.coords[0][1]
            x2=f.geometry.coords[-1][0]
            y2=f.geometry.coords[-1][1]
            l,m=m2l_utils.pts2dircos(x1,y1,x2,y2)
            dip,dip_dir=m2l_utils.dircos2ddd(l,m,0)
            
            Gloop.nodes['Fault_'+f[c_l['o']]]['fault_trend']=int(dip_dir%180)
        if(len(i_contacts_gdf)>0):
            for ind,s in strats.iterrows():
                icontacts_strat=i_contacts_gdf[i_contacts_gdf[c_l['c']].str.replace(" ","_").replace("-","_")==s.name]
                length=0
                for ind2,b in icontacts_strat.iterrows():
                    length=length+b.geometry.length
                if(length>0):
                    Gloop.nodes[s.name]['contact_length']=int(length)
            
        return(Gloop)


    def map2graph(output_path,geology_file,fault_file,mindep_file,c_l,deposits,fault_orientation_clusters,fault_length_clusters,
                            fault_fault_weight=3,
                            fault_weight=1,
                            formation_weight=7,
                            formation_formation_weight=9,
                            fault_formation_weight=5):
        if (not os.path.isdir(output_path)):
            os.mkdir(output_path)
        
        mindep,fault_clean,geology_clean,geology_clean_nona = Map2Graph.clean_inputs(geology_file,fault_file,mindep_file)
        groups,strats,all_contacts,Gloop = Map2Graph.strat_graph(geology_clean,output_path,c_l)
        Map2Graph.fix_Loop_graph(output_path,'pre_loop')
        Topology.use_asud(os.path.join(output_path,'pre_loop.gml'), output_path)
        Gloop,i_contacts_gdf,b_contacts_gdf = Map2Graph.basal_contacts(Gloop,output_path,geology_clean,all_contacts,groups,c_l,formation_weight)
        Gloop = Map2Graph.fault_intersections(Gloop,output_path,fault_clean,c_l,fault_weight,fault_fault_weight)
        Gloop = Map2Graph.fault_formation_intersections(Gloop,output_path,fault_clean,groups,geology_clean_nona,c_l,fault_formation_weight)
        Gloop = Map2Graph.remove_group_info(Gloop)
        Gloop = Map2Graph.group_formation_intersections(Gloop,strats,c_l)
        Gloop = Map2Graph.feature_geometries(Gloop,fault_clean,b_contacts_gdf,i_contacts_gdf,strats,c_l)   
        Gloop = Map2Graph.mineralisation_proximity(Gloop,output_path,deposits,geology_clean,fault_clean,b_contacts_gdf,i_contacts_gdf,mindep,c_l)
        Gloop = Map2Graph.classify_faults(Gloop,fault_orientation_clusters,fault_length_clusters)
        
        nx.write_gml(Gloop, os.path.join(output_path,'pre_loop.gml'))
        Map2Graph.fix_Loop_graph(output_path,'pre_loop')

        Topology.colour_Loop_graph( output_path,'pre_loop') 


    def classify_faults(Gloop,fault_orientation_clusters,fault_length_clusters):
        from sklearn import cluster    
        Af_d=Map2Graph.get_nodes_from_graph(Gloop,'fault')
        Af_d['fault_trend']=Af_d['fault_trend'].astype('float64')
        Af_d['fault_length']=Af_d['fault_length'].astype('float64')

        l,m,n=Map2Graph.convert_dd_lmn(90,Af_d['fault_trend']) 
        faults=pd.DataFrame(l)
        faults.rename(columns={"fault_trend": "l"}, inplace=True)
        faults['m']=m
        faults['n']=n    

        if(len(Af_d)>=fault_orientation_clusters):
            k_means_o = cluster.KMeans(n_clusters=fault_orientation_clusters, max_iter=50, random_state=1)
            k_means_o.fit(faults) 
            labels_o = k_means_o.labels_
            clusters_o=pd.DataFrame(labels_o, columns=['cluster_o'])
        else:
            clusters_o=''
        
        if(len(Af_d)>=fault_length_clusters):
            k_means_l = cluster.KMeans(n_clusters=fault_length_clusters, max_iter=50, random_state=1)
            k_means_l.fit(Af_d['fault_length'].to_numpy().reshape(-1, 1)) 
            labels_l = k_means_l.labels_
            clusters_l=pd.DataFrame(labels_l, columns=['cluster_l'])
        else:
            clusters_l=''

        if(len(clusters_o)):
            Af_d['cluster_o']=clusters_o['cluster_o']
        if(len(clusters_l)):
            Af_d['cluster_l']=clusters_l['cluster_l']

        for n in Gloop.nodes:
            if( "Fault" in n):
                if( Af_d.loc[n].name in Gloop.nodes):
                    #Gloop.nodes[n]['Dip']=Af_o.loc[n]['dip']
                    #Gloop.nodes[n]['DipDirection']=Af_o.loc[n]['DipDirection']
                    #Gloop.nodes[n]['DipPolarity']=Af_o.loc[n]['DipPolarity']
                    if(len(clusters_o) and False):
                        Gloop.nodes[n]['OrientationCluster']=Af_d.loc[n]['cluster_o']
                    else:
                        Gloop.nodes[n]['OrientationCluster']=-1
                    if(len(clusters_l) and False):
                        Gloop.nodes[n]['LengthCluster']=Af_d.loc[n]['cluster_l']
                    else:
                        Gloop.nodes[n]['LengthCluster']=-1
                    
        
        #add centrality measures to fault nodes
        Gf=Map2Graph.get_subgraph_from_graph_nodes(Gloop,'fault')
        Gfcc=nx.closeness_centrality(Gf)
        Gfbc=nx.betweenness_centrality(Gf)

        for n in Gloop.nodes:
            if("Fault" in n):
                if(n in Gfcc.keys()):
                    Gloop.nodes[n]['ClosenessCentrality']=Gfcc[n]
                    Gloop.nodes[n]['BetweennessCentrality']=Gfbc[n]
                else:
                    Gloop.nodes[n]['ClosenessCentrality']=-1
                    Gloop.nodes[n]['BetweennessCentrality']=-1
        return(Gloop)

    def get_nodes_from_graph(Gloop,ntype):

        nodes_all=[]
        for v in Gloop.nodes():
            if(Gloop.nodes[v]['ntype']==ntype):
                nodes_all.append(v)
        nodes=Gloop.subgraph(nodes_all)
        data=pd.DataFrame.from_dict(dict(nodes.nodes(data=True)), orient='index')
        data['name']=data.index

        return(data)

    def get_subgraph_from_graph_edges(Gloop,etype):

        edges_all=[]
        for v in Gloop.edges:
            if(Gloop.edges[v]['etype']==etype):
                edges_all.append(v)
        subgraph=nx.Graph().to_directed()
        for e in edges_all:
            if(not e[0] in subgraph.nodes()):
                subgraph.add_node(e[0],data=Gloop.nodes[e[0]])
            if(not e[1] in subgraph.nodes()):
                subgraph.add_node(e[1],data=Gloop.nodes[e[1]])
            subgraph.add_edge(e[0],e[1])
                
        return(subgraph)
    
    def convert_dd_lmn(dip, dipdir):

            r_dd=np.radians(dipdir)
            r_90_dip=np.radians(90-dip)
            sin_dd=np.sin(r_dd)
            cos_dd=np.cos(r_dd)
            
            cos_90_dip=np.cos(r_90_dip)
            sin_90_dip=np.sin(r_90_dip)
            
            l=sin_dd*cos_90_dip
            m=cos_dd*cos_90_dip
            n=sin_90_dip
            return(l, m, n)

    def get_subgraph_from_graph_nodes(Gloop,gtype):

        nodes_all=[]
        for v in Gloop.nodes():
            if(Gloop.nodes[v]['ntype']==gtype):
                nodes_all.append(v)
        subgraph=Gloop.subgraph(nodes_all)
        return(subgraph)

    def get_metadata_from_graph(Gloop):

        for v in Gloop.nodes():
            if(Gloop.nodes[v]['ntype']=='metadata'):
                run_flags=Gloop.nodes[v]['run_flags']
                c_l=Gloop.nodes[v]['c_l']
                config=Gloop.nodes[v]['config']
        return(c_l,run_flags,config)  

    def get_dst_crs_from_graph(Gloop):

        for v in Gloop.nodes():
            if(Gloop.nodes[v]['ntype']=='dst_crs'):
                data=Gloop.nodes[v]['data']
        return(data)        
        
    def get_bbox_from_graph(Gloop):

        for v in Gloop.nodes():
            if(Gloop.nodes[v]['ntype']=='bbox'):
                data=Gloop.nodes[v]['data']
                bbox=np.array(data.replace('[','').replace(']','').split(','), dtype=np.float32)
        return(bbox)        
    
    def get_dtm_from_graph(Gloop):

        for v in Gloop.nodes():
            if(Gloop.nodes[v]['ntype']=='dtm'):
                data=Gloop.nodes[v]['data']
                shape=Gloop.nodes[v]['shape']
                minx=Gloop.nodes[v]['minx']
                miny=Gloop.nodes[v]['miny']
                maxx=Gloop.nodes[v]['maxx']
                maxy=Gloop.nodes[v]['maxy']
                xscale=Gloop.nodes[v]['xscale']
                yscale=Gloop.nodes[v]['yscale']

                dtm=np.array(data.replace('[','').replace(']','').split(','), dtype=np.float32)
                shape=shape.replace('(','').replace(')','').split(',')
                x=int(shape[0])
                y=int(shape[1])
                dtm=dtm.reshape(x,y)
        return(minx,miny,maxx,maxy,xscale,yscale,dtm)        

    def get_point_data_from_graph(Gloop,datatype):

        for v in Gloop.nodes():
            if(Gloop.nodes[v]['ntype']=='points'):
                data=Gloop.nodes[v]['data']
        df = pd.DataFrame()
        if(datatype=='fault_geom'):
            Param1='Param1'
            Param2='Param2'
            Param3='Param3'
            Param4='Param4'
        elif(datatype=='fault_displacement'):
            Param1='HzRad'
            Param2='Vrad'
            Param3='NDist'
            Param4='Param4'
        elif(datatype=='fault_strat_displacement'):
            Param1='left_fm'
            Param2='right_fm'
            Param3='min_offset'
            Param4='strat_offset'
        elif(datatype=='contact'):
            Param1='Param1'
            Param2='Param2'
            Param3='Param3'
            Param4='Param4'
        elif(datatype=='orientation'):
            Param1='OCluster'
            Param2='LCluster'
            Param3='CCentral'
            Param4='BCentral'
        elif(datatype=='raw_contact'):
            Param1='group'
            Param2='angle'
            Param3='lsx'
            Param4='lsy'
            
        for d in data:
            if(d['type']==datatype):
                df = df.append({  'type': d['type'],
                                'name': d['name'],
                                'X': d['X'],
                                'Y': d['Y'],
                                'Z': d['Z'],
                                Param1:d['Param1'],
                                Param2:d['Param2'],
                                Param3:d['Param3'],
                                Param4:d['Param4']
                            }, ignore_index=True)
        return(df)

    def  granular_map2graph(output_path,geology_file,fault_file,mindep_file,c_l,deposits,
                            fault_fault_weight=3,
                            fault_weight=1,
                            formation_weight=7,
                            formation_formation_weight=9,
                            fault_formation_weight=5):
        if (not os.path.isdir(output_path)):
            os.mkdir(output_path)

        mindep,fault_clean,geology_clean,geology_clean_nona = Map2Graph.clean_inputs(geology_file,fault_file,mindep_file)
        geology_exploded=Map2Graph.granular_explode_geology(geology_clean)
        Map2Graph.granular_strat_graph(output_path,geology_exploded,c_l,formation_weight,formation_formation_weight)
        Map2Graph.fix_Loop_graph(output_path,'granular_pre_loop')
        Gloop=nx.read_gml(os.path.join(output_path,'granular_pre_loop.gml'))
        Gloop = Map2Graph.fault_intersections(Gloop,output_path,fault_clean,c_l,fault_weight,fault_fault_weight)
        Gloop = Map2Graph.granular_fault_formation_intersections(Gloop,fault_clean,geology_exploded,c_l,fault_formation_weight)
        Map2Graph.granular_mineralisation_proximity(Gloop,output_path,fault_clean,geology_exploded,mindep,deposits,c_l)

    def granular_explode_geology(geology):

        geology_exploded=geology.explode()
        geology_exploded.index= pd.RangeIndex(start=0, stop=len(geology_exploded), step=1)
        geology_exploded['idx']=geology_exploded.index
        geology_exploded.crs=geology.crs
        return(geology_exploded)
    
    def granular_strat_graph(output_path,geology_exploded,c_l,formation_weight,formation_formation_weight):

        all_contacts=[]
        for ind,g in geology_exploded.iterrows():
            for ind2,g2 in geology_exploded.iterrows():
                if(not ind >= ind2):
                    if(g.geometry.intersects(g2.geometry)):
                        g1_snapped = snap(g.geometry, g2.geometry,snap_buffer)
                        all_contacts.append([ind,ind2,g1_snapped.buffer(0).intersection(g2.geometry.buffer(0))])

        Gloop=nx.DiGraph()

        index=0

        for ind,p in geology_exploded.iterrows():
            Gloop.add_node(ind,id=index, formation=p[c_l['c']].replace(" ","_").replace("-","_"),LabelGraphics ='[ text "'+geology_exploded.iloc[ind][c_l['c']].replace(' ','_').replace('-','_')+'_'+str(ind)+'" anchor "n" fontStyle "bold" fontSize 14 ]')
            Gloop.nodes[ind]['ntype']='formation'
            Gloop.nodes[ind]['weight']=formation_weight
            index=index+1


        for c in all_contacts:
            ave_age0=int(geology_exploded.iloc[c[0]][c_l['min']])+((int(geology_exploded.iloc[c[0]][c_l['max']])-int(geology_exploded.iloc[c[0]][c_l['min']]))/2)
            ave_age1=int(geology_exploded.iloc[c[1]][c_l['min']])+((int(geology_exploded.iloc[c[1]][c_l['max']])-int(geology_exploded.iloc[c[1]][c_l['min']]))/2)

            if(ave_age0>ave_age1):
                Gloop.add_edge(geology_exploded.iloc[c[1]].name,geology_exploded.iloc[c[0]].name)
                Gloop[geology_exploded.iloc[c[1]].name][geology_exploded.iloc[c[0]].name]['etype']='formation_formation'
                Gloop[geology_exploded.iloc[c[1]].name][geology_exploded.iloc[c[0]].name]['formation1']=geology_exploded.iloc[c[0]][c_l['c']].replace(' ','_').replace('-','_')+'_'+str(geology_exploded.iloc[c[0]]['idx'])
                Gloop[geology_exploded.iloc[c[1]].name][geology_exploded.iloc[c[0]].name]['formation2']=geology_exploded.iloc[c[1]][c_l['c']].replace(' ','_').replace('-','_')+'_'+str(geology_exploded.iloc[c[1]]['idx'])
                Gloop[geology_exploded.iloc[c[1]].name][geology_exploded.iloc[c[0]].name]['weight']=formation_formation_weight
            else:
                Gloop.add_edge(geology_exploded.iloc[c[0]].name,geology_exploded.iloc[c[1]].name)
                Gloop[geology_exploded.iloc[c[0]].name][geology_exploded.iloc[c[1]].name]['etype']='formation_formation'
                Gloop[geology_exploded.iloc[c[0]].name][geology_exploded.iloc[c[1]].name]['formation1']=geology_exploded.iloc[c[1]][c_l['c']].replace(' ','_').replace('-','_')+'_'+str(geology_exploded.iloc[c[1]]['idx'])
                Gloop[geology_exploded.iloc[c[0]].name][geology_exploded.iloc[c[1]].name]['formation2']=geology_exploded.iloc[c[0]][c_l['c']].replace(' ','_').replace('-','_')+'_'+str(geology_exploded.iloc[c[0]]['idx'])
                Gloop[geology_exploded.iloc[c[0]].name][geology_exploded.iloc[c[1]].name]['weight']=formation_formation_weight
                
        nx.write_gml(Gloop, os.path.join(output_path,'granular_pre_loop.gml'))

        """
        geology_clean=geology_exploded.copy()
        igneous_contacts={}
        not_igneous_contacts={}
        i=0
        for c in all_contacts:
            if(c[2].geom_type=='MultiLineString' or c[2].geom_type=='LineString' ): 
                if( c_l['intrusive'] in geology_clean.iloc[c[0]][c_l['r1']] ):
                    igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[0]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[1]][c_l['c']], "geometry": c[2]}
                elif(c_l['intrusive'] in geology_clean.iloc[c[1]][c_l['r1']] ):
                    igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[1]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[0]][c_l['c']], "geometry": c[2]}
                else:
                    if(mini_strat_df.loc[geology_clean.iloc[c[0]][c_l['c']].replace(" ","_").replace("-","_")]['order']>
                    mini_strat_df.loc[geology_clean.iloc[c[1]][c_l['c']].replace(" ","_").replace("-","_")]['order']):
                        not_igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[1]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[0]][c_l['c']], "geometry": c[2]}
                    else:
                        not_igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[0]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[1]][c_l['c']], "geometry": c[2]}
                        
                i=i+1
            elif(c[2].geom_type=='GeometryCollection' ):
                for geom in c[2]:
                    if(geom.geom_type=='MultiLineString' or geom.geom_type=='LineString' ):
                        if( c_l['intrusive'] in geology_clean.iloc[c[0]][c_l['r1']] ):
                            igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[0]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[1]][c_l['c']], "geometry": geom}
                        elif(c_l['intrusive'] in geology_clean.iloc[c[1]][c_l['r1']] ):
                            igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[1]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[0]][c_l['c']], "geometry": geom}
                        else:
                            if(mini_strat_df.loc[geology_clean.iloc[c[0]][c_l['c']].replace(" ","_").replace("-","_")]['order']>
                            mini_strat_df.loc[geology_clean.iloc[c[1]][c_l['c']].replace(" ","_").replace("-","_")]['order']):
                                not_igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[1]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[0]][c_l['c']], "geometry": geom}
                            else:
                                not_igneous_contacts[i] = {"id": i, c_l['c']: geology_clean.iloc[c[0]][c_l['c']], c_l['c']+'2': geology_clean.iloc[c[1]][c_l['c']], "geometry": geom}
                        i=i+1
            
        df = DataFrame.from_dict(igneous_contacts, "index")
        if(len(df)>0):
            i_contacts_gdf = GeoDataFrame(df, crs=geology_clean.crs, geometry='geometry')
        else:
            i_contacts_gdf=DataFrame.from_dict({}, "index")

        df = DataFrame.from_dict(not_igneous_contacts, "index")
        if(len(df)>0):
            b_contacts_gdf = GeoDataFrame(df, crs=geology_clean.crs, geometry='geometry')
        else:
            b_contacts_gdf=DataFrame.from_dict({}, "index")

        return(Gloop,i_contacts_gdf,b_contacts_gdf)
        """

    def granular_mineralisation_proximity(Gloop,output_path,fault,geology_exploded,mindep,commodity,c_l):
        commodities=commodity.split(',')
        if(not 'NONE' in commodity):
            commodities.append('NONE')
        fault_tmp=fault.copy()
        geology_exploded_tmp=geology_exploded.copy()
        geology_exploded_tmp.crs=geology_exploded.crs
        mindep_geology = gpd.sjoin(mindep, geology_exploded, how="left", op="within")
        
        for com in commodities:
            if(not com == 'NONE'):
                mindep_com=mindep_geology[mindep_geology[c_l['mscm']]==com]

                if(len(mindep_com)>0):
                
                    for ind,m in mindep_com.iterrows():
                        fault_tmp['d_'+str(ind)]=fault_tmp.distance(Point(m.geometry))

                    gc=[]
                    for ind,g in geology_exploded_tmp.iterrows():
                        fm_count=0
                        for ind2,m in mindep_com.iterrows():
                            if(m[c_l['c']]==g[c_l['c']]):
                                fm_count=fm_count+1
                        gc.append(fm_count)
                    geology_exploded[com+'_min']=gc    

                    fc=[]
                    for ind,f in fault_tmp.iterrows():
                        fault_count=0
                        for ind2,m in mindep_com.iterrows():
                            if(f['d_'+str(ind2)]<close_f):
                                fault_count=fault_count+1
                        fc.append(fault_count)
                    fault[com+'_min']=fc    

                    for ind,b in geology_exploded.iterrows():
                        Gloop.nodes[str(ind)][com+'_min']=b[com+'_min']

                    for ind,f in fault.iterrows():
                        Gloop.nodes['Fault_'+f[c_l['o']]][com+'_min']=f[com+'_min']

                else:
                    for ind,b in geology_exploded.iterrows():
                        Gloop.nodes[str(ind)][com+'_min']=-1

                    for ind,f in fault.iterrows():
                        Gloop.nodes['Fault_'+f[c_l['o']]][com+'_min']=-1
                    
                    geology_exploded[com+'_min']=-1    
                    fault[com+'_min']=-1    

                
                
        if(len(geology_exploded)>0):
            geology_exploded.to_file(os.path.join(output_path,'exploded_geology_min.shp'))
        if(len(mindep_geology)>0):
            mindep_geology.to_file(os.path.join(output_path,'mindep_geology.shp'))    
        if(len(fault)>0):
            fault.to_file(os.path.join(output_path,'mindep_fault.shp')) 
            
        nx.write_gml(Gloop, os.path.join(output_path,'granular_pre_loop_mindep.gml'))
        Topology.colour_Loop_graph( output_path,'granular_pre_loop_mindep')
        Map2Graph.fix_Loop_graph(output_path,'granular_pre_loop_mindep_colour')


    def granular_fault_formation_intersections(Gloop,fault_clean,geology_exploded,c_l,fault_formation_weight):

        for ind,f in fault_clean.iterrows():
            for ind2,g in geology_exploded.iterrows():
                if(f.geometry.intersects(g.geometry.buffer(snap_buffer).buffer(0))):
                    g1_snapped = snap(g.geometry.buffer(snap_buffer).buffer(0), f.geometry.buffer(0),snap_buffer)
                    x=g1_snapped.intersection(f.geometry)
                    if(x.geom_type!='Point' and x.geom_type!='MultiPoint' ):
                        Gloop.add_edge('Fault_'+f[c_l['o']],str(g.name))
                        Gloop['Fault_'+f[c_l['o']]][str(g.name)]['etype']='fault_formation'
                        Gloop['Fault_'+f[c_l['o']]][str(g.name)]['weight']=fault_formation_weight
                        Gloop['Fault_'+f[c_l['o']]][str(g.name)]['fault']='Fault_'+f[c_l['o']]
                        Gloop['Fault_'+f[c_l['o']]][str(g.name)]['formation']=str(g[c_l['c']].replace(" ","_").replace("-","_"))+'_'+str(g['idx'])


        return(Gloop)

    def get_dijkstra_path(Gloop,source,target):
        Gloopu=Gloop.to_undirected()
        def func(u, v, d):
            node_u_wt = Gloopu.nodes[u].get("weight", 1)
            node_v_wt = Gloopu.nodes[v].get("weight", 1)
            edge_wt = d.get("weight", 1)
            return node_u_wt / 2 + node_v_wt / 2 + edge_wt
        path = nx.dijkstra_path(Gloopu, source=source, target=target,weight=func)
        return(list(path))

    def get_resistance_path_length(Gloop,source,target):
        Gloopu=Gloop.to_undirected()
        return(nx.resistance_distance(Gloopu, nodeA=source, nodeB=target,weight='weight'))