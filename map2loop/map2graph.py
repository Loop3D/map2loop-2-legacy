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
from .m2l_utils import display, print


close_f=2000
close_b=2000
close_i=2000
fault_clip_buffer=10
snap_buffer=.5


class Map2Graph(object):

    def __init__(self):
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
        groups.reset_index(inplace=True)
        groups['order'] = groups.index
        groups.set_index([c_l['g']],inplace=True)
            
        strats=geology_clean.drop_duplicates(subset=[c_l['c']])
        strats[c_l['c']]=strats[c_l['c']].str.replace(' ','_').replace('-','_')
        strats.set_index([c_l['c']],inplace=True)

        Gloop=nx.DiGraph()

        i=[]
        index=0
        for ind,g in groups.iterrows():
            Gloop.add_node(g.name.replace(" ","_").replace("-","_")+"_gp",id=index, isGroup=1,graphics= '[ fill "#FAFAFA" ]',LabelGraphics ='[ text "'+g.name.replace(" ","_").replace("-","_")+"_gp"+'" anchor "n" fontStyle "bold" fontSize 14 ]')
            Gloop.nodes[g.name.replace(" ","_").replace("-","_")+"_gp"]['ntype']='group'
            index=index+1
            i.append(index)
        groups['gid']=i

        for ind,s in strats.iterrows():
            Gloop.add_node(s.name,id=index, gid=str(groups.loc[s[c_l['g']]]['gid']),graphics= '[ fill "#FAFAFA" ]',LabelGraphics ='[ text "'+s.name.replace(" ","_").replace("-","_")+'" anchor "n" fontStyle "bold" fontSize 14 ]')
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
            new_graph=open(os.path.join(output_path,prefix+'.gml'),"w")
            for l in lines:
                l=l.replace('&#34;','"').replace('"[','[').replace(']"',']')
                if('gid "' in l):
                    tmp=l.split('"')
                    new_graph.write('    gid '+str(int(tmp[1])-1)+'\n')                
                else:
                    new_graph.write(l)
            new_graph.close()

        
    
    def basal_contacts(output_path,geology_clean,all_contacts,groups,c_l):
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
        i_contacts_gdf = GeoDataFrame(df, crs=geology_clean.crs, geometry='geometry')

        df = DataFrame.from_dict(not_igneous_contacts, "index")
        b_contacts_gdf = GeoDataFrame(df, crs=geology_clean.crs, geometry='geometry')

        return(i_contacts_gdf,b_contacts_gdf)

    def fault_intersections(output_path,fault_clean,c_l):
        Gloop=nx.read_gml( os.path.join(output_path,'pre_loop.gml'))
        Gfault=nx.DiGraph()
        fault_intersections=[]

        for ind,f in fault_clean.iterrows():
            Gloop.add_node('Fault_'+f[c_l['o']])
            Gloop.nodes['Fault_'+f[c_l['o']]]['ntype']='fault'
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
                            Gloop['Fault_'+f.name]['Fault_'+other]['angle']=int(ang)
                            Gloop['Fault_'+f.name]['Fault_'+other]['topol']=ang_code
                            
                            Gfault.add_edge('Fault_'+f.name,'Fault_'+other)
                            Gfault['Fault_'+f.name]['Fault_'+other]['etype']='fault_fault'
                            Gfault['Fault_'+f.name]['Fault_'+other]['angle']=int(ang)
                            Gfault['Fault_'+f.name]['Fault_'+other]['topol']=ang_code
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
                            Gloop['Fault_'+f.name]['Fault_'+other]['angle']=int(ang)
                            Gloop['Fault_'+f.name]['Fault_'+other]['topol']=ang_code
                            
                            Gfault.add_edge('Fault_'+f.name,'Fault_'+other)
                            Gfault['Fault_'+f.name]['Fault_'+other]['etype']='fault_fault'
                            Gfault['Fault_'+f.name]['Fault_'+other]['angle']=int(ang)
                            Gfault['Fault_'+f.name]['Fault_'+other]['topol']=ang_code

                        ind2=ind2+1

            nx.write_gml(Gfault, os.path.join(output_path,'fault_network.gml'))

            return(Gloop)
    
    def fault_formation_intersections(Gloop,output_path,fault_clean,groups,geology_clean_nona,c_l):
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
            b_contacts_gdf_tmp['d_'+str(ind)]=b_contacts_gdf_tmp.distance(Point(m.geometry))

        for ind,m in mindep_sub.iterrows():
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

        for ind,s in strats.iterrows():
            icontacts_strat=i_contacts_gdf[i_contacts_gdf[c_l['c']].str.replace(" ","_").replace("-","_")==s.name]
            length=0
            for ind2,b in icontacts_strat.iterrows():
                length=length+b.geometry.length
            if(length>0):
                Gloop.nodes[s.name]['contact_length']=int(length)
            
        return(Gloop)


    def map2graph(output_path,geology_file,fault_file,mindep_file,c_l,deposits):
        if (not os.path.isdir(output_path)):
            os.mkdir(output_path)

        mindep,fault_clean,geology_clean,geology_clean_nona = Map2Graph.clean_inputs(geology_file,fault_file,mindep_file)
        groups,strats,all_contacts,Gloop = Map2Graph.strat_graph(geology_clean,output_path,c_l)
        Map2Graph.fix_Loop_graph(output_path,'pre_loop')
        Topology.use_asud(os.path.join(output_path,'pre_loop.gml'), output_path)
        i_contacts_gdf,b_contacts_gdf = Map2Graph.basal_contacts(output_path,geology_clean,all_contacts,groups,c_l)
        Gloop = Map2Graph.fault_intersections(output_path,fault_clean,c_l)
        Gloop = Map2Graph.fault_formation_intersections(Gloop,output_path,fault_clean,groups,geology_clean_nona,c_l)
        Gloop = Map2Graph.remove_group_info(Gloop)
        Gloop = Map2Graph.group_formation_intersections(Gloop,strats,c_l)
        Gloop = Map2Graph.feature_geometries(Gloop,fault_clean,b_contacts_gdf,i_contacts_gdf,strats,c_l)   
        Gloop = Map2Graph.mineralisation_proximity(Gloop,output_path,deposits,geology_clean,fault_clean,b_contacts_gdf,i_contacts_gdf,mindep,c_l)
        
        nx.write_gml(Gloop, os.path.join(output_path,'pre_loop.gml'))
        Topology.colour_Loop_graph( output_path,'pre_loop') 