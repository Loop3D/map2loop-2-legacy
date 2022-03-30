import os
import netCDF4
import numpy as np
from geophys_utils import NetCDFGridUtils
from geophys_utils import get_netcdf_edge_points, points2convex_hull
import matplotlib.pyplot as plt
import time
import map2loop.SpaceTime as st
from math import sqrt, pi
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from geopandas import GeoDataFrame
from math import floor

def plot_profile(petro_data,geophys_data,code,title,axis):
    #fig, ax_left = plt.subplots()
    ax_right = axis.twinx()
    axis.title.set_text(title)

    axis.plot(petro_data[code], color='red')
    ax_right.plot(geophys_data, color='blue')
    #plt.show()
    
def plot_petrophys_map(geol_clip,colour_col,axis):
    geol_clip.plot(ax=axis,column=colour_col, cmap='gist_rainbow')
    axis.title.set_text(colour_col+' Strat')


def plot_den_vs_ms(petrophysics,axis):
    axis.scatter(petrophysics['density'], petrophysics['ms'], c='blue', alpha=0.5)
    axis.title.set_text('density vs log magsus 10e-5')
    axis.set_ylabel('log mag sus')
    axis.set_xlabel('density')
    #plt.show()

def plot_dtm_geol(dtm_profile,geol_clip,mag_long,mag_lat,axis):
    dtm_profile=dtm_profile[~dtm_profile['dtm'].isna()]
    dtm_profile=dtm_profile.to_crs('EPSG:28350')
    base=geol_clip.plot(ax=axis,edgecolor='gray',color='w',figsize=(18, 12))
    axis.title.set_text('dtm')

    dtm_profile.plot(ax=base, column='dtm', markersize=5,cmap='jet')
    #plt.show()

def plot_geol_profile(geol_clip,geol_code,axis):
    base=geol_clip.plot(ax=axis,edgecolor='gray',color='w',figsize=(18, 12))
    axis.title.set_text('geol')
    plt.ylabel('y')
    plt.xlabel('x')

    geol_code.plot(ax=base, column='descriptn', markersize=5)
    #plt.show()

def plot_grid(grid,axis,title):
    axis.imshow(grid,cmap ='gist_rainbow',vmin=np.percentile(grid, 1.5), vmax=np.percentile(grid, 98.5))
    axis.title.set_text(title)
    #plt.show()

def plot_density(density,axis):
    if(len(density)>0):
        density.plot(ax=axis,column='value')
    else:
        plt.text(1.2,.2,"No density values\n in area")
    plt.tight_layout()
    axis.title.set_text('Density Map')

def plot_mag_sus(mag_sus,axis):
    if(len(mag_sus)>0):
        mag_sus.plot(ax=axis,color='blue')
    else:
        plt.text(.2,.2,"No mag sus values\n in area")
    plt.tight_layout()
    axis.title.set_text('Mag Sus Map')

def plot_geol(geol_clip,colour_col,axis):
    geol_clip.plot(ax=axis,column=colour_col, cmap='tab20c')
    axis.title.set_text('geology')

def calc_profile(data,res,bbox):
    
    len=int(sqrt((bbox[0]-bbox[2])**2+(bbox[1]-bbox[3])**2)/res)
    
    long=np.linspace(bbox[0],bbox[2],len)
    lat= np.linspace(bbox[1],bbox[3],len)
    
    xr=(bbox[2]-bbox[0])/(data.shape[0]-1)
    yr=(bbox[3]-bbox[1])/(data.shape[1]-1)

    ll=zip(long,lat)
    values=[]
    for x,y in ll:
        r=int((x-bbox[0])/xr)   
        c=int((y-bbox[1])/yr)
        values.append([data[r,c]])
    
    return(values,long,lat)

def get_density(density_filename,bbox):
    density=gpd.read_file(density_filename,bbox=bbox)
    return(density)

def get_mag_sus(mag_sus_filename,bbox):
    mag_sus=gpd.read_file(mag_sus_filename,bbox=bbox)
    return(mag_sus)

def get_dtm_profile(minlong,minlat,maxlong,maxlat,mag_long,mag_lat,l):
    from owslib.wcs import WebCoverageService
    import rasterio
    from map2loop.m2l_utils import value_from_raster

    url="http://services.ga.gov.au/gis/services/DEM_SRTM_1Second_over_Bathymetry_Topography/MapServer/WCSServer?"
    dtm_file='./dtm.tif'
    bbox = (minlong-.1, minlat-.1, maxlong+.1, maxlat+.1)
    wcs = WebCoverageService(url, version='1.0.0')

    cvg = wcs.getCoverage(identifier='1',  bbox=bbox,
                format='GeoTIFF', crs=4326, width=200, height=200)
    f = open(dtm_file, 'wb')
    bytes_written = f.write(cvg.read())
    f.close()
    print("dtm geotif saved as", dtm_file)
    
    with rasterio.open(dtm_file) as dtm:
        band1 = dtm.read(1)

        dtm_data=[]
        for i in range(len(mag_long)):
            locations=[(mag_long[i],mag_lat[i])]
            height=float(value_from_raster(dtm, locations))
            dtm_data.append([height])


    df=pd.DataFrame(dtm_data)

    dtm_profile = gpd.GeoDataFrame(index=np.arange(mag_long.shape[0]), crs='EPSG:4326', geometry=l) 
    dtm_profile['dtm']=df
    dtm_profile=dtm_profile[~dtm_profile['dtm'].isna()]
    dtm_profile=dtm_profile.to_crs('EPSG:28350')
    return(dtm_profile)

def init_petrophysics(petrophysics,geology):
    for g in list(set(geology['unitname'])):
        petrophysics[g]={'den_samples':0,
                     'density':-99.99,
                     'den_std':-99.99,
                     'den_unit':'kg m3',
                     'den_source':'0',
                     'ms_samples':0,
                     'ms':-99.99,
                     'ms_std':-99.99,
                     'ms_unit':'log(SI unit x 10-5)',
                     'ms_source':'0'
                    }
    petrophysics=pd.DataFrame.from_dict(petrophysics,orient='index')
    return(petrophysics)

def GA_point_petrophysics(petrophysics,mag_sus,density):
    for g in list(set(mag_sus['stratigrap'])):
        g_ms=mag_sus[mag_sus['stratigrap']==g]
        g_ms=g_ms[g_ms['value']>0]
        g_ms_5=g_ms[g_ms['unitofmeas']=='SI unit x 10-5']
        g_ms_3=g_ms[g_ms['unitofmeas']=='SI unit x 10-3']
        g_ms_c=g_ms[g_ms['unitofmeas']=='centimetre gram seconds x 10-6']
        g_ms_3b=g_ms_3.copy()
        g_ms_3b['value']=g_ms_3b['value']*100
        g_ms_cb=g_ms_c.copy()
        g_ms_cb['value']=g_ms_c['value']*4*pi/10 #check this conversion!!
        g_ms=pd.concat([g_ms_5,g_ms_3b,g_ms_cb],sort=False)
        if(len(g_ms)>2 and g in petrophysics.index):
            petrophysics.at[g,'ms']=np.log(g_ms['value']+0.00001).mean()
            petrophysics.at[g,'ms_samples']=len(g_ms)
            if(len(g_ms)==1):
                petrophysics.at[g,'ms_std']=-1            
            else:
                petrophysics.at[g,'ms_std']=np.log(g_ms['value']+0.00001).std()
            petrophysics.at[g,'ms_source']='local'

    for g in list(set(density['stratigrap'])):
        g_den=density[density['stratigrap']==g]

        if(len(g_den)>2 and g in petrophysics.index):
            petrophysics.at[g,'den_samples']=len(g_den)
            petrophysics.at[g,'density']=g_den['value'].mean()
            petrophysics.at[g,'den_std']=g_den['value'].std()
            petrophysics.at[g,'den_source']='local'

            if(len(g_ms)==1):
                petrophysics.at[g,'ms_std']=-1
    return(petrophysics)

def GA_strat_petrophysics(petrophysics,geology,mag_sus_filename_local,density_filename_local):
    geology_strat=list(set(geology['unitname']))
    GA_ms=gpd.read_file(mag_sus_filename_local)
    GA_den=gpd.read_file(density_filename_local)

    for g in geology_strat:
        g_strat=geology[geology['unitname']==g]
        g_den=GA_den[GA_den['stratigrap']==g]
        g_ms=GA_ms[GA_ms['stratigrap']==g]
        g_ms=g_ms[g_ms['value']>0]
        g_ms_5=g_ms[g_ms['unitofmeas']=='SI unit x 10-5']
        g_ms_3=g_ms[g_ms['unitofmeas']=='SI unit x 10-3']
        g_ms_c=g_ms[g_ms['unitofmeas']=='centimetre gram seconds x 10-6']
        g_ms_3b=g_ms_3.copy()
        g_ms_3b['value']=g_ms_3b['value']*100
        g_ms_cb=g_ms_c.copy()
        g_ms_cb['value']=g_ms_c['value']*4*pi/10
        g_ms=pd.concat([g_ms_5,g_ms_3b,g_ms_cb],sort=False)

        if((len(g_den)>0 or len(g_ms)>0 ) and g in petrophysics.index):
            if(petrophysics.loc[g]['den_samples']==0):
                if(len(g_den)>0):
                    petrophysics.at[g,'density']=g_den['value'].mean()
                    petrophysics.at[g,'den_samples']=len(g_den)
                    if(len(g_den)==1):
                        petrophysics.at[g,'den_std']=-1            
                    else:
                        petrophysics.at[g,'den_std']=g_den['value'].std()
                    petrophysics.at[g,'den_source']='strat'
            if(petrophysics.loc[g]['ms_samples']==0):
                if(len(g_ms)>0):
                    petrophysics.at[g,'ms']=np.log(g_ms['value']+0.00001).mean()
                    petrophysics.at[g,'ms_samples']=len(g_ms)
                    if(len(g_ms)==1):
                        petrophysics.at[g,'ms_std']=-1            
                    else:
                        petrophysics.at[g,'ms_std']=np.log(g_ms['value']+0.00001).std()
                    petrophysics.at[g,'ms_source']='strat'
    return(petrophysics,GA_ms,GA_den,g_ms)

def get_litho_from_strat(petrophysics,ens_file):
    strat_dic_gswa=pd.read_csv(ens_file)

    #get strat for which no den or ms available from GA
    no_den=petrophysics[petrophysics['den_samples']==0]
    no_ms=petrophysics[petrophysics['ms_samples']==0]
    no_den_strat=list(set(no_den.index))
    no_ms_strat=list(set(no_ms.index))

    #split list of lithologies
    strat_dic_gswa['exploded']=strat_dic_gswa['DESCRIPTN_DH2LOOP'].str.split(",")

    #get lithologies from strat wiht no petrophysics
    strat_dic_gswa2=strat_dic_gswa.drop_duplicates('UNITNAME')
    strat_dic_gswa2=strat_dic_gswa2.set_index('UNITNAME')
    return(strat_dic_gswa2,no_den_strat,no_ms_strat)

def GA_density_from_litho(petrophysics,strat_dic_gswa2,GA_den,no_den_strat):
    GA_den=GA_den[~GA_den['lithology'].isna()]
    for s in no_den_strat:
        first=True
        if(s in strat_dic_gswa2.index):
            for ls in list(set(strat_dic_gswa2.loc[s]['exploded'])):
                GA_den_l=GA_den[GA_den['lithology'].str.contains(ls)]
                if(first):
                    GA_den_l_all=GA_den_l.copy()
                    first=False
                else:
                    GA_den_l_all=pd.concat([GA_den_l_all,GA_den_l],sort=False)


            if(len(GA_den_l_all)>0):
                petrophysics.at[s,'density']=GA_den_l_all['value'].mean()
                petrophysics.at[s,'den_samples']=len(GA_den_l_all)
                if(len(GA_den_l_all)==1):
                    petrophysics.at[s,'den_std']=-1            
                else:
                    petrophysics.at[s,'den_std']=GA_den_l_all['value'].std()
                petrophysics.at[s,'den_source']='litho'
            else:
                petrophysics.at[s,'density']=2.7
                petrophysics.at[s,'den_samples']=0
                petrophysics.at[s,'den_std']=-1
                petrophysics.at[s,'den_source']='crustal'
                
    return(petrophysics)

def GA_magsus_from_litho(petrophysics,strat_dic_gswa2,GA_ms,no_ms_strat,g_ms):
    GA_ms=GA_ms[~GA_ms['lithology'].isna()]
    for s in no_ms_strat:
        first=True
        if(s in strat_dic_gswa2.index):
            for ls in list(set(strat_dic_gswa2.loc[s]['exploded'])):

                GA_ms_l=GA_ms[GA_ms['lithology'].str.contains(ls) & GA_ms['value']>0 ]
                g_ms_5=GA_ms_l[GA_ms_l['unitofmeas']=='SI unit x 10-5']
                g_ms_3=GA_ms_l[GA_ms_l['unitofmeas']=='SI unit x 10-3']
                g_ms_c=g_ms[g_ms['unitofmeas']=='centimetre gram seconds x 10-6']
                g_ms_3b=g_ms_3.copy()
                g_ms_3b['value']=g_ms_3b['value']*100
                g_ms_cb=g_ms_c.copy()
                g_ms_cb['value']=g_ms_c['value']*4*pi/10
                GA_ms_l=pd.concat([g_ms_5,g_ms_3b,g_ms_cb],sort=False)            
                if(first):
                    GA_ms_l_all=GA_ms_l.copy()
                    first=False
                else:
                    if(len(GA_ms_l)==0):
                        first=True
                    else:    
                        GA_ms_l_all=pd.concat([GA_ms_l_all,GA_ms_l],sort=False)

            if(len(GA_ms_l_all)>0):
                petrophysics.at[s,'ms']=np.log(np.abs(GA_ms_l_all['value'])).mean()
                petrophysics.at[s,'ms_samples']=len(GA_ms_l_all)
                if(len(GA_ms_l_all)==1):
                    petrophysics.at[s,'ms_std']=-1            
                else:
                    petrophysics.at[s,'ms_std']=np.log(np.abs(GA_ms_l_all['value'])).std()
                petrophysics.at[s,'ms_source']='litho'
            else:
                petrophysics.at[s,'ms']=0
                petrophysics.at[s,'ms_samples']=0
                petrophysics.at[s,'ms_std']=-1
                petrophysics.at[s,'ms_source']='crustal'

    return(petrophysics)

def petrophysics_profile(petrophysics,geology,geol_clip,gdf,grav_long,grav_lat):

    petrophysics['unitname']=petrophysics.index
    geol_clip=geol_clip.merge(petrophysics,on='unitname')
    petro_profile_mag=gpd.sjoin(gdf, geol_clip, how="left", predicate="within")
    grav_ll=zip(grav_long,grav_lat)
    l_den = [Point(i[0],i[1]) for i in grav_ll]

    #l_den = [Point(i[0],i[1]) for i in sample_points]

    gdf_den = GeoDataFrame(np.zeros(len(grav_long)), crs='EPSG:4326', geometry=l_den)
    gdf_den=gdf_den.to_crs('EPSG:28350')
    geol_code = gpd.sjoin(gdf_den, geology, how="left", predicate="within")

    petro_profile_den=gpd.sjoin(gdf_den, geol_clip, how="left", predicate="within")
    
    return(petrophysics,geol_clip,petro_profile_mag,petro_profile_den,geol_code)

