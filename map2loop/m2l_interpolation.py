import numpy as np
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from matplotlib import cm
from math import sin, cos, atan, atan2, asin, radians, degrees, sqrt, pow, acos, fabs, tan, isnan
import numpy as np
import geopandas as gpd
from geopandas import GeoDataFrame
import pandas as pd
import os
from shapely.geometry import Polygon, LineString, Point
from map2loop import m2l_utils
import rasterio

######################################
# inspired by https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
# 
# Simple Inverse Distance Weighting interpolation of observations z at x,y locations returned at locations defined by xi,yi arrays. From...
# scipy_idw(x, y, z, xi, yi)
# Args:
# x,y coordinates of points to be interpolated
# z value to be interpolated
# xi,yi grid of points where interpolation of z will be calculated - sci_py version of Simple Inverse Distance Weighting interpolation of observations z at x,y locations returned at locations defined by xi,yi arrays
# 
# simple Inverse Distance Weighting calculation
######################################
def simple_idw(x, y, z, xi, yi):
    dist = distance_matrix(x,y, xi,yi)

    # In IDW, weights are 1 / distance
    weights = 1.0 / (dist)

    # Make weights sum to one
    weights /= weights.sum(axis=0)

    # Multiply the weights for each interpolated point by all observed Z-values
    zi = np.dot(weights.T, z)
    return zi

######################################
# call scipy inverse distance weighting
# 
# Simple Inverse Distance Weighting interpolation of observations z at x,y locations returned at locations defined by xi,yi arrays. From...
# scipy_idw(x, y, z, xi, yi)
# Args:
# x,y coordinates of points to be interpolated
# z value to be interpolated
# xi,yi grid of points where interpolation of z will be calculated - sci_py version of Simple Inverse Distance Weighting interpolation of observations z at x,y locations returned at locations defined by xi,yi arrays
# 
######################################
def scipy_idw(x, y, z, xi, yi):
    interp = Rbf(x, y, z, function='linear')
    return interp(xi, yi)

######################################
# call scipy Radial basis function interpolation
# 
# scipy_rbf(x, y, z, xi, yi)
# Args:
# x,y coordinates of points to be interpolated
# z value to be interpolated
# xi,yi grid of points where interpolation of z will be calculated
# 
# sci_py version of Radial Basis Function interpolation of observations z at x,y locations returned at locations defined by xi,yi arraysplot(x,y,z,grid)
######################################
def scipy_rbf(x, y, z, xi, yi):
    interp = Rbf(x, y, z,function='multiquadric',smooth=.15)
    return interp(xi, yi)

######################################
# calculate all distances between to arrays of points
# Make a distance matrix between pairwise observations
# Note: from <http://stackoverflow.com/questions/1871536>
# (Yay for ufuncs!)
# distance_matrix(x0, y0, x1, y1)
# Args:
# x0,y0 array of point locations
# x1,y1 second array of point locations
# 
# Returns array of distances between all points defined by arrays by x0,y0 and all points defined by arrays x1,y1 from http://stackoverflow.com/questions/1871536
######################################
def distance_matrix(x0, y0, x1, y1):
    obs = np.vstack((x0, y0)).T
    interp = np.vstack((x1, y1)).T

    d0 = np.subtract.outer(obs[:,0], interp[:,0])
    d1 = np.subtract.outer(obs[:,1], interp[:,1])

    return np.hypot(d0, d1)

######################################
# plot an array of data
######################################
def plot(x,y,z,grid):
    plt.figure()
    plt.imshow(grid, extent=(0,100,0,100),origin='lower')
    #plt.hold(True)
    #plt.scatter(x,100-y,c=z)
    plt.colorbar()

######################################
# interpolate three data arrays using various schemes
#
# call_interpolator(calc,x,y,l,m,n,xi,yi,nx,ny,fault_flag)
# Args:
# calc calculation mode, one of 'simple_idw', 'scipy_idw', 'scipy_rbf'
# l,m,n arrays of direction cosines of pole to plane
# xi,yi arrays of locations of interpolated locations (assumes a grid for plotting, otherwise doesn't matter)
# nx,ny number of x,y elemnts in grid
# fault_flag toggle whether calc for near-fault orientations or not
#
# Call interpolator defined by calc for arrays of arbitrary location x,y located observations as triple or double arrays of 3D or 2D direction cosine arrays (l,m,n) and returns grid of nx ,ny interpolated values for points defined by xi,yi locations. Inspired by https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
######################################
def call_interpolator(calc,x,y,l,m,n,xi,yi,nx,ny,fault_flag):
    # Calculate IDW or other interpolators

      
    if(calc=='simple_idw'):
        ZIl = simple_idw(x,y,l,xi,yi)
    if(calc=='scipy_rbf'):
        ZIl = scipy_rbf(x,y,l,xi,yi)
    if(calc=='scipy_idw'):
        ZIl = scipy_idw(x,y,l,xi,yi)
    if(not fault_flag):
        ZIl = ZIl.reshape((ny, nx))
    
    if(calc=='simple_idw'):
        ZIm = simple_idw(x,y,m,xi,yi)
    if(calc=='scipy_rbf'):
        ZIm = scipy_rbf(x,y,m,xi,yi)
    if(calc=='scipy_idw'):
        ZIm = scipy_idw(x,y,m,xi,yi)
    if(not fault_flag):
        ZIm = ZIm.reshape((ny, nx))
    
    if(type(n) is not int):
        if(calc=='simple_idw'):
            ZIn = simple_idw(x,y,n,xi,yi)
        if(calc=='scipy_rbf'):
            ZIn = scipy_rbf(x,y,n,xi,yi)
        if(calc=='scipy_idw'):
            ZIn = scipy_idw(x,y,n,xi,yi)
        if(not fault_flag):
            ZIn = ZIn.reshape((ny, nx))   
    else:
        ZIn=0
    return(ZIl,ZIm,ZIn)
    
######################################
# Interpolate dipd,dipdirection data from shapefile     
#
# interpolate_orientations(structure_file,tmp_path,bbox,c_l,use_gcode,scheme,gridx,gridy,fault_flag)
# Args:
# structure_file path to orientation layer
# tmp_path directory of temporary outputs from m2l
# bbox bounding box of region of interest
# c_l dictionary of codes and labels specific to input geo information layers
# use_gcode list of groups whose orientation data will be interpolated 
# scheme interpolation scheme one of 'simple_idw', 'scipy_idw', 'scipy_rbf'
# gridx,gridy number of cols & rows in interpolation grid
# fault_flag toggle whether calc for near-fault orientations or not
# 
# Interpolate orientation layer to produce regular grid of l,m,n direction cosines
# Can choose between various RBF and IDW options
# The purpose of these interpolations and associated code is to help in three cases:
# -- Providing estimated dips and contacts in fault-bounded domains where no structural data are available
# -- Needed to estimate true thickness of formations
# -- Useful for poulating parts of maps where little structural data is available
######################################
def interpolate_orientations(structure_file,output_path,bbox,c_l,this_gcode,calc,gridx,gridy,fault_flag):
    structure = gpd.read_file(structure_file,bbox=bbox)
    
    if(len(this_gcode)==1):       
        is_gp=structure[c_l['g']] == thisgcode # subset orientations to just those with this group
        gp_structure = structure[is_gp]
        #print('single group')
        #display(gp_structure)
    else:
        #print('first code',this_gcode[0])
        is_gp=structure[c_l['g']] == this_gcode[0] # subset orientations to just those with this group
        gp_structure = structure[is_gp]
        gp_structure_all = gp_structure.copy()
        #print('first group')
        #display(gp_structure)

        for i in range (1,len(this_gcode)):
            #print('next code',this_gcode[i])
            is_gp=structure[c_l['g']] == this_gcode[i] # subset orientations to just those with this group
            temp_gp_structure = structure[is_gp]
            gp_structure_all = pd.concat([gp_structure_all, temp_gp_structure], ignore_index=True)
            #print('next group')
            #display(gp_structure)

    npts = len(gp_structure_all)
    
    if(fault_flag):
        nx, ny = len(gridx),len(gridy)
    else:
        nx, ny = gridx,gridy
        xi = np.linspace(bbox[0],bbox[2], nx)
        yi = np.linspace(bbox[1],bbox[3], ny)
        xi, yi = np.meshgrid(xi, yi)
        xi, yi = xi.flatten(), yi.flatten()

    x = np.zeros(npts)
    y = np.zeros(npts)
    dip = np.zeros(npts)
    dipdir = np.zeros(npts)
    
    i=0
    for a_pt in gp_structure_all.iterrows():
        x[i]=a_pt[1]['geometry'].x+(np.random.ranf()*0.01)
        y[i]=a_pt[1]['geometry'].y+(np.random.ranf()*0.01)
        dip[i] = a_pt[1][c_l['d']]
        
        if(c_l['otype']=='strike'):
            dipdir[i] = a_pt[1][c_l['dd']]+90
        else:
            dipdir[i] = a_pt[1][c_l['dd']]
        i=i+1


    
    l=np.zeros(npts)
    m=np.zeros(npts)
    n=np.zeros(npts)
    
    for i in range(0,npts):
        l[i],m[i],n[i]=m2l_utils.ddd2dircos(dip[i],dipdir[i])

    if(fault_flag):
        ZIl,ZIm,ZIn=call_interpolator(calc,x,y,l,m,n,gridx,gridy,nx,ny,fault_flag)
    else:
        ZIl,ZIm,ZIn=call_interpolator(calc,x,y,l,m,n,xi,yi,nx,ny,fault_flag)
    
    # Comparisons...
    if(not fault_flag):
        plot(x,-y,l,ZIl)
        plt.title('l')
        plot(x,-y,m,ZIm)
        plt.title('m')
        plot(x,-y,n,ZIn)
        plt.title('n')
    
        plt.show()
    
    if(fault_flag):
        f=open(output_path+'f_input.csv','w')
        fi=open(output_path+'f_interpolation_'+calc+'.csv','w')
        fl=open(output_path+'f_interpolation_l.csv','w')
        fm=open(output_path+'f_interpolation_m.csv','w')
        fn=open(output_path+'f_interpolation_n.csv','w')
    else:
        f=open(output_path+'input.csv','w')
        fi=open(output_path+'interpolation_'+calc+'.csv','w')
        fl=open(output_path+'interpolation_l.csv','w')
        fm=open(output_path+'interpolation_m.csv','w')
        fn=open(output_path+'interpolation_n.csv','w')
    
    f.write("x,y,dip,dipdirection\n")
    fi.write("x,y,dip,dipdirection\n")
    fl.write("x,y,l\n")
    fm.write("x,y,m\n")
    fn.write("x,y,n\n")
    
    for i in range (0,npts):
        ostr="{},{},{},{}\n"\
              .format(x[i],y[i],int(dip[i]),int(dipdir[i]))
        #ostr=str(x[i])+","+str(y[i])+","+str(int(dip[i]))+","+str(int(dipdir[i]))+'\n'
        f.write(ostr)
    
    if(fault_flag):
        for i in range (0,len(gridx)):
                L=ZIl[i]/(sqrt((pow(ZIl[i],2.0))+(pow(ZIm[i],2.0))+(pow(ZIn[i],2.0))))
                M=ZIm[i]/(sqrt((pow(ZIl[i],2.0))+(pow(ZIm[i],2.0))+(pow(ZIn[i],2.0))))
                N=ZIn[i]/(sqrt((pow(ZIl[i],2.0))+(pow(ZIm[i],2.0))+(pow(ZIn[i],2.0))))
                
                dip,dipdir=m2l_utils.dircos2ddd(L,M,N)
                if(isnan(dip) or isnan(dipdir)):
                    dip=dipdir=L=M=N=0
                    print("Warning, no interpolated value for element No. ",i)
                ostr="{},{},{},{}\n"\
                      .format(gridx[i],gridy[i],int(dip),int(dipdir))
                #ostr=str(gridx[i])+","+str(gridy[i])+","+str(int(dip))+","+str(int(dipdir))+'\n'
                fi.write(ostr)
                
                ostr="{},{},{}\n"\
                      .format(gridx[i],gridy[i],L)               
                #ostr=str(gridx[i])+","+str(gridy[i])+","+str(L)+'\n'
                fl.write(ostr)
                ostr="{},{},{}\n"\
                      .format(gridx[i],gridy[i],M)               
                #ostr=str(gridx[i])+","+str(gridy[i])+","+str(M)+'\n'
                fm.write(ostr)
                ostr="{},{},{}\n"\
                      .format(gridx[i],gridy[i],N)               
                #ostr=str(gridx[i])+","+str(gridy[i])+","+str(N)+'\n'
                fn.write(ostr)       
    else:
        for xx in range (0,gridx):
            for yy in range (0,gridy):
                yyy=xx
                xxx=gridy-2-yy
                L=ZIl[xxx,yyy]/(sqrt((pow(ZIl[xxx,yyy],2.0))+(pow(ZIm[xxx,yyy],2.0))+(pow(ZIn[xxx,yyy],2.0))))
                M=ZIm[xxx,yyy]/(sqrt((pow(ZIl[xxx,yyy],2.0))+(pow(ZIm[xxx,yyy],2.0))+(pow(ZIn[xxx,yyy],2.0))))
                N=ZIn[xxx,yyy]/(sqrt((pow(ZIl[xxx,yyy],2.0))+(pow(ZIm[xxx,yyy],2.0))+(pow(ZIn[xxx,yyy],2.0))))
                
                dip,dipdir=m2l_utils.dircos2ddd(L,M,N)
                if(isnan(dip) or isnan(dipdir)):
                    dip=dipdir=L=M=N=0
                    print("Warning, no interpolated value for grid point No. ",xx,',',yy)
                ostr="{},{},{},{}\n"\
                      .format(bbox[0]+(xx*((bbox[2]-bbox[0])/gridx)),bbox[1]+((gridy-1-yy)*((bbox[3]-bbox[1])/gridy)),int(dip),int(dipdir))
                #ostr=str(bbox[0]+(xx*((bbox[2]-bbox[0])/gridx)))+","+str(bbox[1]+((gridy-1-yy)*((bbox[3]-bbox[1])/gridy)))+","+str(int(dip))+","+str(int(dipdir))+'\n'
                fi.write(ostr)
                
                ostr="{},{},{}\n"\
                      .format(xx,yy,L)               
                #ostr=str(xx)+","+str(yy)+","+str(L)+'\n'
                fl.write(ostr)
                ostr="{},{},{}\n"\
                      .format(xx,yy,M)                  
                #ostr=str(xx)+","+str(yy)+","+str(M)+'\n'
                fm.write(ostr)
                ostr="{},{},{}\n"\
                      .format(xx,yy,N)                  
                #ostr=str(xx)+","+str(yy)+","+str(N)+'\n'
                fn.write(ostr)
    
    f.close()
    fi.close()
    fl.close()
    fm.close()
    fn.close()
    
    if(fault_flag):
        print("orientations interpolated as dip dip direction",output_path+'f_interpolation_'+calc+'.csv')
        print("orientations interpolated as l,m,n dir cos",output_path+'f_interpolation_l.csv etc.')
    else:    
        fig, ax = plt.subplots(figsize=(10, 10),)
        q = ax.quiver(xi, yi, -ZIm, ZIl,headwidth=0)
        plt.show()
        print("orientations interpolated as dip dip direction",output_path+'interpolation_'+calc+'.csv')
        print("orientations interpolated as l,m,n dir cos",output_path+'interpolation_l.csv etc.')

######################################
# Interpolate 2D contact data from shapefile
# 
# interpolate_contacts(geology_file,tmp_path,dtm,bbox,c_l,use_gcode,scheme,gridx,gridy,fault_flag)
# Args:
# geology_file path to basal contacts layer
# tmp_path directory of temporary outputs from m2l
# dtm rasterio format elevation grid
# bbox bounding box of region of interest
# c_l dictionary of codes and labels specific to input geo information layers
# use_gcode list of groups whose contact data will be interpolated 
# scheme interpolation scheme one of 'simple_idw', 'scipy_idw', 'scipy_rbf'
# gridx,gridy number of cols & rows in interpolation grid
# fault_flag toggle whether calc for near-fault orientations or not
# 
# Interpolate basal contacts layer to produce regular grid of l,m direction cosines
######################################
def interpolate_contacts(geology_file,output_path,dtm,dtb,dtb_null,cover_map,bbox,c_l,use_gcode,calc,gridx,gridy,fault_flag):
    geol_file = gpd.read_file(geology_file,bbox=bbox)
    #print(len(geol_file))
    #geol_file.plot( color='black',edgecolor='black') 
    
    # Setup: Generate data...
    npts = 0
    decimate=1
    if(fault_flag):
        nx, ny = len(gridx),len(gridy)       
    else:           
        nx, ny= gridx,gridy
        xi = np.linspace(bbox[0],bbox[2], nx)
        yi = np.linspace(bbox[1],bbox[3], ny)
        xi, yi = np.meshgrid(xi, yi)
        xi, yi = xi.flatten(), yi.flatten()
            
    x = np.zeros(20000)   ############## FUDGE ################
    y = np.zeros(20000)   # should go through geology file to see how many contact 
    l = np.zeros(20000)   # segments will be made then define arrays?
    m = np.zeros(20000)   #####################################
    
    if(fault_flag):
        f=open(output_path+'f_raw_contacts.csv','w')
    else:           
        f=open(output_path+'raw_contacts.csv','w')
        
    f.write("X,Y,Z,angle,lsx,lsy,formation,group\n")
    j=0
    i=0
    for indx,acontact in geol_file.iterrows():   #loop through distinct linestrings in MultiLineString
        if(acontact.geometry.type=='MultiLineString'):
            #print(i)
            for line in acontact.geometry: # loop through line segments
                #print(i,len(acontact.geometry))
                if(m2l_utils.mod_safe(i,decimate)  ==0 and acontact[c_l['g']] in use_gcode):
                    #if(acontact['id']==170): 
                        #display(npts,line.coords[0][0],line.coords[1][0]) 
                    dlsx=line.coords[0][0]-line.coords[1][0]
                    dlsy=line.coords[0][1]-line.coords[1][1]
                    if(not line.coords[0][0]==line.coords[1][0] or not line.coords[0][1]==line.coords[1][1]):               
                        lsx=dlsx/sqrt((dlsx*dlsx)+(dlsy*dlsy))
                        lsy=dlsy/sqrt((dlsx*dlsx)+(dlsy*dlsy))
                        x[i]=line.coords[1][0]+(dlsx/2)
                        y[i]=line.coords[1][1]+(dlsy/2)
                        angle=degrees(atan2(lsx,lsy))
                        l[i]=lsx
                        m[i]=lsy
                        locations=[(x[i],y[i])] #doesn't like point right on edge?
                        height=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                        if(str(acontact[c_l['g']])=='None'):
                            ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+"\n"
                        else:
                            ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+","+acontact[c_l['g']].replace(" ","_").replace("-","_")+"\n"
                        f.write(ostr)
                        npts=npts+1
                i=i+1
        else:
            #display(acontact.geometry,acontact.geometry.coords)
            #for line in acontact: # loop through line segments in LineString
            if(  m2l_utils.mod_safe(i,decimate)  ==0 and acontact[c_l['g']] in use_gcode):
                dlsx=acontact.geometry.coords[0][0]-acontact.geometry.coords[1][0]
                dlsy=acontact.geometry.coords[0][1]-acontact.geometry.coords[1][1]
                if(not acontact.geometry.coords[0][0]==acontact.geometry.coords[1][0] 
                   or not acontact.geometry.coords[0][1]==acontact.geometry.coords[1][1]):
                    lsx=dlsx/sqrt((dlsx*dlsx)+(dlsy*dlsy))
                    lsy=dlsy/sqrt((dlsx*dlsx)+(dlsy*dlsy))
                    x[i]=acontact.geometry.coords[1][0]+(dlsx/2)
                    y[i]=acontact.geometry.coords[1][1]+(dlsy/2)
                    angle=degrees(atan2(lsx,lsy))
                    l[i]=lsx
                    m[i]=lsy
                    locations=[(x[i],y[i])] #doesn't like point right on edge?
                    height=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                    if(str(acontact[c_l['g']])=='None'):
                        ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+"\n"
                    else:
                        ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+","+acontact[c_l['g']].replace(" ","_").replace("-","_")+"\n"
                    #print(ostr)
                    f.write(ostr)
                    #print(npts,dlsx,dlsy)
                    npts=npts+1
                i=i+1
        j=j+1
    f.close()
    #print("i",i,"npts",npts)

    for i in range(0,npts):
        x[i]=x[i]+(np.random.ranf()*0.01)
        y[i]=y[i]+(np.random.ranf()*0.01)

    if(fault_flag):
        ZIl,ZIm,ZIn=call_interpolator(calc,x[:npts],y[:npts],l[:npts],m[:npts],0,gridx,gridy,nx,ny,fault_flag)        
    else:
        ZIl,ZIm,ZIn=call_interpolator(calc,x[:npts],y[:npts],l[:npts],m[:npts],0,xi,yi,nx,ny,fault_flag)        
    
    # Comparisons...
    if(not fault_flag):
        plot(x,-y,l,ZIl)
        plt.title('l')
        plot(x,-y,m,ZIm)
        plt.title('m')

    if(fault_flag):
        fi=open(output_path+'f_interpolation_contacts_'+calc+'.csv','w')
        fl=open(output_path+'f_interpolation_contacts_l.csv','w')
        fm=open(output_path+'f_interpolation_contacts_m.csv','w')       
    else:
        fi=open(output_path+'interpolation_contacts_'+calc+'.csv','w')
        fl=open(output_path+'interpolation_contacts_l.csv','w')
        fm=open(output_path+'interpolation_contacts_m.csv','w')
    
    fi.write("x,y,angle\n")
    fl.write("x,y,l\n")
    fm.write("x,y,m\n")
    
    
    if(fault_flag):
        for i in range (0,len(gridx)):
                L=ZIl[i]/(sqrt((pow(ZIl[i],2.0))+(pow(ZIm[i],2.0))))
                M=ZIm[i]/(sqrt((pow(ZIl[i],2.0))+(pow(ZIm[i],2.0))))
                S=degrees(atan2(L,M))
                
                if(isnan(S)):
                    S=0
                    print("Warning, no interpolated value for element No. ",i)
                
                ostr=str(gridx[i])+","+str(gridy[i])+","+str(int(S))+'\n'
                fi.write(ostr)
                
                ostr=str(gridx[i])+","+str(gridy[i])+","+str(L)+'\n'
                fl.write(ostr)
                ostr=str(gridx[i])+","+str(gridy[i])+","+str(M)+'\n'
                fm.write(ostr)
    else:
        for xx in range (0,gridx):
            for yy in range (0,gridy):
                yyy=xx
                xxx=gridy-2-yy
                L=ZIl[xxx,yyy]/(sqrt((pow(ZIl[xxx,yyy],2.0))+(pow(ZIm[xxx,yyy],2.0))))
                M=ZIm[xxx,yyy]/(sqrt((pow(ZIl[xxx,yyy],2.0))+(pow(ZIm[xxx,yyy],2.0))))
                S=degrees(atan2(L,M))
        
                ostr=str(bbox[0]+(xx*((bbox[2]-bbox[0])/(gridx))))+","+str(bbox[1]+((gridy-2-yy)*((bbox[3]-bbox[1])/(gridy))))+","+str(int(S))+'\n'
                fi.write(ostr)
                
                ostr=str(xx)+","+str(yy)+","+str(L)+'\n'
                fl.write(ostr)
                ostr=str(xx)+","+str(yy)+","+str(M)+'\n'
                fm.write(ostr)


    fi.close()
    fl.close()
    fm.close()
    if(fault_flag):
        print("contacts interpolated as strike",output_path+'f_interpolation_contacts_'+calc+'.csv')
        print("contacts interpolated as l,m dir cos",output_path+'f_interpolation_contacts_l.csv etc.')
    else:
        fig, ax = plt.subplots(figsize=(10, 10))
        q = ax.quiver(xi, yi, ZIl, ZIm,headwidth=0)
        plt.show()
        print("contacts interpolated as strike",output_path+'interpolation_contacts_'+calc+'.csv')
        print("contacts interpolated as l,m dir cos",output_path+'interpolation_contacts_l.csv etc.')

######################################
# save all contacts as vectors (used for debugging)
#
# save_contact_vectors(geology_file,tmp_path,dtm,bbox,c_l,calc,decimate)
# Args:
# geology_file file path to geology polygons
# tmp_path directory of temporary outputs from m2l
# dtm rasterio format dtm raster
# bbox bounding box of model
# c_l dictionary of codes and labels specific to input geo information layers
# calc NOT USED
# decimate simple decimation factor for saving vectors
######################################
def save_contact_vectors(geology_file,tmp_path,dtm,dtb,dtb_null,cover_map,bbox,c_l,calc,decimate):
    geol_file = gpd.read_file(geology_file,bbox=bbox)
    print(len(geol_file))
    #geol_file.plot( color='black',edgecolor='black') 
    
    npts = 0
    i=0
    for indx,acontact in geol_file.iterrows():   #loop through distinct linestrings in MultiLineString
        if(acontact.geometry.type=='MultiLineString'):
            for line in acontact.geometry: # loop through line segments
                if(m2l_utils.mod_safe(i,decimate)  ==0):
                    npoint=1
                i=i+1
        else:
            if(  m2l_utils.mod_safe(i,decimate)  ==0):
                npoint=1
            i=i+1

    x = np.zeros(i+1)
    y = np.zeros(i+1)
    l = np.zeros(i+1)
    m = np.zeros(i+1)
    
    f=open(tmp_path+'raw_contacts.csv','w')
    f.write("X,Y,Z,angle,lsx,lsy,formation,group\n")
    j=0
    i=0
    for indx,acontact in geol_file.iterrows():   #loop through distinct linestrings in MultiLineString
        if(acontact.geometry.type=='MultiLineString'):
            #print(i)
            for line in acontact.geometry: # loop through line segments
                #print(i,len(acontact.geometry))
                if(m2l_utils.mod_safe(i,decimate)  ==0):
                    #if(acontact['id']==170): 
                        #display(npts,line.coords[0][0],line.coords[1][0]) 
                    dlsx=line.coords[0][0]-line.coords[1][0]
                    dlsy=line.coords[0][1]-line.coords[1][1]
                    if(not line.coords[0][0]==line.coords[1][0] or not line.coords[0][1]==line.coords[1][1]):               
                        lsx=dlsx/sqrt((dlsx*dlsx)+(dlsy*dlsy))
                        lsy=dlsy/sqrt((dlsx*dlsx)+(dlsy*dlsy))
                        x[i]=line.coords[1][0]+(dlsx/2)
                        y[i]=line.coords[1][1]+(dlsy/2)
                        angle=degrees(atan2(lsx,lsy))
                        l[i]=lsx
                        m[i]=lsy
                        locations=[(x[i],y[i])] #doesn't like point right on edge?
                        height=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)

                        if(str(acontact[c_l['g']])=='None'):
                            ostr="{},{},{},{},{},{},{},{}\n"\
                                  .format(x[i],y[i],height,angle%180,lsx,lsy,acontact[c_l['c']].replace(" ","_").replace("-","_"),acontact[c_l['c']].replace(" ","_").replace("-","_"))
                            #ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+"\n"
                        else:
                            ostr="{},{},{},{},{},{},{},{}\n"\
                                  .format(x[i],y[i],height,angle%180,lsx,lsy,acontact[c_l['c']].replace(" ","_").replace("-","_"),acontact[c_l['g']].replace(" ","_").replace("-","_"))
                            #ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+","+acontact[c_l['g']].replace(" ","_").replace("-","_")+"\n"
                        f.write(ostr)
                        npts=npts+1
                i=i+1
        else:
            #display(acontact.geometry,acontact.geometry.coords)
            #for line in acontact: # loop through line segments in LineString
            if(  m2l_utils.mod_safe(i,decimate)  ==0):
                dlsx=acontact.geometry.coords[0][0]-acontact.geometry.coords[1][0]
                dlsy=acontact.geometry.coords[0][1]-acontact.geometry.coords[1][1]
                if(not acontact.geometry.coords[0][0]==acontact.geometry.coords[1][0] 
                   or not acontact.geometry.coords[0][1]==acontact.geometry.coords[1][1]):
                    lsx=dlsx/sqrt((dlsx*dlsx)+(dlsy*dlsy))
                    lsy=dlsy/sqrt((dlsx*dlsx)+(dlsy*dlsy))
                    x[i]=acontact.geometry.coords[1][0]+(dlsx/2)
                    y[i]=acontact.geometry.coords[1][1]+(dlsy/2)
                    angle=degrees(atan2(lsx,lsy))
                    l[i]=lsx
                    m[i]=lsy
                    locations=[(x[i],y[i])] #doesn't like point right on edge?
                    height=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                    if(str(acontact[c_l['g']])=='None'):
                        ostr="{},{},{},{},{},{},{},{}\n"\
                                  .format(x[i],y[i],height,angle%180,lsx,lsy,acontact[c_l['c']].replace(" ","_").replace("-","_"),acontact[c_l['c']].replace(" ","_").replace("-","_"))
                        #ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+"\n"
                    else:
                        ostr="{},{},{},{},{},{},{},{}\n"\
                                  .format(x[i],y[i],height,angle%180,lsx,lsy,acontact[c_l['c']].replace(" ","_").replace("-","_"),acontact[c_l['g']].replace(" ","_").replace("-","_"))
                        #ostr=str(x[i])+","+str(y[i])+","+str(height)+","+str(angle%180)+","+str(lsx)+","+str(lsy)+","+acontact[c_l['c']].replace(" ","_").replace("-","_")+","+acontact[c_l['g']].replace(" ","_").replace("-","_")+"\n"
                    #print(ostr)
                    f.write(ostr)
                    #print(npts,dlsx,dlsy)
                    npts=npts+1
                i=i+1
        j=j+1
    f.close()
    print(npts,'points saved to',tmp_path+'raw_contacts.csv')


####################################
# combine interpolated contact information (to provide l,m with interpolated dip,dipdirection data (to provide n) 
#
# join_contacts_and_orientations(combo_file,geology_file,tmp_path,dtm_reproj_file,c_l,lo,mo,no,lc,mc,xy,dst_crs,bbox,fault_flag)
# combo_file path to temporary combined information geology_file path to basal contacts layer
# tmp_path directory of temporary outputs from m2l
# dtm_reproj_file path to reprojected dtm file
# c_l dictionary of codes and labels specific to input geo information layers
# lo,mo,no 3D direction cosines of interpolated orientations
# lc,mc 2D direction cosines of interpolated contacts
# xy interpolated orientations (used to get x,y locations only) dst_crs Coordinate Reference System of destination geotif (any length-based projection)
# bbox bounding box of region of interest
# fault_flag toggle whether calc for near-fault orientations or not
# 
# Combine interpolation orientations with interpolated basal contacts layers to produce regular grid of interpolated dip, dip direction estimates
# Uses normalised direction cosines (l,m,n):
# -- l,m from RBF of basal contact orientations -- signs of l & m from misorientation with RBF of orientation data and -- n from RBF of orientation data
# 
# Useful for adding data where no orientations are available (e.g. in fault bounded domains) and for calculating true thickness of layers. Assumes a 2D plane of data, but if 3D RBF was calulated and projected contact info was used it should apply with topography too.
####################################    
def join_contacts_and_orientations(combo_file,geology_file,output_path,dtm_reproj_file,dtb,dtb_null,cover_map,c_l,lo,mo,no,lc,mc,xy,dst_crs,bbox,fault_flag):
    f=open(combo_file,'w')
    f.write('x,y,dip,dipdirection,misorientation,dotproduct\n')

    for i in range(0,len(lc)):
        #print(mc[i,2],lc[i,2],lo[i,2],mo[i,2],no[i,2])
        scale=sqrt(1-pow(no[i,2],2)) #scaling contact dircos to *include* dip info
        lcscaled=scale*-mc[i,2] #includes 90 rotation to account for orthogonality of contact and dip direction
        mcscaled=scale*lc[i,2]
        scale2=sqrt(pow(lo[i,2],2)+pow(mo[i,2],2)) #scaling dip dipdir dircos to *exclude* dip info
        if(scale2>0.0):
            loscaled=lo[i,2]/scale2
            moscaled=mo[i,2]/scale2
        else:
            loscaled=0
            moscaled=0
        dotproduct=(-mc[i,2]*loscaled)+(lc[i,2]*moscaled) #includes 90 rotation to account for orthogonality of contact and dip direction
        if(dotproduct<0):
            lcscaled=-lcscaled
            mcscaled=-mcscaled

        misorientation=degrees(acos(dotproduct))
        dip,dipdir=m2l_utils.dircos2ddd(lcscaled,mcscaled,no[i,2])
        ostr="{},{},{},{},{},{}\n"\
              .format(xy[i,0],xy[i,1],int(dip),int(dipdir),int(misorientation),dotproduct)
        #ostr=str(xy[i,0])+','+str(xy[i,1])+','+str(int(dip))+','+str(int(dipdir))+','+str(int(misorientation))+','+str(dotproduct)+'\n'
        f.write(ostr)
    f.close()   
    
    
    geology = gpd.read_file(geology_file,bbox=bbox)
    geology.crs=dst_crs
    geology = m2l_utils.explode(geology)
    
    data = pd.read_csv(combo_file)
    
    geometry = [Point(xy) for xy in zip(data['x'], data['y'])]
        
    gdf = GeoDataFrame(data, crs=dst_crs, geometry=geometry)
    
    gdf.crs=dst_crs
    print(gdf.crs,geology.crs)    
    structure_code = gpd.sjoin(gdf, geology, how="left", op="within")
    dtm = rasterio.open(dtm_reproj_file)
    if(fault_flag):
        f=open(output_path+'f_combo_full.csv','w')
    else:
        f=open(output_path+'combo_full.csv','w')
    f.write('X,Y,Z,azimuth,dip,polarity,formation\n')
    last_code=''
    for indx,a_point in structure_code.iterrows(): 
        
        locations=[(a_point['x'],a_point['y'])]
        height=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
        ostr=str(a_point['x'])+','
        ostr=ostr+str(a_point['y'])+','
        ostr=ostr+str(height)+','+str(int(a_point['dipdirection']))+','
        ostr=ostr+str(int(a_point['dip']))+',1,'
        ostr=ostr+str(a_point[c_l['c']]).replace("-","_").replace(" ","_")+'\n'
        
        if(not str(a_point[c_l['c']])=='nan' ):
            f.write(ostr)
        last_code=a_point[c_l['c']]
    f.close()  
    if(fault_flag):
        print("contacts and orientations interpolated as dip dip direction",output_path+'f_combo_full.csv')
    else:
        print("contacts and orientations interpolated as dip dip direction",output_path+'combo_full.csv')
 
######################################
# Interpolate dipd,dipdirection data from shapefile usin fold axial traces as additional constraints     
# interpolate_orientations_with_fat(structure_file,output_path,bbox,c_l,this_gcode,calc,gridx,gridy)
# structure_file path to orientation layer
# output_path directory for outputs from m2l
# bbox bounding box of region of interest
# c_l dictionary of codes and labels specific to input geo information layers
# this_gcode list of groups whose orientation data will be interpolated 
# calc interpolation scheme one of 'simple_idw', 'scipy_idw', 'scipy_rbf'
# gridx,gridy number of cols & rows in interpolation grid
# 
# Interpolate orientation layer to produce regular grid of l,m,n direction cosines
# Can choose between various RBF and IDW options
# The purpose of these interpolations and associated code is to help in three cases:
# -- Providing estimated dips and contacts in fault-bounded domains where no structural data are available
# -- Needed to estimate true thickness of formations
# -- Useful for poulating parts of maps where little structural data is available
######################################
def interpolate_orientations_with_fat(structure_file,output_path,bbox,c_l,this_gcode,calc,gridx,gridy):
    structure = gpd.read_file(structure_file,bbox=bbox)
    fat_orientations=pd.read_csv(output_path+'fold_axial_trace_orientations2.csv',",")
    
    
    
    if(len(this_gcode)==1):       
        is_gp=structure[c_l['g']] == thisgcode # subset orientations to just those with this group
        gp_structure = structure[is_gp]
        print('single group')
        display(gp_structure)
    else:
        print('first code',this_gcode[0])
        is_gp=structure[c_l['g']] == this_gcode[0] # subset orientations to just those with this group
        gp_structure = structure[is_gp]
        gp_structure_all = gp_structure.copy()
        print('first group')
        display(gp_structure)

        for i in range (1,len(this_gcode)):
            print('next code',this_gcode[i])
            is_gp=structure[c_l['g']] == this_gcode[i] # subset orientations to just those with this group
            temp_gp_structure = structure[is_gp]
            gp_structure_all = pd.concat([gp_structure_all, temp_gp_structure], ignore_index=True)
            print('next group')
            display(gp_structure)

    npts = len(gp_structure_all)+len(fat_orientations)
    
    nx, ny = gridx,gridy

    xi = np.linspace(bbox[0],bbox[2], nx)
    yi = np.linspace(bbox[1],bbox[3], ny)
    xi, yi = np.meshgrid(xi, yi)
    xi, yi = xi.flatten(), yi.flatten()
    x = np.zeros(npts)
    y = np.zeros(npts)
    dip = np.zeros(npts)
    dipdir = np.zeros(npts)
    
    i=0
    for ind,a_pt in gp_structure_all.iterrows():
        x[i]=a_pt['geometry'].x
        y[i]=a_pt['geometry'].y
        dip[i] = a_pt[c_l['d']]
        if(c_l['otype']=='strike'):
            dipdir[i] = float(a_pt[c_l['dd']])+90
        else:
            dipdir[i] = a_pt[c_l['dd']]
        i=i+1

    for ind,a_pt in fat_orientations.iterrows():
        x[i]=a_pt['X']
        y[i]=a_pt['Y']
        dip[i] = a_pt['dip']
        dipdir[i] = a_pt['azimuth']
        i=i+1
    
    l=np.zeros(npts)
    m=np.zeros(npts)
    n=np.zeros(npts)
    
    for i in range(0,npts):
        l[i],m[i],n[i]=m2l_utils.ddd2dircos(dip[i],dipdir[i])

    ZIl,ZIm,ZIn=call_interpolator(calc,x,y,l,m,n,xi,yi,nx,ny)
    
    # Comparisons...
    plot(x,-y,l,ZIl)
    plt.title('l')
    plot(x,-y,m,ZIm)
    plt.title('m')
    plot(x,-y,n,ZIn)
    plt.title('n')
    
    plt.show()
    
    f=open(output_path+'input.csv','w')
    fi=open(output_path+'interpolation_'+calc+'.csv','w')
    fl=open(output_path+'interpolation_l.csv','w')
    fm=open(output_path+'interpolation_m.csv','w')
    fn=open(output_path+'interpolation_n.csv','w')
    
    f.write("x,y,dip,dipdirection\n")
    fi.write("x,y,dip,dipdirection\n")
    fl.write("x,y,l\n")
    fm.write("x,y,m\n")
    fn.write("x,y,n\n")
    
    for i in range (0,npts):
        ostr="{},{},{},{}\n"\
              .format(x[i],y[i],int(dip[i]),int(dipdir[i]))
        #ostr=str(x[i])+","+str(y[i])+","+str(int(dip[i]))+","+str(int(dipdir[i]))+'\n'
        f.write(ostr)
    
    for xx in range (0,gridx):
        for yy in range (0,gridy):
            yyy=xx
            xxx=gridy-2-yy
            L=ZIl[xxx,yyy]/(sqrt((pow(ZIl[xxx,yyy],2.0))+(pow(ZIm[xxx,yyy],2.0))+(pow(ZIn[xxx,yyy],2.0))))
            M=ZIm[xxx,yyy]/(sqrt((pow(ZIl[xxx,yyy],2.0))+(pow(ZIm[xxx,yyy],2.0))+(pow(ZIn[xxx,yyy],2.0))))
            N=ZIn[xxx,yyy]/(sqrt((pow(ZIl[xxx,yyy],2.0))+(pow(ZIm[xxx,yyy],2.0))+(pow(ZIn[xxx,yyy],2.0))))
            
            dip,dipdir=m2l_utils.dircos2ddd(L,M,N)
            ostr="{},{},{},{}\n"\
                  .format(bbox[0]+(xx*((bbox[2]-bbox[0])/gridx)),bbox[1]+((gridy-1-yy)*((bbox[3]-bbox[1])/gridy)),int(dip),int(dipdir))
            #ostr=str(bbox[0]+(xx*((bbox[2]-bbox[0])/gridx)))+","+str(bbox[1]+((gridy-1-yy)*((bbox[3]-bbox[1])/gridy)))+","+str(int(dip))+","+str(int(dipdir))+'\n'
            fi.write(ostr)
            
            ostr="{},{},{}\n".format(xx,yy,L)
            #ostr=str(xx)+","+str(yy)+","+str(L)+'\n'
            fl.write(ostr)
            ostr="{},{},{}\n".format(xx,yy,M)
            #ostr=str(xx)+","+str(yy)+","+str(M)+'\n'
            fm.write(ostr)
            ostr="{},{},{}\n".format(xx,yy,N)
            #ostr=str(xx)+","+str(yy)+","+str(N)+'\n'
            fn.write(ostr)
    
    f.close()
    fi.close()
    fl.close()
    fm.close()
    fn.close()
    
    fig, ax = plt.subplots(figsize=(10, 10),)
    q = ax.quiver(xi, yi, -ZIm, ZIl,headwidth=0)
    plt.show()
    print("orientations interpolated as dip dip direction",output_path+'interpolation_'+calc+'.csv')
    print("orientations interpolated as l,m,n dir cos",output_path+'interpolation_l.csv etc.')

####################################################
# For each fault string:
# process_fault_throw_and_near_orientations(tmp_path,output_path,dtm_reproj_file,c_l,use_gcode,use_gcode2,dst_crs,bbox,scheme)
# Args:
#
#    incementally advance along polyline every at each inter-node (no point in doing more?)
#    find local stratigraphy 10m to left and right of fault
# Once full fault has been traversed:
#
#    Find list of contacts left
#    Find equivalent contacts on right
#    use interpolated orientations to estimate minimum true offset assuming vertical displacement and store
#    if no equivalent found, flag as domain fault and find min strat offset for contact, use cumulative minimum thickness estimate and store with flag (not implemented)
#    estimate median & sd of minimum fault offset and store with flag (not implemented)
# Local Orientations: Since much of the code is the same, we benefit by calculating local orientation data either side of 
# fault so that geomodeller/gempy have satisfied fault compartment orientation data
###################################################

def process_fault_throw_and_near_orientations(tmp_path,output_path,dtm_reproj_file,dtb,dtb_null,cover_map,c_l,use_gcode,use_gcode2,dst_crs,bbox,scheme):
    fault_file=tmp_path+'faults_clip.shp'
    geology_file=tmp_path+'geol_clip.shp'

    faults = gpd.read_file(fault_file)
    geology = gpd.read_file(geology_file)
    dtm = rasterio.open(dtm_reproj_file)

    all_long_faults=np.genfromtxt(output_path+'fault_dimensions.csv',delimiter=',',dtype='U100')
    fault_names=all_long_faults[1:,:1]
    m_step=10.0 #outstep from fault 
    # if we have two faults at low angle this can be a problem as the point crosses the next fault
    xi=[]
    yi=[]
    fdc=[]
    all_coordsdist=[]
    all_coords_x=[]
    all_coords_y=[]

    fftc=open(output_path+'fault_tip_contacts.csv','w')
    fftc.write('X,Y,Z,formation\n')
    
    # loop through all faults
    
    for index,fault in faults.iterrows():
        #if(fault[c_l['o']]!=1071):
            #continue
        if(fault.geometry.type=='LineString'):
            if('Fault_'+str(fault[c_l['o']]) in fault_names): # in is dangerous as Fault_1 is in Fault_10
                #print('LineString','Fault_'+str(fault[c_l['o']]),len(fault.geometry.coords))
                lcoords=[]
                rcoords=[]
                index=[]
                
                # make a list of points just offset from mid-points between fault nodes and join
                # geology polygon information to points
                j=0                
                for i in range (0,len(fault.geometry.coords)-1):
                    for inc in np.arange(0.01, 1, 0.01): 
                        midx=fault.geometry.coords[i][0]+((fault.geometry.coords[i+1][0]-fault.geometry.coords[i][0])*inc)            
                        midy=fault.geometry.coords[i][1]+((fault.geometry.coords[i+1][1]-fault.geometry.coords[i][1])*inc)
                        l,m=m2l_utils.pts2dircos(fault.geometry.coords[i][0],fault.geometry.coords[i][1],fault.geometry.coords[i+1][0],fault.geometry.coords[i+1][1])
                        lcoords.append([(midx+(m_step*m),midy-(m_step*l))])
                        rcoords.append([(midx-(m_step*m),midy+(m_step*l))])
                        index.append([(j)])
                        j=j+1
                lgeom=[Point(xy) for xy in lcoords]        
                rgeom=[Point(xy) for xy in rcoords]

                lgdf = GeoDataFrame(index, crs=dst_crs, geometry=lgeom)
                rgdf = GeoDataFrame(index, crs=dst_crs, geometry=rgeom)
                lcode = gpd.sjoin(lgdf, geology, how="left", op="within")        
                rcode = gpd.sjoin(rgdf, geology, how="left", op="within")
                # display(lcode)
                
                # add points to list if they have different geology code than previous node on left side
                
                first=True
                lcontact=[]
                lastlcode=''
                
                for ind,indl in lcode.iterrows():
                    
                    if(ind<len(lcode) and not isnan(indl['index_right'])):
                        ntest1=str(indl[c_l['ds']])
                        ntest2=str(indl[c_l['r1']])
                        
                        if(not ntest1 == 'None' and not ntest2 == 'None'):
                            if(ind==1 or (not lastlcode==indl[c_l['c']] and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']]) ))):
                                all_coords_x.append(indl.geometry.x)
                                all_coords_y.append(indl.geometry.y)
                                if(first):
                                    first=False
                                    firstlx=indl.geometry.x
                                    firstly=indl.geometry.y
                                    firstlc=indl[c_l['c']].replace(" ","_").replace("-","_")
                                lastlx=indl.geometry.x
                                lastly=indl.geometry.y
                                lastlc=indl[c_l['c']].replace(" ","_").replace("-","_")
                                
    
                            if(lastlcode=='' and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']] ))):
                                lastlcode=indl[c_l['c']]
    
                            #print('l',ind,indl[c_l['c']],indl[c_l['ds']],indl[c_l['r1']])
                            if(not ntest1 == 'None' and not ntest2 == 'None' and not str(indl[c_l['c']])=='nan'):
                                if((not indl[c_l['c']]==lastlcode) and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']] ))):
                                    lcontact.append([(ind,lastlcode,indl[c_l['c']],indl.geometry)])
                                    lastlcode=indl[c_l['c']]
                    
                #add points to list if they have different geology code than previous node on right side
                
                first=True
                rcontact=[]
                lastrcode=''
                for ind,indr in rcode.iterrows():
                    if(ind<len(rcode) and not isnan(indr['index_right'])):
                        ntest1=str(indr[c_l['ds']])
                        ntest2=str(indr[c_l['r1']])
                        
                        if(not ntest1 == 'None' and not ntest2 == 'None'):                            
                            if(ind==1 or (not lastrcode==indr[c_l['c']] and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']]) ))):
                                all_coords_x.append(indr.geometry.x)
                                all_coords_y.append(indr.geometry.y)
                                if(first):
                                    first=False
                                    firstrx=indr.geometry.x
                                    firstry=indr.geometry.y
                                    firstrc=indr[c_l['c']].replace(" ","_").replace("-","_")
                                lastrx=indr.geometry.x
                                lastry=indr.geometry.y
                                lastrc=indr[c_l['c']].replace(" ","_").replace("-","_")
    
                            if(lastrcode=='' and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']] ))):
                                lastrcode=indr[c_l['c']]                        
                            #print('r',ind,indr[c_l['c']],indr[c_l['ds']],indr[c_l['r1']])
    
                            #print(lastrcode,ntest1,ntest2,str(indr[c_l['c']]),indr[c_l['c']],c_l['sill'],c_l['intrusive'])
                            
                            if(not ntest1 == 'None' and not ntest2 == 'None' and not str(indr[c_l['c']])=='nan' ):
                                if((not indr[c_l['c']]==lastrcode) and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']] ))):
                                    rcontact.append([(ind,lastrcode,indr[c_l['c']],indr.geometry)]) 
                                    lastrcode=indr[c_l['c']]             

                locations=[(firstlx,firstly)]
                first_height_l=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                ostr="{},{},{},{}\n"\
                              .format(firstlx,firstly,first_height_l,firstlc)
                fftc.write(ostr)
                locations=[(firstrx,firstry)]
                first_height_r=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                ostr="{},{},{},{}\n"\
                              .format(firstrx,firstry,first_height_r,firstrc)
                fftc.write(ostr)
                locations=[(lastlx,lastly)]
                last_height_l=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                ostr="{},{},{},{}\n"\
                              .format(lastlx,lastly,last_height_l,lastlc)
                fftc.write(ostr)
                locations=[(lastrx,lastry)]
                last_height_r=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                ostr="{},{},{},{}\n"\
                              .format(lastrx,lastry,last_height_r,lastrc)
                fftc.write(ostr)   
                                    
                # loop through left and right sides to find equivalent contact pairs along fault
                
                if(len(lcontact)>0 and len(rcontact)>0):                
                    for lc in lcontact:
                        for rc in rcontact:
                            #display('l',lc[0][3].x,'r',rc[0][3].x)
                            if(lc[0][1]==rc[0][1] and lc[0][2]==rc[0][2] and not lc[0][1]==''):
                                dist=m2l_utils.ptsdist(lc[0][3].x,lc[0][3].y,rc[0][3].x,rc[0][3].y)
                                if(lc[0][0]<rc[0][0]):
                                    dist=-dist
                                #print('***',lc,rc)

                                xi.append((lc[0][3].x))
                                yi.append((lc[0][3].y))
                                l,m=m2l_utils.pts2dircos(lc[0][3].x,lc[0][3].y,rc[0][3].x,rc[0][3].y)
                                if(not (l==0.0 and m==0.0)):
                                    fdc.append((l,m,'Fault_'+str(fault[c_l['o']])))
                                    all_coordsdist.append((dist))
        else:
            if('Fault_'+str(fault[c_l['o']]) in fault_names): # in is dangerous as Fault_1 is in Fault_10
                for fls in fault.geometry:              
                    fault_ls=LineString(fls)        
                    lcoords=[]
                    rcoords=[]
                    index=[]
                    #display("MLS DEBUG",fault.geometry.type)
                    
                    j=0                
                    for i in range (0,len(fault_ls.coords)-1):
                        for inc in np.arange(0.01, 1, 0.01): 
                            midx=fault_ls.coords[i][0]+((fault_ls.coords[i+1][0]-fault_ls.coords[i][0])*inc)            
                            midy=fault_ls.coords[i][1]+((fault_ls.coords[i+1][1]-fault_ls.coords[i][1])*inc)
                            l,m=m2l_utils.pts2dircos(fault_ls.coords[i][0],fault_ls.coords[i][1],fault_ls.coords[i+1][0],fault_ls.coords[i+1][1])
                            lcoords.append([(midx+(m_step*m),midy-(m_step*l))])
                            rcoords.append([(midx-(m_step*m),midy+(m_step*l))])
                            index.append([(j)])
                            j=j+1

                    lgeom=[Point(xy) for xy in lcoords]        
                    rgeom=[Point(xy) for xy in rcoords]
                    lgdf = GeoDataFrame(index, crs=dst_crs, geometry=lgeom)
                    rgdf = GeoDataFrame(index, crs=dst_crs, geometry=rgeom)
                    lcode = gpd.sjoin(lgdf, geology, how="left", op="within")        
                    rcode = gpd.sjoin(rgdf, geology, how="left", op="within")

                    #add points to list if they have different geology code than previous node on left side

                    first=True
                    lcontact=[]
                    lastlcode=''
                    for ind,indl in lcode.iterrows():

                        if(ind<len(lcode) and not isnan(indl['index_right'])):
                            ntest1=str(indl[c_l['ds']])
                            ntest2=str(indl[c_l['r1']])
                            
                            if(not ntest1 == 'None' and not ntest2 == 'None'):
                                if(ind==1 or (not lastlcode==indl[c_l['c']] and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']]) ))):
                                    all_coords_x.append(indl.geometry.x)
                                    all_coords_y.append(indl.geometry.y)
                                    if(first):
                                        first=False
                                        firstlx=indl.geometry.x
                                        firstly=indl.geometry.y
                                        firstlc=indl[c_l['c']].replace(" ","_").replace("-","_")
                                    lastlx=indl.geometry.x
                                    lastly=indl.geometry.y
                                    lastlc=indl[c_l['c']].replace(" ","_").replace("-","_")
    
                                if(lastlcode=='' and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']] ))):
                                    lastlcode=indl[c_l['c']]
    
                                if(not ntest1 == 'None' and not ntest2 == 'None' and not str(indl[c_l['c']])=='nan'):
                                    if((not indl[c_l['c']]==lastlcode) and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']] ))):
                                        lcontact.append([(ind,lastlcode,indl[c_l['c']],indl.geometry)])
                                        lastlcode=indl[c_l['c']]

                    #add points to list if they have different geology code than previous node on right side

                    first=True
                    rcontact=[]
                    lastrcode=''
                    for ind,indr in rcode.iterrows():
                        if(ind<len(rcode) and not isnan(indr['index_right'])):
                            ntest1=str(indr[c_l['ds']])
                            ntest2=str(indr[c_l['r1']])
                            
                            if(not ntest1 == 'None' and not ntest2 == 'None'):
                                if(ind==1 or (not lastrcode==indr[c_l['c']] and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']]) ))):
                                    all_coords_x.append(indr.geometry.x)
                                    all_coords_y.append(indr.geometry.y)
                                    if(first):
                                        first=False
                                        firstrx=indr.geometry.x
                                        firstry=indr.geometry.y
                                        firstrc=indr[c_l['c']].replace(" ","_").replace("-","_")
                                    lastrx=indr.geometry.x
                                    lastry=indr.geometry.y
                                    lastrc=indr[c_l['c']].replace(" ","_").replace("-","_")
    
                                if(lastrcode=='' and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']] ))):
                                    lastrcode=indr[c_l['c']]                        
    
                                if(not ntest1 == 'None' and not ntest2 == 'None' and not str(indr[c_l['c']])=='nan' ):
                                    if((not indr[c_l['c']]==lastrcode) and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']] ))):
                                        rcontact.append([(ind,lastrcode,indr[c_l['c']],indr.geometry)]) 
                                        lastrcode=indr[c_l['c']]

                    locations=[(firstlx,firstly)]
                    first_height_l=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                    ostr="{},{},{},{}\n"\
                                  .format(firstlx,firstly,first_height_l,firstlc)
                    fftc.write(ostr)
                    locations=[(firstrx,firstry)]
                    first_height_r=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                    ostr="{},{},{},{}\n"\
                                  .format(firstrx,firstry,first_height_r,firstrc)
                    fftc.write(ostr)
                    locations=[(lastlx,lastly)]
                    last_height_l=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                    ostr="{},{},{},{}\n"\
                                  .format(lastlx,lastly,last_height_l,lastlc)
                    fftc.write(ostr)
                    locations=[(lastrx,lastry)]
                    last_height_r=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                    ostr="{},{},{},{}\n"\
                                  .format(lastrx,lastry,last_height_r,lastrc)
                    fftc.write(ostr)   
        
                    # loop through left and right sides to find equivalent contact pairs along fault

                    if(len(lcontact)>0 and len(rcontact)>0):
                        for lc in lcontact:
                            for rc in rcontact:
                                #display('l',lc[0][3].x,'r',rc[0][3].x)
                                if(lc[0][1]==rc[0][1] and lc[0][2]==rc[0][2] and not lc[0][1]==''):
                                    dist=m2l_utils.ptsdist(lc[0][3].x,lc[0][3].y,rc[0][3].x,rc[0][3].y)
                                    if(lc[0][0]<rc[0][0]):
                                        dist=-dist
                                    #print('***',lc,rc)

                                    xi.append((lc[0][3].x))
                                    yi.append((lc[0][3].y))
                                    l,m=m2l_utils.pts2dircos(lc[0][3].x,lc[0][3].y,rc[0][3].x,rc[0][3].y)
                                    if(not (l==0.0 and m==0.0)):
                                        fdc.append((l,m,'Fault_'+str(fault[c_l['o']])))
                                        all_coordsdist.append((dist))
                                        
        
                
                            
    fftc.close()
    structure_file=tmp_path+'structure_clip.shp'

    # first calculate interpolation for fault displacement calcs
    interpolate_orientations(structure_file,tmp_path,bbox,c_l,use_gcode,scheme,xi,yi,True)
    # then for near-fault calcs
    interpolate_orientations(structure_file,tmp_path+'ex_',bbox,c_l,use_gcode,scheme,all_coords_x,all_coords_y,True)    

    basal_contacts_file=tmp_path+'basal_contacts.shp'

    

    interpolate_contacts(basal_contacts_file,tmp_path,dtm,dtb,dtb_null,cover_map,bbox,c_l,use_gcode2,scheme,xi,yi,True)
    interpolate_contacts(basal_contacts_file,tmp_path+'ex_',dtm,dtb,dtb_null,cover_map,bbox,c_l,use_gcode2,scheme,all_coords_x,all_coords_y,True)    

    combo_file=tmp_path+'f_combo.csv'
    ex_combo_file=tmp_path+'ex_f_combo.csv'

    lc=np.loadtxt(tmp_path+'f_interpolation_contacts_l.csv',skiprows =1,delimiter =',',dtype=float)
    mc=np.loadtxt(tmp_path+'f_interpolation_contacts_m.csv',skiprows =1,delimiter =',',dtype=float)
    lo=np.loadtxt(tmp_path+'f_interpolation_l.csv',skiprows =1,delimiter =',',dtype=float)
    mo=np.loadtxt(tmp_path+'f_interpolation_m.csv',skiprows =1,delimiter =',',dtype=float)
    no=np.loadtxt(tmp_path+'f_interpolation_n.csv',skiprows =1,delimiter =',',dtype=float)
    xy=np.loadtxt(tmp_path+'f_interpolation_'+scheme+'.csv',skiprows =1,delimiter =',',dtype=float)

    ex_lc=np.loadtxt(tmp_path+'ex_f_interpolation_contacts_l.csv',skiprows =1,delimiter =',',dtype=float)
    ex_mc=np.loadtxt(tmp_path+'ex_f_interpolation_contacts_m.csv',skiprows =1,delimiter =',',dtype=float)
    ex_lo=np.loadtxt(tmp_path+'ex_f_interpolation_l.csv',skiprows =1,delimiter =',',dtype=float)
    ex_mo=np.loadtxt(tmp_path+'ex_f_interpolation_m.csv',skiprows =1,delimiter =',',dtype=float)
    ex_no=np.loadtxt(tmp_path+'ex_f_interpolation_n.csv',skiprows =1,delimiter =',',dtype=float)
    ex_xy=np.loadtxt(tmp_path+'ex_f_interpolation_'+scheme+'.csv',skiprows =1,delimiter =',',dtype=float)

    join_contacts_and_orientations(combo_file,geology_file,tmp_path,dtm_reproj_file,dtb,dtb_null,cover_map,c_l,lo,mo,no,lc,mc,xy,dst_crs,bbox,True)
    join_contacts_and_orientations(ex_combo_file,geology_file,tmp_path+'ex_',dtm_reproj_file,dtb,dtb_null,cover_map,c_l,ex_lo,ex_mo,ex_no,ex_lc,ex_mc,ex_xy,dst_crs,bbox,True)

    # loop though all nodes and calulate true displacement based on local estimated dip
    
    ddd=pd.read_csv(tmp_path+'f_combo_full.csv')
    f=open(output_path+'fault_displacements3.csv','w')
    f.write('X,Y,fname,apparent_displacement,vertical_displacement\n')

    for i in range (len(fdc)):
        l,m,n=m2l_utils.ddd2dircos(ddd.iloc[i]['dip'],ddd.iloc[i]['azimuth'])
        lnorm=l/sqrt(pow(l,2)+pow(m,2))
        mnorm=m/sqrt(pow(l,2)+pow(m,2))
        dotproduct=fabs((fdc[i][0]*lnorm)+(fdc[i][1]*mnorm))
        ostr=str(xi[i])+','+str(yi[i])+','+str(fdc[i][2])+','+str(int(all_coordsdist[i]))+','+str(abs(int(all_coordsdist[i]*tan(radians(dotproduct*ddd.iloc[i]['dip'])))))+'\n'
        f.write(ostr)

    f.close()
    print('fault displacement estimates saved as',output_path+'fault_displacements3.csv')
    print('near-fault orientations saved as',tmp_path+'ex_f_combo_full.csv')


def call_interpolator_grid(calc,x,y,l,m,n,xi,yi):
    # Calculate IDW or other interpolators

      
    if(calc=='simple_idw'):
        ZIl = simple_idw(x,y,l,xi,yi)
    if(calc=='scipy_rbf'):
        ZIl = scipy_rbf(x,y,l,xi,yi)
    if(calc=='scipy_idw'):
        ZIl = scipy_idw(x,y,l,xi,yi)

    
    if(calc=='simple_idw'):
        ZIm = simple_idw(x,y,m,xi,yi)
    if(calc=='scipy_rbf'):
        ZIm = scipy_rbf(x,y,m,xi,yi)
    if(calc=='scipy_idw'):
        ZIm = scipy_idw(x,y,m,xi,yi)

    
    if(type(n) is not int):
        if(calc=='simple_idw'):
            ZIn = simple_idw(x,y,n,xi,yi)
        if(calc=='scipy_rbf'):
            ZIn = scipy_rbf(x,y,n,xi,yi)
        if(calc=='scipy_idw'):
            ZIn = scipy_idw(x,y,n,xi,yi)
   
    else:
        ZIn=0
    return(ZIl,ZIm,ZIn)

def interpolate_orientation_grid(structures,calc,xcoords,ycoords,c_l):

    npts=len(structures)
    x = np.zeros(npts)
    y = np.zeros(npts)
    dip = np.zeros(npts)
    dipdir = np.zeros(npts)
    
    i=0
    for a_pt in structures.iterrows():
        x[i]=a_pt[1]['geometry'].x+(np.random.ranf()*0.01)
        y[i]=a_pt[1]['geometry'].y+(np.random.ranf()*0.01)
        dip[i] = a_pt[1][c_l['d']]
        
        if(c_l['otype']=='strike'):
            dipdir[i] = a_pt[1][c_l['dd']]+90
        else:
            dipdir[i] = a_pt[1][c_l['dd']]
            
        if(c_l['bo']==c_l['btype']):
            dipdir[i] = a_pt[1][c_l['dd']]+180
            dip[i]=-dip[i]
        i=i+1
  
    l=np.zeros(npts)
    m=np.zeros(npts)
    n=np.zeros(npts)
    
    for i in range(0,npts):
        l[i],m[i],n[i]=m2l_utils.ddd2dircos(dip[i],dipdir[i])

    ZIl,ZIm,ZIn=call_interpolator_grid(calc,x,y,l,m,n,xcoords,ycoords)

    l2=ZIl/np.sqrt(ZIl**2+ZIm**2+ZIn**2)
    m2=ZIm/np.sqrt(ZIl**2+ZIm**2+ZIn**2)
    n2=ZIn/np.sqrt(ZIl**2+ZIm**2+ZIn**2)

    dip = 90.0-np.degrees(np.arcsin(n2))
    dip_direction=np.where(m2>0, (360+np.degrees(np.arctan(l2/m2)))%360, 1)
    dip_direction=np.where(m2<0, (540+np.degrees(np.arctan(l2/m2)))%360, dip_direction)
    dip_direction=np.where(m2==0, 90, dip_direction)
    
    return(l2,m2,n2,dip,dip_direction)

def is_odd(num):
    return num & 0x1

def pts2dircos_arr(p1x,p1y,p2x,p2y):
    dx=p1x-p2x
    dy=p1y-p2y
    return(dx,dy)

def interpolate_contacts_grid(contacts,calc,xcoords_group,ycoords_group):
    decimate=1
    i=0
    listarray = []
    for indx,acontact in contacts.iterrows():   #loop through distinct linestrings in MultiLineString
        if(acontact.geometry.type=='MultiLineString'):
            for line in acontact.geometry:
                if(m2l_utils.mod_safe(i,decimate)  ==0):
                    listarray.append([line.coords[0][0], line.coords[0][1]])
                i=i+1
        else:
            if(  m2l_utils.mod_safe(i,decimate)  ==0):
                listarray.append([acontact.geometry.coords[0][0], acontact.geometry.coords[0][1]])
            i=i+1

    coords = np.array(listarray)
    if(is_odd(len(coords))):
        coords2=coords[:len(coords)-1,:].reshape((int(len(coords)/2),4))
    else:
        coords2=coords[:len(coords),:].reshape((int(len(coords)/2),4))

    dx=coords2[:,0:1]-coords2[:,2:3]
    dy=coords2[:,1:2]-coords2[:,3:4]
    x=coords2[:,0:1]+(dx/2)
    y=coords2[:,1:2]+(dy/2)

    dx2=dx**2
    dy2=dy**2
    mask=np.logical_and(dx2>0,dy2>0)


    dx=dx[mask]
    dy=dy[mask]
    dx2=dx2[mask]
    dy2=dy2[mask]
    x=x[mask]
    y=y[mask]
    scale=np.sqrt(dx2+dy2)
    l=dx/scale
    m=dy/scale
    #m=np.where(l<0, -m, m)
    #l=np.where(l<0, -l, l)
    
    if(len(x)>2):
        ZIl,ZIm,ZIn=call_interpolator_grid(calc,x,y,l,m,0,xcoords_group,ycoords_group)
        l2=ZIl/np.sqrt(ZIl**2+ZIm**2)
        m2=ZIm/np.sqrt(ZIl**2+ZIm**2)
        S=np.degrees(np.arctan2(l2,m2))
                        
        return(l2,m2,S)
    else:
        return(0,0,0)
    
    

def interpolation_grids(geology_file,structure_file,basal_contacts,bbox,spacing,dst_crs,scheme,super_groups,c_l):
    geology = gpd.read_file(geology_file,bbox=bbox)
    structures_code = gpd.read_file(structure_file,bbox=bbox)
    contacts = gpd.read_file(basal_contacts,bbox=bbox)

    geology[c_l['g']].fillna(geology[c_l['g2']], inplace=True)
    geology[c_l['g']].fillna(geology[c_l['c']], inplace=True)
    
    if(spacing<0):
        spacing=-(bbox[2]-bbox[0])/spacing
    
    x=(bbox[2]-bbox[0])/spacing
    y=(bbox[3]-bbox[1])/spacing

    xcoords=np.arange(bbox[0],bbox[2],spacing)
    ycoords=np.arange(bbox[1],bbox[3],spacing)
    xcoords, ycoords = np.meshgrid(xcoords, ycoords)
    xcoords, ycoords = xcoords.flatten(), ycoords.flatten()

    xycoords=np.vstack((xcoords,ycoords)).transpose()
    xycoords=xycoords.reshape(len(xcoords),2)

    nodes=gpd.GeoDataFrame(xycoords, geometry=[Point(xy) for xy in zip(xcoords, ycoords)])
    nodes.crs=dst_crs
    nodes_code = gpd.sjoin(nodes, geology, how="left", op="within")
    #structures_code = gpd.sjoin(structures, geology, how="left", op="within")

    first_supergroup=True
    for groups in super_groups:
        first=True
        for group in groups:
                        
            if(first):
                all_nodes=nodes_code[nodes_code[c_l['g']]==group] 
                all_structures = structures_code[structures_code[c_l['g']]==group]
                all_contacts=contacts[contacts[c_l['g']]==group.replace(' ','_').replace('-','_')]
                first=False
            else:
                another_node=nodes_code[nodes_code[c_l['g']]==group] 
                all_nodes=pd.concat([all_nodes,another_node],sort=False)
                                      
                another_contact=contacts[contacts[c_l['g']]==group.replace(' ','_').replace('-','_')]
                all_contacts=pd.concat([all_contacts,another_contact],sort=False)
                      
                another_structure = structures_code[structures_code[c_l['g']]==group]
                all_structures=pd.concat([all_structures,another_structure],sort=False)
            
        xcoords_group=all_nodes.geometry.x
        ycoords_group=all_nodes.geometry.y
        
        if(len(xcoords_group)>0):
            if(len(all_structures)>2):
                l,m,n,d,dd=interpolate_orientation_grid(all_structures,scheme,xcoords_group,ycoords_group,c_l)
                xy_lmn=np.vstack((xcoords_group,ycoords_group,l,m,n,d,dd)).transpose()
                xy_lmn=xy_lmn.reshape(len(l),7)
            else:
                print(groups,'has no structures')

                xy_lmn=np.zeros((5,len(xcoords_group)))    
                xy_lmn=np.vstack((xcoords_group,ycoords_group,xy_lmn)).transpose()
                xy_lmn=xy_lmn.reshape(len(xcoords_group),7)
                                             
            if(len(all_contacts)>0):
                l,m,S=interpolate_contacts_grid(all_contacts,scheme,xcoords_group,ycoords_group)  
                if(type(l) is not int):                   
                    xy_lm_contacts=np.vstack((xcoords_group,ycoords_group,l,m,S)).transpose()
                    xy_lm_contacts=xy_lm_contacts.reshape(len(l),5)
                else:
                    xy_lm_contacts=np.zeros((3,len(xcoords_group)))   
                    xy_lm_contacts=np.vstack((xcoords_group,ycoords_group,xy_lm_contacts)).transpose()
                    xy_lm_contacts=xy_lm_contacts.reshape(len(xcoords_group),5)
            else:
                print(groups,'has no contacts')

                xy_lm_contacts=np.zeros((3,len(xcoords_group)))   
                xy_lm_contacts=np.vstack((xcoords_group,ycoords_group,xy_lm_contacts)).transpose()
                xy_lm_contacts=xy_lm_contacts.reshape(len(xcoords_group),5)
    
            if(first_supergroup):
                first_supergroup=False
                xy_lmn_all=np.copy(xy_lmn)
                xy_lm_contacts_all=np.copy(xy_lm_contacts)
            else:
                xy_lmn_all=np.vstack((xy_lmn_all,xy_lmn))
                xy_lm_contacts_all=np.vstack((xy_lm_contacts_all,xy_lm_contacts))
    
    # sort to get back to x,y grid ordering
    dt = [('x', xy_lmn_all.dtype),('y', xy_lmn_all.dtype),('l', xy_lmn_all.dtype),('m', xy_lmn_all.dtype),('n', xy_lmn_all.dtype),('d', xy_lmn_all.dtype),('dd', xy_lmn_all.dtype)]
    #assert xy_lmn_all.flags['C_CONTIGUOUS']
    orientation_interp = xy_lmn_all.ravel().view(dt)
    orientation_interp.sort(order=['x','y','l','m','n','d','dd'])
        
    dt = [('x', xy_lm_contacts_all.dtype),('y', xy_lm_contacts_all.dtype),('l', xy_lm_contacts_all.dtype),('m', xy_lm_contacts_all.dtype),('S', xy_lm_contacts_all.dtype)]
    #assert xy_lm_contacts_all.flags['C_CONTIGUOUS']
    contact_interp = xy_lm_contacts_all.ravel().view(dt)
    contact_interp.sort(order=['x','y','l','m','S'])
    
    scale=np.sqrt(1-(orientation_interp['n']**2))
    lscaled=-scale*contact_interp['m']
    mscaled=scale*contact_interp['l']
    scale2=np.sqrt(orientation_interp['l']**2+orientation_interp['m']**2)
    loscaled=np.where(scale2>0, contact_interp['l']/scale2, 0)
    moscaled=np.where(scale2>0, contact_interp['m']/scale2, 0)
    dotproduct=(-contact_interp['m']*loscaled)+(contact_interp['l']*moscaled)
    
    #lscaled=np.where(dotproduct<0, -lscaled,lscaled)    
    #mscaled=np.where(dotproduct<0, -mscaled,mscaled)
    dotproduct=np.where(dotproduct>1,1,dotproduct)
    dotproduct=np.where(dotproduct<-1,-1,dotproduct)
    misorientation=np.degrees(np.arccos(dotproduct))
    
    dip = 90.0-np.degrees(np.arcsin(orientation_interp['n']))
    mscaled=np.where(mscaled==0,1e-5,mscaled)
    #dip_direction=np.where(mscaled>0, (360+np.degrees(np.arctan2(lscaled,mscaled)))%360, 1)
    dip_direction=(360+np.degrees(np.arctan2(lscaled,mscaled)))%360
    #dip_direction=np.where(mscaled<0, (540+np.degrees(np.arctan2(lscaled,mscaled)))%360, dip_direction)
    #dip_direction=np.where(mscaled==0, 90, dip_direction)

    combo_interp=np.vstack((contact_interp['x'],contact_interp['y'],lscaled,mscaled,orientation_interp['n'],
                            dip,dip_direction)).transpose()

    return(orientation_interp,contact_interp,combo_interp)


def process_fault_throw_and_near_faults_from_grid(tmp_path,output_path,dtm_reproj_file,dtb,dtb_null,cover_map,c_l,dst_crs,bbox,scheme,dip_grid,dip_dir_grid,xgrid,ygrid,spacing):
    fault_file=tmp_path+'faults_clip.shp'
    geology_file=tmp_path+'geol_clip.shp'

    faults = gpd.read_file(fault_file)
    faults=faults.dropna(subset=['geometry'])
    geology = gpd.read_file(geology_file)
    dtm = rasterio.open(dtm_reproj_file)

    all_long_faults=np.genfromtxt(output_path+'fault_dimensions.csv',delimiter=',',dtype='U100')
    fault_names=all_long_faults[1:,:1]
    m_step=5 #outstep from fault
    decimate_near=50
    xi=[]
    yi=[]
    fdc=[]
    all_coordsdist=[]
    all_coords_x=[]
    all_coords_y=[]

    fftc=open(output_path+'fault_tip_contacts.csv','w')
    fftc.write('X,Y,Z,formation\n')
    
    # loop through all faults
    
    for index,fault in faults.iterrows():
        #if(fault[c_l['o']]!=1071):
            #continue
        if(not str(fault.geometry.type) =='None'):
            if(fault.geometry.type=='LineString'):
                if('Fault_'+str(fault[c_l['o']]) in fault_names): # in is dangerous as Fault_1 is in Fault_10
                    #print('LineString','Fault_'+str(fault[c_l['o']]),len(fault.geometry.coords))
                    lcoords=[]
                    rcoords=[]
                    index=[]
                    
                    # make a list of points just offset from mid-points between fault nodes and join
                    # geology polygon information to points
                    j=0                
                    for i in range (0,len(fault.geometry.coords)-1):
                        l,m=m2l_utils.pts2dircos(fault.geometry.coords[i][0],fault.geometry.coords[i][1],fault.geometry.coords[i+1][0],fault.geometry.coords[i+1][1])
                        dx=m_step*m
                        dy=m_step*l
                        for inc in np.arange(0.1, 1, 0.01): 
                            midx=fault.geometry.coords[i][0]+((fault.geometry.coords[i+1][0]-fault.geometry.coords[i][0])*inc)            
                            midy=fault.geometry.coords[i][1]+((fault.geometry.coords[i+1][1]-fault.geometry.coords[i][1])*inc) 
                            lcoords.append([(midx+dx,midy-dy)])
                            rcoords.append([(midx-dx,midy+dy)])
                            index.append([(j)])
                            j=j+1
                    lgeom=[Point(xy) for xy in lcoords]        
                    rgeom=[Point(xy) for xy in rcoords]

                    lgdf = GeoDataFrame(index, crs=dst_crs, geometry=lgeom)
                    rgdf = GeoDataFrame(index, crs=dst_crs, geometry=rgeom)
                    lcode = gpd.sjoin(lgdf, geology, how="left", op="within")        
                    rcode = gpd.sjoin(rgdf, geology, how="left", op="within")
                    #display(lcode)
                    
                    # add lots of points left and right of fault to make sure ellipse range is happy in geomodeller
                    lgroups=[]
                    for ind, indl in lcode.iterrows():
                        if(not str(indl[c_l['c']]) == 'nan' and not str(indl[c_l['r1']]) == 'nan'):
                            if(m2l_utils.mod_safe(ind, decimate_near) == 0 or ind == len(lcode)-1):
                                if((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']])):
                                    locations = [
                                        (indl.geometry.x, indl.geometry.y)]
                                    last_height_l = m2l_utils.value_from_dtm_dtb(
                                        dtm, dtb, dtb_null, cover_map, locations)
                                    if(not indl[c_l['g']] in lgroups):
                                        ostr = "{},{},{},{}\n"\
                                            .format(indl.geometry.x, indl.geometry.y, last_height_l, indl[c_l['c']].replace(" ", "_").replace("-", "_"))
                                        fftc.write(ostr)
                                        lgroups.append(indl[c_l['g']])
                    rgroups=[]
                    for ind, indr in rcode.iterrows():
                        if(not str(indr[c_l['c']]) == 'nan' and not str(indr[c_l['r1']]) == 'nan'):
                            if((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']])):
                                if(m2l_utils.mod_safe(ind, decimate_near) == 0 or ind == len(rcode)-1):
                                    locations = [
                                        (indr.geometry.x, indr.geometry.y)]
                                    last_height_r = m2l_utils.value_from_dtm_dtb(
                                        dtm, dtb, dtb_null, cover_map, locations)
                                    if(not indr[c_l['g']] in rgroups):
                                        ostr = "{},{},{},{}\n"\
                                            .format(indr.geometry.x, indr.geometry.y, last_height_r, indr[c_l['c']].replace(" ", "_").replace("-", "_"))
                                        fftc.write(ostr)
                                        rgroups.append(indr[c_l['g']])

                    # add points to list if they have different geology code than previous node on left side
                    
                    first=True
                    lcontact=[]
                    lastlcode=''
                    
                    for ind,indl in lcode.iterrows():                       
                        if(ind<len(lcode) and not isnan(indl['index_right'])):
                            ntest1=str(indl[c_l['ds']])
                            ntest2=str(indl[c_l['r1']])
                            
                            if(not ntest1 == 'None' and not ntest2 == 'None'):
                                if(ind==1 or (not lastlcode==indl[c_l['c']] and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']]) ))):
                                    all_coords_x.append(indl.geometry.x)
                                    all_coords_y.append(indl.geometry.y)
                                    if(first):
                                        first=False
                                        firstlx=indl.geometry.x
                                        firstly=indl.geometry.y
                                        firstlc=indl[c_l['c']].replace(" ","_").replace("-","_")
                                    lastlx=indl.geometry.x
                                    lastly=indl.geometry.y
                                    lastlc=indl[c_l['c']].replace(" ","_").replace("-","_")
                                    
        
                                if(lastlcode=='' and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']] ))):
                                    lastlcode=indl[c_l['c']]
        
                                #print('l',ind,indl[c_l['c']],indl[c_l['ds']],indl[c_l['r1']])
                                if(not ntest1 == 'None' and not ntest2 == 'None' and not str(indl[c_l['c']])=='nan'):
                                    if((not indl[c_l['c']]==lastlcode) and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']] ))):
                                        lcontact.append([(ind,lastlcode,indl[c_l['c']],indl.geometry)])
                                        lastlcode=indl[c_l['c']]
                                        locations=[(lastlx,lastly)]
                                        last_height_l=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                                        ostr="{},{},{},{}\n"\
                                                    .format(lastlx,lastly,last_height_l,indl[c_l['c']].replace(" ","_").replace("-","_"))
                                        #fftc.write(ostr)

                                        
                        
                    #add points to list if they have different geology code than previous node on right side
                    
                    first=True
                    rcontact=[]
                    lastrcode=''
                    for ind,indr in rcode.iterrows():
                        if(ind<len(rcode) and not isnan(indr['index_right'])):
                            ntest1=str(indr[c_l['ds']])
                            ntest2=str(indr[c_l['r1']])
                            
                            if(not ntest1 == 'None' and not ntest2 == 'None'):                            
                                if(ind==1 or (not lastrcode==indr[c_l['c']] and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']]) ))):
                                    all_coords_x.append(indr.geometry.x)
                                    all_coords_y.append(indr.geometry.y)
                                    if(first):
                                        first=False
                                        firstrx=indr.geometry.x
                                        firstry=indr.geometry.y
                                        firstrc=indr[c_l['c']].replace(" ","_").replace("-","_")
                                    lastrx=indr.geometry.x
                                    lastry=indr.geometry.y
                                    lastrc=indr[c_l['c']].replace(" ","_").replace("-","_")
        
                                if(lastrcode=='' and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']] ))):
                                    lastrcode=indr[c_l['c']]                        
                                #print('r',ind,indr[c_l['c']],indr[c_l['ds']],indr[c_l['r1']])
        
                                #print(lastrcode,ntest1,ntest2,str(indr[c_l['c']]),indr[c_l['c']],c_l['sill'],c_l['intrusive'])
                                
                                if(not ntest1 == 'None' and not ntest2 == 'None' and not str(indr[c_l['c']])=='nan' ):
                                    if((not indr[c_l['c']]==lastrcode) and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']] ))):
                                        rcontact.append([(ind,lastrcode,indr[c_l['c']],indr.geometry)]) 
                                        lastrcode=indr[c_l['c']]             
                                        locations=[(lastrx,lastry)]
                                        last_height_r=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                                        ostr="{},{},{},{}\n"\
                                                    .format(lastrx,lastry,last_height_r,indr[c_l['c']].replace(" ","_").replace("-","_"))
                                        #fftc.write(ostr)
                    
                                
                    # loop through left and right sides to find equivalent contact pairs along fault
                    
                    if(len(lcontact)>0 and len(rcontact)>0):                
                        for lc in lcontact:
                            for rc in rcontact:
                                #display('l',lc[0][3].x,'r',rc[0][3].x)
                                if(lc[0][1]==rc[0][1] and lc[0][2]==rc[0][2] and not lc[0][1]==''):
                                    dist=m2l_utils.ptsdist(lc[0][3].x,lc[0][3].y,rc[0][3].x,rc[0][3].y)
                                    if(lc[0][0]<rc[0][0]):
                                        dist=-dist
                                    #print('***',lc,rc)

                                    xi.append((lc[0][3].x))
                                    yi.append((lc[0][3].y))
                                    l,m=m2l_utils.pts2dircos(lc[0][3].x,lc[0][3].y,rc[0][3].x,rc[0][3].y)
                                    if(not (l==0.0 and m==0.0)):
                                        fdc.append((l,m,'Fault_'+str(fault[c_l['o']])))
                                        all_coordsdist.append((dist))
            else:
                if('Fault_'+str(fault[c_l['o']]) in fault_names): # in is dangerous as Fault_1 is in Fault_10
                    for fls in fault.geometry:              
                        fault_ls=LineString(fls)        
                        lcoords=[]
                        rcoords=[]
                        index=[]
                        #display("MLS DEBUG",fault.geometry.type)
                        
                        j=0                
                        for i in range (0,len(fault_ls.coords)-1):
                            for inc in np.arange(0.01, 1, 0.01): 
                                midx=fault_ls.coords[i][0]+((fault_ls.coords[i+1][0]-fault_ls.coords[i][0])*inc)            
                                midy=fault_ls.coords[i][1]+((fault_ls.coords[i+1][1]-fault_ls.coords[i][1])*inc)
                                l,m=m2l_utils.pts2dircos(fault_ls.coords[i][0],fault_ls.coords[i][1],fault_ls.coords[i+1][0],fault_ls.coords[i+1][1])
                                lcoords.append([(midx+(m_step*m),midy-(m_step*l))])
                                rcoords.append([(midx-(m_step*m),midy+(m_step*l))])
                                index.append([(j)])
                                j=j+1

                        lgeom=[Point(xy) for xy in lcoords]        
                        rgeom=[Point(xy) for xy in rcoords]
                        lgdf = GeoDataFrame(index, crs=dst_crs, geometry=lgeom)
                        rgdf = GeoDataFrame(index, crs=dst_crs, geometry=rgeom)
                        lcode = gpd.sjoin(lgdf, geology, how="left", op="within")        
                        rcode = gpd.sjoin(rgdf, geology, how="left", op="within")

                        #add points to list if they have different geology code than previous node on left side

                        first=True
                        lcontact=[]
                        lastlcode=''
                        for ind,indl in lcode.iterrows():

                            if(ind<len(lcode) and not isnan(indl['index_right'])):
                                ntest1=str(indl[c_l['ds']])
                                ntest2=str(indl[c_l['r1']])
                                
                                if(not ntest1 == 'None' and not ntest2 == 'None'):
                                    if(ind==1 or (not lastlcode==indl[c_l['c']] and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']]) ))):
                                        all_coords_x.append(indl.geometry.x)
                                        all_coords_y.append(indl.geometry.y)
                                        if(first):
                                            first=False
                                            firstlx=indl.geometry.x
                                            firstly=indl.geometry.y
                                            firstlc=indl[c_l['c']].replace(" ","_").replace("-","_")
                                        lastlx=indl.geometry.x
                                        lastly=indl.geometry.y
                                        lastlc=indl[c_l['c']].replace(" ","_").replace("-","_")
        
                                    if(lastlcode=='' and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']] ))):
                                        lastlcode=indl[c_l['c']]
        
                                    if(not ntest1 == 'None' and not ntest2 == 'None' and not str(indl[c_l['c']])=='nan'):
                                        if((not indl[c_l['c']]==lastlcode) and ((not c_l['sill'] in indl[c_l['ds']]) or (not c_l['intrusive'] in indl[c_l['r1']] ))):
                                            lcontact.append([(ind,lastlcode,indl[c_l['c']],indl.geometry)])
                                            lastlcode=indl[c_l['c']]
                                            locations=[(lastlx,lastly)]
                                            last_height_l=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                                            ostr="{},{},{},{}\n"\
                                                        .format(lastlx,lastly,last_height_l,indl[c_l['c']].replace(" ","_").replace("-","_"))
                                            fftc.write(ostr)
                                            

                        #add points to list if they have different geology code than previous node on right side

                        first=True
                        rcontact=[]
                        lastrcode=''
                        for ind,indr in rcode.iterrows():
                            if(ind<len(rcode) and not isnan(indr['index_right'])):
                                ntest1=str(indr[c_l['ds']])
                                ntest2=str(indr[c_l['r1']])
                                
                                if(not ntest1 == 'None' and not ntest2 == 'None'):
                                    if(ind==1 or (not lastrcode==indr[c_l['c']] and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']]) ))):
                                        all_coords_x.append(indr.geometry.x)
                                        all_coords_y.append(indr.geometry.y)
                                        if(first):
                                            first=False
                                            firstrx=indr.geometry.x
                                            firstry=indr.geometry.y
                                            firstrc=indr[c_l['c']].replace(" ","_").replace("-","_")
                                        lastrx=indr.geometry.x
                                        lastry=indr.geometry.y
                                        lastrc=indr[c_l['c']].replace(" ","_").replace("-","_")
        
                                    if(lastrcode=='' and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']] ))):
                                        lastrcode=indr[c_l['c']]                        
        
                                    if(not ntest1 == 'None' and not ntest2 == 'None' and not str(indr[c_l['c']])=='nan' ):
                                        if((not indr[c_l['c']]==lastrcode) and ((not c_l['sill'] in indr[c_l['ds']]) or (not c_l['intrusive'] in indr[c_l['r1']] ))):
                                            rcontact.append([(ind,lastrcode,indr[c_l['c']],indr.geometry)]) 
                                            lastrcode=indr[c_l['c']]
                                            locations=[(lastrx,lastry)]
                                            last_height_r=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                                            ostr="{},{},{},{}\n"\
                                                        .format(lastrx,lastry,last_height_r,indr[c_l['c']].replace(" ","_").replace("-","_"))
                                            fftc.write(ostr)
        
                        locations=[(firstlx,firstly)]
                        first_height_l=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                        ostr="{},{},{},{}\n"\
                                    .format(firstlx,firstly,first_height_l,firstlc.replace(" ","_").replace("-","_"))
                        fftc.write(ostr)
                        locations=[(firstrx,firstry)]
                        first_height_r=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                        ostr="{},{},{},{}\n"\
                                    .format(firstrx,firstry,first_height_r,firstrc.replace(" ","_").replace("-","_"))
                        fftc.write(ostr)
                        locations=[(lastlx,lastly)]
                        last_height_l=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                        ostr="{},{},{},{}\n"\
                                    .format(lastlx,lastly,last_height_l,lastlc.replace(" ","_").replace("-","_"))
                        fftc.write(ostr)
                        locations=[(lastrx,lastry)]
                        last_height_r=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                        ostr="{},{},{},{}\n"\
                                    .format(lastrx,lastry,last_height_r,lastrc.replace(" ","_").replace("-","_"))
                        fftc.write(ostr)   
                        if(len(lcode)>5):                        
                            locations=[(lcode.iloc[len(lcode)-3].geometry.x,lcode.iloc[len(lcode)-3].geometry.y)]
                            last_height_l=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                            ostr="{},{},{},{}\n"\
                                        .format(lcode.iloc[len(lcode)-3].geometry.x,lcode.iloc[len(lcode)-3].geometry.y,last_height_l,str(lcode.iloc[len(lcode)-3][c_l['c']]).replace(" ","_").replace("-","_"))
                            if(not str(lcode.iloc[len(lcode)-3][c_l['c']])=='nan'):
                                fftc.write(ostr)                     
                        if(len(rcode)>5):    
                            locations=[(rcode.iloc[len(rcode)-3].geometry.x,rcode.iloc[len(rcode)-3].geometry.y)]
                            last_height_r=m2l_utils.value_from_dtm_dtb(dtm,dtb,dtb_null,cover_map,locations)
                            ostr="{},{},{},{}\n"\
                                        .format(lcode.iloc[len(rcode)-3].geometry.x,rcode.iloc[len(rcode)-3].geometry.y,last_height_r,str(rcode.iloc[len(rcode)-3][c_l['c']]).replace(" ","_").replace("-","_"))
                            if(not str(rcode.iloc[len(rcode)-3][c_l['c']])=='nan'):
                                fftc.write(ostr)                     
            
                        # loop through left and right sides to find equivalent contact pairs along fault

                        if(len(lcontact)>0 and len(rcontact)>0):
                            for lc in lcontact:
                                for rc in rcontact:
                                    #display('l',lc[0][3].x,'r',rc[0][3].x)
                                    if(lc[0][1]==rc[0][1] and lc[0][2]==rc[0][2] and not lc[0][1]==''):
                                        dist=m2l_utils.ptsdist(lc[0][3].x,lc[0][3].y,rc[0][3].x,rc[0][3].y)
                                        if(lc[0][0]<rc[0][0]):
                                            dist=-dist
                                        #print('***',lc,rc)

                                        xi.append((lc[0][3].x))
                                        yi.append((lc[0][3].y))
                                        l,m=m2l_utils.pts2dircos(lc[0][3].x,lc[0][3].y,rc[0][3].x,rc[0][3].y)
                                        if(not (l==0.0 and m==0.0)):
                                            fdc.append((l,m,'Fault_'+str(fault[c_l['o']])))
                                            all_coordsdist.append((dist))
                                            
            
                
                            
    fftc.close()
    
    fault_orien=pd.read_csv(output_path+'fault_orientations.csv',",")
    fault_orien.set_index('formation',  inplace = True)
    f=open(output_path+'fault_displacements3.csv','w')
    f.write('X,Y,fname,apparent_displacement,vertical_displacement,downthrow_dir\n')

    for i in range (len(fdc)):
        r=int((yi[i]-bbox[1])/spacing)
        c=int((xi[i]-bbox[0])/spacing)

        l,m,n=m2l_utils.ddd2dircos(dip_grid[r,c],dip_dir_grid[r,c])
        lnorm=l/sqrt(pow(l,2)+pow(m,2))
        mnorm=m/sqrt(pow(l,2)+pow(m,2))
        
        # product of (rotation sign between bed dip direction and fault dip direction) 
        # and appararent throw sign gives downthrown side relative to dip direction of fault
        lf,mf,nf=m2l_utils.ddd2dircos(fault_orien.loc[fdc[i][2]]['dip'],fault_orien.loc[fdc[i][2]]['DipDirection'])
        lfnorm=lf/sqrt(pow(lf,2)+pow(mf,2))
        mfnorm=mf/sqrt(pow(lf,2)+pow(mf,2))
        
        dotproduct=fabs((fdc[i][0]*lnorm)+(fdc[i][1]*mnorm))
        crossproduct= asin(lnorm * mfnorm - mnorm * lfnorm)
        if(crossproduct*all_coordsdist[i]<0): 
            down_dipdir=fault_orien.loc[fdc[i][2]]['DipDirection']
        else:
            down_dipdir=(fault_orien.loc[fdc[i][2]]['DipDirection']+180.0)%360

            
        ostr=str(xi[i])+','+str(yi[i])+','+str(fdc[i][2])+','+str(int(all_coordsdist[i]))+','+str(abs(int(all_coordsdist[i]*tan(radians(dotproduct*dip_grid[r,c])))))+','+str(down_dipdir)+'\n'
        f.write(ostr)

    f.close()
    print('fault displacement estimates saved as',output_path+'fault_displacements3.csv')
    print('near-fault orientations saved as',tmp_path+'ex_f_combo_full.csv')
    print('near-fault orientations saved as',tmp_path+'ex_f_combo_full.csv')  
    
    # when no fault displacment data are available, set to 1m displacment so they still are calculated
    
    if(os.path.exists(fault_file)):
        fault_ori=pd.read_csv(output_path+'fault_orientations.csv')
        fault_disp=pd.read_csv(output_path+'fault_displacements3.csv')
        
        fault_disp_found=fault_disp['fname'].unique()
        
        f=open(output_path+'fault_displacements3.csv','a')
        
        for ind,fault in fault_ori.iterrows():
            if(not fault['formation'] in fault_disp_found):
                ostr=str(fault['X'])+','+str(fault['Y'])+','+str(fault['formation'])+','+str(1)+','+str(1)+','+str(1)+'\n'
                f.write(ostr)
    
        f.close()