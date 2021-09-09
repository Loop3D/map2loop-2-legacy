import geopandas as gpd
import matplotlib.pyplot as plt


def plot_geochron(geochron_filename,bbox):
    geochron=gpd.read_file(geochron_filename,bbox=bbox)
    geochron_sort=geochron.sort_values(by='ANALYSIS_A')
    fig = plt.figure(figsize=[7,7])
    plt.gca().set_ylim(0, 4000)
    xscale=200
    i=0
    for ind,date in geochron_sort.iterrows():
        if('Ar-Ar' in date['ANALYSIS_T']):
            c='r'
            label='Ar-Ar'
        elif('U-Pb' in date['ANALYSIS_T']):
            c='b'
            label='U-Pb'
        elif('K-Ar' in date['ANALYSIS_T']):
            c='g'
            label='K-Ar'
        else:
            c='k'
            label='other'

        if('igneous' in date['ANALYSIS_E']):
            s='o'
        elif('metam' in date['ANALYSIS_E']):
            s='d'
        elif('xeno' in date['ANALYSIS_E']):
            s='s'
        elif('cool' in date['ANALYSIS_E']):
            s='v'
        elif('maximum age' in date['ANALYSIS_E']):
            s='1'
        elif('depositional' in date['ANALYSIS_E']):
            s='p'
        elif('minimum' in date['ANALYSIS_E']):
            s='+'
        else:
            s='.'

        line = plt.Line2D((i*xscale, i*xscale),(date['ANALYSIS_A']-date['ANALYSIS_U'], date['ANALYSIS_A']+date['ANALYSIS_U']), lw=1.0,color='#000000')
        plt.gca().add_line(line)
        plt.plot([i*xscale], [date['ANALYSIS_A']], c+s)

        i=i+1
    ax = plt.axis()
    plt.axis((ax[0],ax[1],ax[3],ax[2]))
    plt.title('Geochronology')
    plt.savefig("geochron.pdf")
    #geochron_sort.plot(column='ANALYSIS_E',figsize=(7,7),edgecolor='#000000',linewidth=0.2)
    #plt.savefig('geochron_map.pdf')  
    plt.show()  
    
def plot_geology_codes(geology_filename,bbox):
    geol=gpd.read_file(geology_filename,bbox=bbox)
    geol_sort=geol.sort_values(by='MAX_AGE_MA')
    geol_sort_unique=geol_sort.drop_duplicates(subset=['CODE'])

    plt.axes()
    xscale=200
    i=0
    print(len(geol_sort_unique))
    for ind,poly in geol_sort_unique.iterrows():
        if(not str(poly['MIN_AGE_MA'])=='None' and not str(poly['MAX_AGE_MA'])=='None'  ):
            ave_age=float(poly['MIN_AGE_MA'])+(float(poly['MAX_AGE_MA'])-float(poly['MIN_AGE_MA']))/2.0
            plt.plot([i*xscale], ave_age, 'ro')        
            line = plt.Line2D((i*xscale, i*xscale),(float(poly['MAX_AGE_MA']), float(poly['MIN_AGE_MA'])), lw=1.5)
            plt.gca().add_line(line)
            i=i+1
    ax = plt.axis()
    plt.axis((ax[0],ax[1],ax[3],ax[2]))
    plt.title('Strat Ages')
    plt.savefig('Strat_ages.pdf')  

    plt.show()   
     
def plot_geology_parents(geology_filename,bbox):
    colors=['b','g','r','c','m','y','k','w']
    fig = plt.figure(figsize=[7, 7])        
    geol=gpd.read_file(geology_filename,bbox=bbox)

    geol_unique=geol.drop_duplicates(subset=['CODE'])

    geol_unique_code=geol_unique.set_index('CODE')
    parents=geol_unique['PARENTNAME'].unique()
    print('parents=',len(parents),'codes=',len(geol_unique))
    minmax_parent={}
    for parent in parents:
        min=1e10
        max=0
        for code in geol_unique['CODE']:
            if(geol_unique_code.loc[code]['PARENTNAME']==parent):
                if(not str(geol_unique_code.loc[code]['MIN_AGE_MA'])=='None' and not str(geol_unique_code.loc[code]['MAX_AGE_MA'])=='None'  ):
                    if(float(geol_unique_code.loc[code]['MIN_AGE_MA'])<min):
                        min=float(geol_unique_code.loc[code]['MIN_AGE_MA'])        
                    if(float(geol_unique_code.loc[code]['MAX_AGE_MA'])>max):
                        max=float(geol_unique_code.loc[code]['MAX_AGE_MA'])
        if( not '_Top' in parent):
            minmax_parent[parent.replace("_","")]=[min,max]
    #display(minmax_parent)
    minmax_parent_sort=sorted(minmax_parent.items(), key = lambda kv:(kv[1], kv[0]))

    xscale=500
    i=0
    c=0
    for parent in minmax_parent_sort:
        key=parent[0]
        rectangle = plt.Rectangle((i, minmax_parent[key][0]), xscale, minmax_parent[key][1]-minmax_parent[key][0], fc='C'+str(c),edgecolor='k',label=key)
        plt.gca().add_patch(rectangle)
        i=i+xscale
        c=c+1

    plt.axis('scaled')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax = plt.axis()
    plt.axis((ax[0],ax[1],ax[3],ax[2]))
    plt.gca().set_ylim(4000, 0)
    plt.savefig('parents.pdf',orientation = 'portrait')  
    #geol.plot(column='PARENTNAME',figsize=(7,7),edgecolor='#000000',linewidth=0.2,legend=True,legend_kwds={'loc': [0,-.4]})

    #plt.savefig('parents_map.pdf',papertype = 'a4', orientation = 'portrait')  
    plt.title('Stratigraphy')

    plt.show()
    
def plot_orogenic(orogenic_filename,bbox):
    fig = plt.figure(figsize=[7,7])
    geol=gpd.read_file(orogenic_filename,bbox=bbox)

    geol_unique=geol.drop_duplicates(subset=['EVENTNAME'])
    minmax={}
    for ind,geol in geol_unique.iterrows():
        min=1e10
        max=0
        if(float(geol['min_age'])<min):
            min=float(geol['min_age'])        
        if(float(geol['max_age'])>max):
            max=float(geol['max_age'])
        minmax[geol['EVENTNAME']]=[min,max]
    #display(minmax_parent)
    minmax_sort=sorted(minmax.items(), key = lambda kv:(kv[1], kv[0]))

    xscale=500
    i=0
    c=0
    for parent in minmax_sort:
        key=parent[0]
        rectangle = plt.Rectangle((i, minmax[key][0]), xscale, minmax[key][1]-minmax[key][0], fc='C'+str(c),edgecolor='k',label=key)
        plt.gca().add_patch(rectangle)
        i=i+xscale
        c=c+1

    if(len(minmax_sort)>0):
        plt.axis('scaled')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax = plt.axis()
        plt.axis((ax[0],ax[1],ax[3],ax[2]))
        plt.gca().set_ylim(4000, 0)
        plt.savefig('orogenic.pdf', orientation = 'portrait')  
        #geol.plot(column='PARENTNAME',figsize=(7,7),edgecolor='#000000',linewidth=0.2,legend=True,legend_kwds={'loc': [0,-.4]})
    
        #plt.savefig('parents_map.pdf',papertype = 'a4', orientation = 'portrait')  
        plt.title('Orogenic Events')
    
        plt.show()
        
def plot_faults(linear_filename,bbox):
    import fracpaq as fpq 
    import matplotlib.pylab as plt 
    import sys 
    import numpy as np 
    linear=gpd.read_file(linear_filename,bbox=bbox)
    faults=linear[linear['FEATURE'].str.contains("Fault")]

    xnodelist=[]
    ynodelist=[]
    for inf,fault in faults.iterrows():
        xtmplist=[]
        ytmplist=[]
        for coords in fault['geometry'].coords:
            xtmplist.append(coords[0])
            ytmplist.append(coords[1])
        xnodelist.append(xtmplist)
        ynodelist.append(ytmplist)

    nodelist=[xnodelist,ynodelist]

    #   defaults 
    CM2INCHES = 0.3937 
    bGrid = False 
    bEqualArea = False 
    nBinWidth = 5  
    sColour = 'C0' 
    xSize = 15.0 * CM2INCHES 
    ySize = 15.0 * CM2INCHES   

    #   get nodes and calculate angles from nodes 

    xnodelist = nodelist[0] 
    ynodelist = nodelist[1]
    segangle = fpq.getSegAngles(xnodelist, ynodelist)
    nSegs = len(segangle)

    if(nSegs>0):

        #   bin the data and find maximum per bin 
        nBins = int(round(360/nBinWidth))
        segangleDoubled = np.zeros(len(segangle))
        segangleDoubled = np.copy(segangle)
        segangleDoubled = np.concatenate([segangleDoubled, segangleDoubled+180.0])
        n, b = plt.histogram(segangleDoubled, nBins)
        nMax = max(n) 

        #   plot the segment angle distribution 
        plt.figure(figsize=(xSize, ySize))
        plt.subplot(111, projection='polar')
        plt.bar(np.deg2rad(np.arange(0, 360, nBinWidth)), n, 
                width=np.deg2rad(nBinWidth), bottom=0.0, color='.8', edgecolor='k')
        plt.xticks(np.radians(range(0, 360, 45)), 
                   ['0', '45', '90', '135', '180', '215', '270', '315'])
        if(nMax>6):        
            plt.rgrids(range(0, int(round(nMax*1.1)), int(round((nMax*1.1)/5))), angle=330)
        plt.ylim(0, int(round(nMax*1.1)))
        plt.title('Segment strikes, n={}'.format(nSegs) )
        plt.savefig("faults_rose.pdf")
        #faults.plot(figsize=(7,7),edgecolor='#ff0000',linewidth=0.3)
        #print('Plotted %5d segments & angles' % nSegs)
        #plt.savefig('faults.pdf')  
        #plt.title('Fault traces')

        plt.show()

def plot_folds(linear_filename,bbox):
    import fracpaq as fpq 
    import matplotlib.pylab as plt 
    import sys 
    import numpy as np 

    linear=gpd.read_file(linear_filename,bbox=bbox)
    folds=linear[linear['FEATURE'].str.contains("Fold")]
    xnodelist=[]
    ynodelist=[]
    for inf,fold in folds.iterrows():
        xtmplist=[]
        ytmplist=[]
        for coords in fold['geometry'].coords:
            xtmplist.append(coords[0])
            ytmplist.append(coords[1])
        xnodelist.append(xtmplist)
        ynodelist.append(ytmplist)
    nodelist=[xnodelist,ynodelist]

    #   defaults 
    CM2INCHES = 0.3937 
    bGrid = False 
    bEqualArea = False 
    nBinWidth = 5  
    sColour = 'C0' 
    xSize = 15.0 * CM2INCHES 
    ySize = 15.0 * CM2INCHES   

    #   get nodes and calculate angles from nodes 

    xnodelist = nodelist[0] 
    ynodelist = nodelist[1]
    segangle = fpq.getSegAngles(xnodelist, ynodelist)
    nSegs = len(segangle)
    if(nSegs>0):
        #   bin the data and find maximum per bin 
        nBins = int(round(360/nBinWidth))
        segangleDoubled = np.zeros(len(segangle))
        segangleDoubled = np.copy(segangle)
        segangleDoubled = np.concatenate([segangleDoubled, segangleDoubled+180.0])
        n, b = plt.histogram(segangleDoubled, nBins)
        nMax = max(n) 

        #   plot the segment angle distribution 
        plt.figure(figsize=(xSize, ySize))
        plt.subplot(111, projection='polar')
        plt.bar(np.deg2rad(np.arange(0, 360, nBinWidth)), n, 
                width=np.deg2rad(nBinWidth), bottom=0.0, color='.8', edgecolor='k')

        plt.xticks(np.radians(range(0, 360, 45)), 
                   ['0', '45', '90', '135', '180', '215', '270', '315'])
        if(nMax>6):        
            plt.rgrids(range(0, int(round(nMax*1.1)), int(round((nMax*1.1)/5))), angle=330)
        plt.ylim(0, int(round(nMax*1.1)))
        plt.title('Segment strikes, n=%i' % nSegs)
        plt.savefig("folds_rose.pdf")
        #folds.plot(figsize=(7,7),edgecolor='#ff0000',linewidth=0.3)
        #print('Plotted %5d segments & angles' % nSegs)
        #plt.savefig('folds.pdf')  
        plt.title('Fold axial traces')
        plt.show()


def plot_mag(mag_netcdf_path,bbox):
    import os
    import netCDF4
    import numpy as np
    from geophys_utils import NetCDFGridUtils
    from geophys_utils import get_netcdf_edge_points, points2convex_hull
    import matplotlib.pyplot as plt
    
    netcdf_dataset = netCDF4.Dataset(mag_netcdf_path, 'r')

    max_bytes = 500000000

    #netcdf_grid_utils = NetCDFGridUtils(netcdf_dataset)
    lats = netcdf_dataset.variables['lat'][:]
    lons = netcdf_dataset.variables['lon'][:]
    latselect = np.logical_and(lats>bbox[1],lats<bbox[3])
    lonselect = np.logical_and(lons>bbox[0],lons<bbox[2])
    data = netcdf_dataset.variables['Band1'][latselect,lonselect]
    #plt.contourf(data[::-1],cmap ='gist_ncar') # flip latitudes so they go south -> north
    plt.imshow(data,cmap ='gist_rainbow',vmin=-1000,vmax=1000)
    plt.title('Magnetics')
    plt.show()

def plot_grav(grav_netcdf_path,bbox):
    import os
    import netCDF4
    import numpy as np
    from geophys_utils import NetCDFGridUtils
    from geophys_utils import get_netcdf_edge_points, points2convex_hull
    import matplotlib.pyplot as plt
    
    netcdf_dataset = netCDF4.Dataset(grav_netcdf_path, 'r')

    max_bytes = 500000000

    #netcdf_grid_utils = NetCDFGridUtils(netcdf_dataset)
    lats = netcdf_dataset.variables['lat'][:]
    lons = netcdf_dataset.variables['lon'][:]
    latselect = np.logical_and(lats>bbox[1],lats<bbox[3])
    lonselect = np.logical_and(lons>bbox[0],lons<bbox[2])
    data = netcdf_dataset.variables['Band1'][latselect,lonselect]
    plt.imshow(data,cmap ='gist_rainbow')
    plt.title('Gravity')
    plt.show()

def plot_eggs(eggs_filename,bbox):
    eggs=gpd.read_file(eggs_filename,bbox=bbox)
    if(len(eggs)>0):
        eggs.plot(figsize=(7,7),column='LOCATION_Z_IN_METRES')
    else:
        plt.text(.2,.2,"No eggs values\n in area")
    plt.title('EGGS')
    plt.show()

def plot_density(density_filename,bbox):
    density=gpd.read_file(density_filename,bbox=bbox)
    if(len(density)>0):
        density.plot(figsize=(7,7),column='value')
    else:
        plt.text(.2,.2,"No density values\n in area")
    plt.tight_layout()
    plt.title('Density')
    plt.show()

def plot_mag_sus(mag_sus_filename,bbox):
    mag_sus=gpd.read_file(mag_sus_filename,bbox=bbox)
    if(len(mag_sus)>0):
        mag_sus.plot(figsize=(7,7),column='lithologyg')
    else:
        plt.text(.2,.2,"No mag sus values\n in area")
    plt.tight_layout()
    plt.title('Magnetic Susceptibility')
    plt.show()

def plot_mag_mag_sus(mag_netcdf_path,mag_sus_filename,bbox):
    import os
    import netCDF4
    import numpy as np
    import geopandas as gpd
    from geophys_utils import NetCDFGridUtils
    from geophys_utils import get_netcdf_edge_points, points2convex_hull
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2)
    
    netcdf_dataset = netCDF4.Dataset(mag_netcdf_path, 'r')

    max_bytes = 500000000

    #netcdf_grid_utils = NetCDFGridUtils(netcdf_dataset)
    lats = netcdf_dataset.variables['lat'][:]
    lons = netcdf_dataset.variables['lon'][:]
    latselect = np.logical_and(lats>bbox[1],lats<bbox[3])
    lonselect = np.logical_and(lons>bbox[0],lons<bbox[2])
    data = netcdf_dataset.variables['Band1'][latselect,lonselect]
    #plt.contourf(data[::-1],cmap ='gist_ncar') # flip latitudes so they go south -> north
    ax[0].imshow(data,cmap ='gist_rainbow',vmin=-1000,vmax=1000)
    
    mag_sus=gpd.read_file(mag_sus_filename,bbox=bbox)
    if(len(mag_sus)>0):
        mag_sus.plot(ax=ax[1],figsize=(7,7),column='lithologyg')
    else:
        plt.subplot(122)
        plt.text(.2,.2,"No mag sus values\n in area")
    plt.tight_layout()
    plt.title('Magnetics and Magnetic Susceptibility')

    plt.show()
    
def plot_grav_density(grav_netcdf_path,density_filename,bbox):
    import os
    import netCDF4
    import numpy as np
    import geopandas as gpd
    from geophys_utils import NetCDFGridUtils
    from geophys_utils import get_netcdf_edge_points, points2convex_hull
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2)
    
    netcdf_dataset = netCDF4.Dataset(grav_netcdf_path, 'r')

    max_bytes = 500000000

    #netcdf_grid_utils = NetCDFGridUtils(netcdf_dataset)
    lats = netcdf_dataset.variables['lat'][:]
    lons = netcdf_dataset.variables['lon'][:]
    latselect = np.logical_and(lats>bbox[1],lats<bbox[3])
    lonselect = np.logical_and(lons>bbox[0],lons<bbox[2])
    data = netcdf_dataset.variables['Band1'][latselect,lonselect]
    ax[0].imshow(data,cmap ='gist_rainbow')
    
    density=gpd.read_file(density_filename,bbox=bbox)
    if(len(density)>0):
        density.plot(ax=ax[1],figsize=(7,7),column='value')
    else:
        plt.subplot(122)
        plt.text(.2,.2,"No density values\n in area")
    plt.tight_layout()
    plt.title('Gravity and Density')
    plt.show()
    
def plot_plutons(geology_filename,bbox):
    import pandas as pd
    geol=gpd.read_file(geology_filename,bbox=bbox)
    plutons=geol[ geol['ROCKTYPE1'].str.contains( 'igneous')]
    #display(plutons)
    plutons_unique=plutons.drop_duplicates(subset=['CODE'])
    #plutons["MAX_AGE_MA"] = pd.to_numeric(plutons["MAX_AGE_MA"], downcast="float")
    #plutons["MIN_AGE_MA"] = pd.to_numeric(plutons["MIN_AGE_MA"], downcast="float")
    plutons_sort_unique=plutons_unique.sort_values(by='MAX_AGE_MA')
    #display(plutons_sort_unique)
    #return

    plt.axes()
    plt.gca().set_ylim(0, 4000)
    xscale=200
    i=0
    #print(len(plutons_sort_unique))
    for ind,poly in plutons_sort_unique.iterrows():
        if(not str(poly['MIN_AGE_MA'])=='None' and not str(poly['MAX_AGE_MA'])=='None'  ):
            ave_age=float(poly['MIN_AGE_MA'])+(float(poly['MAX_AGE_MA'])-float(poly['MIN_AGE_MA']))/2.0
            #print("mafic",poly['ROCKTYPE1'])
            if("volcanic" in poly['ROCKTYPE1']):
                s='v'
            else:
                s='o'  
            if("mafic" in poly['ROCKTYPE1']):
                c='b'
            else:
                c='r'  
                
            line = plt.Line2D((i*xscale, i*xscale),(float(poly['MAX_AGE_MA']), float(poly['MIN_AGE_MA'])), lw=1.5,color='#000000')
            plt.gca().add_line(line)
            plt.plot([i*xscale], ave_age, c+s)                            
            i=i+1
    ax = plt.axis()
    plt.axis((ax[0],ax[1],ax[3],ax[2]))
    plt.title('Igneous')
    plt.savefig('Igneous.pdf')  

    plt.show() 