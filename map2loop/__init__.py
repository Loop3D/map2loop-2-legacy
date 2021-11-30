# if you add a variable in here it can be changed from within python code e.g.
# import map2loop
# map2loop._clut_path = newpath

__version__ = "1.1.75"

geology_loopdata = {
    'WA':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=loop:geol_500k&bbox=&srs=EPSG:28350',
    'QLD':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=qld_geol_dissolved_join_fix_mii_clip_wgs84&bbox=&srsName=EPSG:28355',
    'SA':
    'http://geo.loop-gis.org/geoserver/GSSA/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=GSSA:2m_surface_geology_28354_relage&bbox=&srs=EPSG:28354',
    'VIC': 
    'http://geo.loop-gis.org/geoserver/GSV/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=GSV:geolunit_250k_py_28355&bbox=&srs=EPSG:28355',
    'NSW': '',
    'ACT': '',
    'TAS': ''
}
structure_loopdata = {
    'WA':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=waroxi_wa_28350_bed&bbox=&srs=EPSG:28350',
    'QLD':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=outcrops_28355&bbox=&srsName=EPSG:28355',
    'SA':
    'http://geo.loop-gis.org/geoserver/GSSA/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=GSSA:sth_flinders_28354&bbox=&srs=EPSG:28354',
    'VIC': 
    'http://geo.loop-gis.org/geoserver/GSV/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=GSV:struc_28355&bbox=&srs=EPSG:28355',
    'NSW': '',
    'ACT': '',
    'TAS': ''
}
fault_loopdata = {
    'WA':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=linear_500k&bbox=&srs=EPSG:28350',
    'QLD':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=qld_faults_folds_28355&bbox=&srsName=EPSG:28355',
    'SA':
    'http://geo.loop-gis.org/geoserver/GSSA/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=GSSA:2m_linear_structures_28354&bbox=&srs=EPSG:28354',
    'VIC': 
    'http://geo.loop-gis.org/geoserver/GSV/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=GSV:geolstructure_250k_ln_28355&bbox=&srs=EPSG:28355',
    'NSW': '',
    'ACT': '',
    'TAS': ''
}
fold_loopdata = {
    'WA':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=linear_500k&bbox=&srs=EPSG:28350',
    'QLD':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=qld_faults_folds_28355&bbox=&srsName=EPSG:28355',
    'SA':
    'http://geo.loop-gis.org/geoserver/GSSA/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=GSSA:2m_linear_structures_28354&bbox=&srs=EPSG:28354',
    'VIC':     
    'http://geo.loop-gis.org/geoserver/GSV/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=GSV:geolstructure_250k_ln_28355&bbox=&srs=EPSG:28355',
    'NSW': '',
    'ACT': '',
    'TAS': ''
}
mindep_loopdata = {
    'WA':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=loop:mindeps_2018_28350&bbox=&srs=EPSG:28350',
    'QLD':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=qld_mindeps_28355&bbox=&srsName=EPSG:28355',
    'SA': '',
    'VIC':
    'http://geo.loop-gis.org/geoserver/loop/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=loop:mindeps_2018_28350&bbox=&srs=EPSG:28350',
    'NSW': '',
    'ACT': '',
    'TAS': ''
}
metafiles = {
    'WA':
    'https://gist.githubusercontent.com/yohanderose/8f843de0dde531f009a3973cbdadcc9f/raw/cad97e91a71b4e4791629a6d384e2700095eb3b6/WA_meta.hjson',
    'QLD':
    'https://gist.githubusercontent.com/yohanderose/7ad2ae1b36e4e0b3f27dff17eeae0cc2/raw/ec8d33e52d913e2df4658ca94e7eb467e35adccd/QLD_meta.hjson',
    'SA':
    'https://gist.githubusercontent.com/yohanderose/9401ce1ac72969e2e7fe78e64ca57212/raw/bb7afca4436b6a844053b0ae06978bf59495ed87/SA_meta.hjson',
    'VIC': './source_data/VIC.hjson',
    'NSW': '',
    'ACT': '',
    'TAS': ''
}
clut_paths = {
    'WA':
    'https://gist.githubusercontent.com/yohanderose/8f7e2d57db9086fbe1a7c651b9e25996/raw/9144994d162662ec015321699c3658a8dbf57def/WA_clut.csv',
    'QLD': '',
    'SA':
    'https://gist.githubusercontent.com/yohanderose/c7927fa6f9a2b9e983eeec6941ad9f56/raw/253b33afee74fcb1d487ef1c12093134f31093f4/SA_clut.csv',
    'VIC': '',
    'NSW': '',
    'ACT': '',
    'TAS': ''
}
