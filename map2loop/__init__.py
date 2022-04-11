# if you add a variable in here it can be changed from within python code e.g.
# import map2loop
# map2loop._clut_path = newpath

__version__ = "1.2.0"

geology_loopdata = {
    'WA':
    'http://13.211.217.129:8080/geoserver/loop/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=loop:500k_geol_28350&bbox={BBOX_STR}&srs=EPSG:28350&outputFormat=shape-zip',
    'QLD':
    'http://13.211.217.129:8080/geoserver/QLD/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=QLD:qld_geol_asud&bbox={BBOX_STR}&srsName=EPSG:28355&outputFormat=shape-zip',
    'SA':
    'http://13.211.217.129:8080/geoserver/SA/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=SA:2m_surface_geology_28354_relage&bbox={BBOX_STR}&srs=EPSG:28354&outputFormat=shape-zip',
    'VIC': 
    'http://13.211.217.129:8080/geoserver/VIC/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=VIC:geolunit_250k_py_28355&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
    'NSW': 
    'http://13.211.217.129:8080/geoserver/NSW/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=NSW:ge_rockunit_lao2&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
    'ACT': '',
    'TAS':
    'http://13.211.217.129:8080/geoserver/TAS/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=TAS:geol_poly_250_28355&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',

}
structure_loopdata = {
    'WA':
    'http://13.211.217.129:8080/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=loop:waroxi_wa_28350&bbox={BBOX_STR}&srs=EPSG:28350&outputFormat=shape-zip',
    'QLD':
    'http://13.211.217.129:8080/geoserver/QLD/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=QLD:outcrops_28355&bbox={BBOX_STR}&srsName=EPSG:28355&outputFormat=shape-zip',
    'SA':
    'http://13.211.217.129:8080/geoserver/SA/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=SA:sth_flinders_28354&bbox={BBOX_STR}&srs=EPSG:28354&outputFormat=shape-zip',
    'VIC': 
    'http://13.211.217.129:8080/geoserver/VIC/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=VIC:struc_28355&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
    'NSW': 
    'http://13.211.217.129:8080/geoserver/NSW/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=NSW:lao_struct_pt&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
    'ACT': '',
    'TAS': 
    'http://13.211.217.129:8080/geoserver/TAS/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=TAS:geol_struc_250_28355&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
}
fault_loopdata = {
    'WA':
    'http://13.211.217.129:8080/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=loop:500k_faults_28350&bbox={BBOX_STR}&srs=EPSG:28350&outputFormat=shape-zip',
    'QLD':
    'http://13.211.217.129:8080/geoserver/QLD/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=QLD:qld_faults_folds_28355&bbox={BBOX_STR}&srsName=EPSG:28355&outputFormat=shape-zip',
    'SA':
    'http://13.211.217.129:8080/geoserver/SA/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=SA:2m_linear_structures&bbox={BBOX_STR}&srs=EPSG:28354&outputFormat=shape-zip',
    'VIC': 
    'http://13.211.217.129:8080/geoserver/VIC/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=VIC:geolstructure_250k_ln_28355&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
    'NSW':  
    'http://13.211.217.129:8080/geoserver/NSW/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=NSW:faults_joined_left_contains_drop_dups&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
    'ACT': '',
    'TAS': 
    'http://13.211.217.129:8080/geoserver/TAS/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=TAS:geol_line_250_28355&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
}
fold_loopdata = {
    'WA':
    'http://13.211.217.129:8080/geoserver/loop/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=loop:500k_faults_28350&bbox={BBOX_STR}&srs=EPSG:28350&outputFormat=shape-zip',
    'QLD':
    'http://13.211.217.129:8080/geoserver/QLD/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=QLD:qld_faults_folds_28355&bbox={BBOX_STR}&srsName=EPSG:28355&outputFormat=shape-zip',
    'SA':
    'http://13.211.217.129:8080/geoserver/SA/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=SA:2m_linear_structures&bbox={BBOX_STR}&srs=EPSG:28354&outputFormat=shape-zip',
    'VIC':     
    'http://13.211.217.129:8080/geoserver/VIC/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=VIC:geolstructure_250k_ln_28355&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
    'NSW': 
    'http://13.211.217.129:8080/geoserver/NSW/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=NSW:folds_lao&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
    'ACT': '',
    'TAS': 
    'http://13.211.217.129:8080/geoserver/TAS/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=TAS:geol_line_250_28355&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
}
mindep_loopdata = {
    'WA':
    'http://13.211.217.129:8080/geoserver/loop/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=loop:mindeps_2018_28350&bbox={BBOX_STR}&srs=EPSG:28350&outputFormat=shape-zip',
    'QLD':
    'http://13.211.217.129:8080/geoserver/QLD/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=qld_mindeps_28355&bbox={BBOX_STR}&srsName=EPSG:28355&outputFormat=shape-zip',
    'SA': '',
    'VIC':
    'http://13.211.217.129:8080/geoserver/VIC/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=VIC:mindeps_2018_28350&bbox={BBOX_STR}&srs=EPSG:28350&outputFormat=shape-zip',
    'NSW': 
    'http://13.211.217.129:8080/geoserver/NSW/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=NSW:nsw_mindeps&bbox={BBOX_STR}&srs=EPSG:28355&outputFormat=shape-zip',
    'ACT': '',
    'TAS': 
    'http://13.211.217.129:8080/geoserver/VIC/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=VIC:mindeps_2018_28350&bbox={BBOX_STR}&srs=EPSG:28350&outputFormat=shape-zip',
}
metafiles = {
    'WA':
    'https://gist.githubusercontent.com/yohanderose/8f843de0dde531f009a3973cbdadcc9f/raw/cad97e91a71b4e4791629a6d384e2700095eb3b6/WA_meta.hjson',
    'QLD':
    'https://gist.githubusercontent.com/yohanderose/7ad2ae1b36e4e0b3f27dff17eeae0cc2/raw/ec8d33e52d913e2df4658ca94e7eb467e35adccd/QLD_meta.hjson',
    'SA':
    'https://gist.githubusercontent.com/yohanderose/9401ce1ac72969e2e7fe78e64ca57212/raw/bb7afca4436b6a844053b0ae06978bf59495ed87/SA_meta.hjson',
    'VIC': 
    'https://gist.githubusercontent.com/markjessell/17acdbc9061530af76109b0afe519673/raw/80816c511ebad0152a826edcd51f5bea32c2b3db/VIC.hjson',
    'NSW': 
    'https://gist.githubusercontent.com/markjessell/0bcb02e6f3ae6da9f00f299b5f0a7f77/raw/09c8018e2d92c7f0051a45be4f9e19dc7cd9a30c/NSW.hjson',
    'ACT': '',
    'TAS': 
    'https://gist.githubusercontent.com/markjessell/e6521f0d0b7a29014deb04996033554f/raw/296c364b618a49654e18eba10196fa3cc2aabd2f/TAS.hjson'
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
