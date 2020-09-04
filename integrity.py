import sys, arcpy
import numpy as np
import pandas as pd
from wdpa.qa import arcgis_table_to_df, find_wdpa_rows, poly_checks, INPUT_FIELDS_META, forbidden_character
from wdpa.export import output_errors_to_excel

# input
input_poly = sys.argv[1]
input_pt = sys.argv[2]
input_meta = sys.argv[3]
output_path = sys.argv[4]

#input_poly = r'C:\Users\clairev\Desktop\WDPA_May2020_Public.gdb\WDPA_poly_May_2020'
#input_pt = r'C:\Users\clairev\Desktop\WDPA_May2020_Public.gdb\WDPA_point_May_2020'
#input_meta = r'C:\Users\clairev\Desktop\WDPA_May2020_Public.gdb\WDPA_source_May2020'
#output_path = r'C:\Users\clairev\Desktop'

# poly and points
INPUT_FIELDS = ['WDPAID', 'WDPA_PID', 'METADATAID', 'NAME', 'ISO3', 'DESIG']

# make dfs
df_poly = arcgis_table_to_df(input_poly, INPUT_FIELDS)
df_pt = arcgis_table_to_df(input_pt, INPUT_FIELDS)
df_meta = arcgis_table_to_df(input_meta, INPUT_FIELDS_META)

# Load the list of WDPA_PIDs to make sure have not been removed from WDPA in the update.
ID_not_to_delete = pd.read_csv('//gis-rs/gis/f03_centre_initiatives/Protected_Planet_Initiative/Updating_WDPA_and_WD_OECM/Preparing WDPA_WDOECM Release/IDs_not_to_delete/IDs_not_to_delete.csv')

result = dict()
# check duplicate WDPAID and WDPA_PID across point and polygon
result['overlap_wdpaid'] = df_poly[df_poly['WDPAID'].isin(
    np.intersect1d(df_poly['WDPAID'].values, df_pt['WDPAID'].values))]

result['overlap_wdpa_pid'] = df_poly[df_poly['WDPA_PID'].isin(
    np.intersect1d(df_poly['WDPA_PID'].values, df_pt['WDPA_PID'].values))]

# check matching metadata ID in source table
indata_meta = np.union1d(df_poly['METADATAID'].values, df_pt['METADATAID'].values)
inref_meta = df_meta['METADATAID'].values

result['metaid_only_in_data'] = df_poly[df_poly['METADATAID'].isin(
    np.setdiff1d(indata_meta, inref_meta))]

result['metaid_only_in_metadata'] =  df_meta[df_meta['METADATAID'].isin(
    np.setdiff1d(inref_meta, indata_meta))]

# check "IDs not to delete" remain in WDPA
poly_point = np.union1d(df_poly['WDPAID'].values, df_pt['WDPAID'].values)
ID_to_keep = ID_not_to_delete['WDPA_PID'].values

result['ivd_WDPA_PID_deleted'] = ID_not_to_delete[ID_not_to_delete['WDPA_PID'].isin(
    np.setdiff1d(ID_not_to_delete['WDPA_PID'].values, poly_point))]

##################################################
#check for forbidden characters in source table  #
##################################################
# ## Factory Function ##
# def forbidden_character(df_meta, check_field, return_metaid=False):
#     # Import regular expression package and the forbidden characters
#     forbidden_characters = ['<','>','?','*','\r','\n']
#     forbidden_characters_esc = [re.escape(s) for s in forbidden_characters]
#
#     pattern = '|'.join(forbidden_characters_esc)
#
#     # Obtain the WDPA_PIDs with forbidden characters
#     # remove those with nas
#     df_meta = df_meta.dropna()
#     invalid_metaid = df_meta[df_meta[check_field].str.contains(pattern, case=False)]['METADATAID'].values
#
#     if return_metaid:
#         return invalid_metaid
#
#     return len(invalid_metaid) > 0
# #check in DATA_TITLE
# def forbidden_character_data_title(df_meta, return_metaid=False):
#     check_field = 'DATA_TITLE'
#     result['forbidden_character_data_title'] = return forbidden_character(df_meta, check_field, return_metaid)
# #check in RESP_PARTY
# def forbidden_character_resp_party(df_meta, return_metaid=False):
#     check_field = 'RESP_PARTY'
#     result['forbidden_character_resp_party'] = return forbidden_character(df_meta, check_field, return_metaid)
# #check in VERIFIER
# def forbidden_character_verifier(df_meta, return_metaid=False):
#     check_field = 'VERIFIER'
#     result['forbidden_character_verifier'] = return forbidden_character(df_meta, check_field, return_metaid)
# #check in CHAR_SET
# def forbidden_character_char_set(df_meta, return_metaid=False):
#     check_field = 'CHAR_SET'
#     result['forbidden_character_char_set'] = return forbidden_character(df_meta, check_field, return_metaid)
# #check in LINEAGE
# def forbidden_character_lineage(df_meta, return_metaid=False):
#     check_field = 'LINEAGE'
#     result['forbidden_character_lineage'] = return forbidden_character(df_meta, check_field, return_metaid)
# #check in CITATION
# def forbidden_character_citation(df_meta, return_metaid=False):
#     check_field = 'CITATION'
#     result['forbidden_character_citation'] = return forbidden_character(df_meta, check_field, return_metaid)
# #check in DISCLAIMER
# def forbidden_character_disclaimer(df_meta, return_metaid=False):
#     check_field = 'DISCLAIMER'
#     result['forbidden_character_disclaimer'] = return forbidden_character(df_meta, check_field, return_metaid)
# #check in REF_SYSTEM
# def forbidden_character_ref_system(df_meta, return_metaid=False):
#     check_field = 'REF_SYSTEM'
#     result['forbidden_character_ref_system'] = return forbidden_character(df_meta, check_field, return_metaid)

output_errors_to_excel(result, output_path, [{'name': key} for key in result.keys()],'integrity', 'meta')
