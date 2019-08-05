# Load packages and modules
import sys, arcpy
from wdpa.qas import arcgis_table_to_df, find_wdpa_rows, poly_checks, INPUT_FIELDS_POLY, invalid_data_import
from wdpa.stijn import output_errors_to_excel

# Load input
input_poly = sys.argv[1]
output_path = sys.argv[2]

# Let us welcome our guest of honour
arcpy.AddMessage('\nAll hail the WDPA\n')

# Convert Polygon table to pandas DataFrame
arcpy.AddMessage('Converting to pandas DataFrame')
poly_df = arcgis_table_to_df(input_poly, INPUT_FIELDS_POLY)
result = dict()

# Verify whether data import is correct. If not, exit the programme.
if invalid_data_import(poly_df, INPUT_FIELDS_POLY):
        arcpy.AddMessage('ERROR: the list of fields in the Point table is incorrect')
        sys.exit()

# Run the checks
arcpy.AddMessage('--- Running QA checks ---')
for poly_check in poly_checks: # poly_checks is a dictionary with checks' descriptive names and function names
    arcpy.AddMessage('Running:' + poly_check['name'])
    # checks are not currently optimised, thus return all pids regardless
    wdpa_pid = poly_check['func'](poly_df, True)

    # For each check, obtain the rows that contain errors
    if wdpa_pid.size > 0:
        result[poly_check['name']] = find_wdpa_rows(poly_df, wdpa_pid)

output_errors_to_excel(result, output_path, poly_checks, 'poly')
arcpy.AddMessage('The QA checks have finished. \n\nWritten by Stijn den Haan and Yichuan Shi\nAugust 2019')