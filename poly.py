# Load packages and modules
import sys, arcpy, os
from wdpa.qa import arcgis_table_to_df, find_wdpa_rows, poly_checks, INPUT_FIELDS_POLY
from wdpa.export import output_errors_to_excel

# Load input
input_poly = sys.argv[1]
output_path = sys.argv[2]

inputfilename = os.path.basename(input_poly)

# Let us welcome our guest of honour
arcpy.AddMessage('\nAll hail the WDPA\n')

# Convert Polygon table to pandas DataFrame
arcpy.AddMessage('Converting to pandas DataFrame')
poly_df = arcgis_table_to_df(input_poly, INPUT_FIELDS_POLY)
result = dict()

# Run the checks
arcpy.AddMessage('--- Running QA checks on Polygons ---')
for poly_check in poly_checks: # poly_checks is a dictionary with checks' descriptive names and function names
    arcpy.AddMessage('Running:' + poly_check['name'])
    # checks are not currently optimised, thus return all pids regardless
    wdpa_pid = poly_check['func'](poly_df, True)

    # For each check, obtain the rows that contain errors
    if wdpa_pid.size > 0:
        result[poly_check['name']] = find_wdpa_rows(poly_df, wdpa_pid)

# Write output to file
arcpy.AddMessage('Writing output to Excel')
output_errors_to_excel(result, output_path, poly_checks, inputfilename, 'poly')
arcpy.AddMessage('\nThe QA checks on POLYGONS have finished. \n\nWritten by Stijn den Haan and Yichuan Shi\nAugust 2019')
