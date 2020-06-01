# Load packages and modules
import sys, arcpy
from wdpa.qa import arcgis_table_to_df, find_wdpa_rows, pt_checks, INPUT_FIELDS_PT
from wdpa.export import output_errors_to_excel

# Load input
input_pt = sys.argv[1]
output_path = sys.argv[2]

inputfilename = os.path.basename(input_pt)

# Let us welcome our guest of honour
arcpy.AddMessage('\nAll hail the WDPA\n')

# Convert Point table to pandas DataFrame
arcpy.AddMessage('Converting to pandas DataFrame')
pt_df = arcgis_table_to_df(input_pt, INPUT_FIELDS_PT)
result = dict()

# Run the checks
arcpy.AddMessage('--- Running QA checks on Points ---')
for pt_check in pt_checks: # pt_checks is a dictionary with checks' descriptive names and function names
    arcpy.AddMessage('Running:' + pt_check['name'])
    # checks are not currently optimised, thus return all pids regardless
    wdpa_pid = pt_check['func'](pt_df, True)

    # For each check, obtain the rows that contain errors
    if wdpa_pid.size > 0:
        result[pt_check['name']] = find_wdpa_rows(pt_df, wdpa_pid)

# Write output to file
arcpy.AddMessage('Writing output to Excel')
output_errors_to_excel(result, output_path, pt_checks, inputfilename, 'point')
arcpy.AddMessage('\nThe QA checks on POINTS have finished. \n\nWritten by Stijn den Haan and Yichuan Shi\nAugust 2019')
