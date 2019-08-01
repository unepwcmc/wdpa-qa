import sys, arcpy
from wdpa.qas import arcgis_table_to_df, find_wdpa_rows, poly_checks, INPUT_FIELDS_POLY
from wdpa.stijn import output_errors_to_excel

# load input
input_poly = sys.argv[1]
output_path = sys.argv[2]

# poly
arcpy.AddMessage('Converting to pandas dataframe')
poly_df = arcgis_table_to_df(input_poly, INPUT_FIELDS_POLY)
result = dict()

arcpy.AddMessage('--- Running QA checks ---')
for poly_check in poly_checks:
    arcpy.AddMessage('Running:' + poly_check['name'])
    # checks are not currently optimised, thus return all pids regardless
    wdpa_pid = poly_check['func'](poly_df, True)

    if wdpa_pid.size > 0:
        result[poly_check['name']] = find_wdpa_rows(poly_df, wdpa_pid)

output_errors_to_excel(result, output_path, poly_checks)
