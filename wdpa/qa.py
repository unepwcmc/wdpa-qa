###################################################################################
#### RAMBO: a Quality Assurance Tool for the World Database on Protected Areas ####
#### Python script containing all quality assurance checks for the WDPA        ####
###################################################################################

'''
Author: Stijn den Haan
Supervisor: Yichuan Shi
Bioinformatics internship • UNEP-WCMC • 10 June - 9 August 2019

This Python script contains all quality assurance checks for the WDPA that are part of RAQTOW.
These checks are subsequently called by the 'main' scripts poly.py and point.py,
to execute the checks on the WDPA feature class attribute table provided.

## Definitions ##

**Offending fields** are WDPA fields (columns) that contain values that do not adhere to the rules set in the WDPA manual or
do not adhere to general logical rules, e.g. the marine area of the protected area being larger than the total protected area.
- Offending fields are subdivided into several types:
    - *Duplicate*: records holding exactly the same values for all fields. Notably, the WDPA_PID field should not contain duplicates.
    - *Inconsistent*: multiple records (rows) about the same protected area (same WDPAID) contain conflicting field information
        - Example: records with the same `WDPAID` have different values present in field `NAME`, e.g. 'De Veluwe' vs 'De VeLUwe'.
    - *Invalid*: a record has an incorrect value for a particular field where only a particular set of values is allowed.
        - Example: `DESIG_TYPE` = 'Individual' while only 'National', 'International', and 'Regional' are allowed values for this field.
    - *Area invalid*: a record has an incorrect value for one or several area fields.
        - Example: `GIS_M_AREA` is larger than `GIS_AREA`.
    - *Forbidden character*: a record contains a field that has a forbidden character. These can affect downstream analyses on the WDPA.
        - Example: asterisk ('*') present in `NAME`.
    - *NaN values*: a record contains a field that is NA, NaN, or None, which can be the result of e.g. division by zero.

In this document, we use:
- **field** to refer to a column of the database;
    - Example: `ISO3`
- **value** to refer to each individual entry present in a field - i.e. the intersection of the field and row.
    - Example: 12345 present in field `WDPAID` on row 12
'''

###########################################
##### 0. Load packages and WDPA fields ####
###########################################

#### Load packages ####

import numpy as np
import pandas as pd
import arcpy
import datetime
import os
import re

#### Load fields present in the WDPA tables ####

# Polygon data

INPUT_FIELDS_POLY = ['WDPAID', 'WDPA_PID', 'PA_DEF', 'NAME', 'ORIG_NAME', 'DESIG',
                     'DESIG_ENG', 'DESIG_TYPE', 'IUCN_CAT', 'INT_CRIT', 'MARINE', 'REP_M_AREA',
                     'GIS_M_AREA', 'REP_AREA', 'GIS_AREA', 'NO_TAKE', 'NO_TK_AREA', 'STATUS', 'STATUS_YR',
                     'GOV_TYPE', 'OWN_TYPE', 'MANG_AUTH', 'MANG_PLAN', 'VERIF', 'METADATAID', 'SUB_LOC',
                     'PARENT_ISO3', 'ISO3', ]

# Point data

INPUT_FIELDS_PT = ['WDPAID', 'WDPA_PID', 'PA_DEF', 'NAME', 'ORIG_NAME', 'DESIG',
                      'DESIG_ENG', 'DESIG_TYPE', 'IUCN_CAT', 'INT_CRIT', 'MARINE', 'REP_M_AREA',
                      'REP_AREA', 'NO_TAKE', 'NO_TK_AREA', 'STATUS', 'STATUS_YR', 'GOV_TYPE',
                      'OWN_TYPE', 'MANG_AUTH', 'MANG_PLAN', 'VERIF', 'METADATAID', 'SUB_LOC',
                      'PARENT_ISO3', 'ISO3', ]

# Source Table

INPUT_FIELDS_META = ['METADATAID','DATA_TITLE','RESP_PARTY','VERIFIER','YEAR',
                       'UPDATE_YR', 'LANGUAGE','CHAR_SET','REF_SYSTEM', 'SCALE',
                       'LINEAGE', 'CITATION','DISCLAIMER', ]



#####################################################
#### 1. Convert ArcGIS table to pandas DataFrame ####
#####################################################

# Use this for the Polygons, Points, and the Source Table

# Source: https://gist.github.com/d-wasserman/e9c98be1d0caebc2935afecf0ba239a0
def arcgis_table_to_df(in_fc, input_fields, query=''):
    '''
    Function will convert an arcgis table into a pandas DataFrame with an OBJECTID index, and the selected
    input fields using an arcpy.da.SearchCursor.
    For in_fc, specify the name of the geodatabase (.gdb) and feature class attribute table

    ## Arguments ##
    in_fc -- feature class attribute table - inside geodatabase - to import.
             Specify: <nameOfGeodatabase>/<nameOfFeatureClassAttributeTable>
    input_fields -- list of all fields that must be imported from the dataset
    query -- optional where_clause of arcpy.da.SearchCursor. Leave default for normal usage.

    ## Example ##
    arcgis_table_to_df(in_fc='WDPA_Jun2019_Public.gdb/WDPA_Jun2019_errortest',
    input_fields=input_fields_poly,
    query='')
    '''

    OIDFieldName = arcpy.Describe(in_fc).OIDFieldName # obtain OBJECTID field.
    final_fields = [OIDFieldName] + input_fields # Make a list of all fields that need to be extracted
    data = [row for row in arcpy.da.SearchCursor(in_fc,final_fields,where_clause=query)] # for all fields, obtain all rows
    fc_dataframe = pd.DataFrame(data,columns=final_fields) # Put data into pandas DataFrame
    fc_dataframe = fc_dataframe.set_index(OIDFieldName,drop=True) # set OBJECTID as index, but no longer use it as column
    fc_dataframe.replace('', np.nan, inplace=True) # set '' to np.nan

    return fc_dataframe


#########################################
##### 1.1 Obtain allowed ISO3 values ####
#########################################

# Download from GitHub and store in a pandas DataFrame
column_with_iso3 = ['alpha-3']
url = 'https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv'
iso3_df = pd.read_csv(url, usecols = column_with_iso3)
iso3 = np.append(iso3_df['alpha-3'].values, 'ABNJ')

#######################################
#### 2. Utility & hardcoded checks ####
#######################################

'''
The utility returns a subset of the WDPA DataFrame based on a list of WDPA_PIDs provided.
The hardcoded checks are not Factory Functions that can handle different inputs. Instead,
these are specific checks that have a set of input variables that cannot change.

'''

#############################################################################
#### 2.0. Utility to extract rows from the WDPA, based on WDPA_PID input ####
#############################################################################

def find_wdpa_rows(wdpa_df, wdpa_pid):
    '''
    Return a subset of DataFrame based on wdpa_pid list

    ## Arguments ##
    wdpa_df --  wdpa DataFrame
    wdpa_pid -- a list of WDPA_PIDs
    '''

    return wdpa_df[wdpa_df['WDPA_PID'].isin(wdpa_pid)]

#######################################
#### 2.1. Find duplicate WDPA_PIDs ####
#######################################

def duplicate_wdpa_pid(wdpa_df, return_pid=False):
    '''
    Return True if WDPA_PID is duplicate in the DataFrame.
    Return list of WDPA_PID, if duplicates are present
    and return_pid is set True.
    '''

    if return_pid:
        ids = wdpa_df['WDPA_PID'] # make a variable of the field to find
        return ids[ids.duplicated()].unique() # return duplicate WDPA_PIDs

    return wdpa_df['WDPA_PID'].nunique() != wdpa_df.index.size # this returns True if there are WDPA_PID duplicates

###########################################################################
#### 2.2. Invalid: MARINE designation based on GIS_AREA and GIS_M_AREA ####
###########################################################################

def area_invalid_marine(wdpa_df, return_pid=False):
    '''
    Assign a new 'MARINE' value based on GIS calculations, called marine_GIS_value
    Return True if marine_GIS_value is unequal to MARINE
    Return list of WDPA_PIDs where MARINE is invalid, if return_pid is set True
    '''

    # set min and max for 'coastal' designation (MARINE = 1)
    coast_min = 0.1
    coast_max = 0.9

    # create new column with proportion marine vs total GIS area
    wdpa_df['marine_GIS_proportion'] = wdpa_df['GIS_M_AREA'] / wdpa_df['GIS_AREA']

    def assign_marine_gis_value(wdpa_df):
        if wdpa_df['marine_GIS_proportion'] <= coast_min:
            return '0'
        elif coast_min < wdpa_df['marine_GIS_proportion'] < coast_max:
            return '1'
        elif wdpa_df['marine_GIS_proportion'] >= coast_max:
            return '2'

    # calculate the marine_value
    wdpa_df['marine_GIS_value'] = wdpa_df.apply(assign_marine_gis_value, axis=1)

    # find invalid WDPA_PIDs
    invalid_wdpa_pid = wdpa_df[wdpa_df['marine_GIS_value'] != wdpa_df['MARINE']]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

############################################
#### 2.3. Invalid: GIS_AREA >> REP_AREA ####
############################################

def area_invalid_too_large_gis(wdpa_df, return_pid=False):
    '''
    Return True if GIS_AREA is too large compared to REP_AREA - based on thresholds specified below.
    Return list of WDPA_PIDs where GIS_AREA is too large compared to REP_AREA, if return_pid=True
    '''

    # Set maximum allowed absolute difference between GIS_AREA and REP_AREA (in km²)
    MAX_ALLOWED_SIZE_DIFF_KM2 = 50

    # Create two Series:
    # One to calculate the mean and stdev without outliers
    # One to use as index, to find WDPA_PIDs with a too large GIS_AREA

    # Compare GIS_AREA to REP_AREA, replace outliers with NaN, then obtain mean and stdev
    # Settings
    calc =      (wdpa_df['REP_AREA'] + wdpa_df['GIS_AREA']) / wdpa_df['REP_AREA']
    condition = [calc > 100,
                calc < 0]
    choice =    [np.nan,np.nan]

    # Produce column without outliers
    relative_size_stats = pd.Series(
        np.select(condition, choice, default = calc))

    # Calculate the maximum allowed values for relative_size using mean and stdev
    max_gis = relative_size_stats.mean() + (2*relative_size_stats.std())

    # Series: compare REP_AREA to GIS_AREA
    relative_size = pd.Series((wdpa_df['REP_AREA'] + wdpa_df['GIS_AREA']) / wdpa_df['REP_AREA'])

    # Find the rows with an incorrect GIS_AREA
    invalid_wdpa_pid= wdpa_df[(relative_size > max_gis) & (abs(wdpa_df['GIS_AREA']-wdpa_df['REP_AREA']) > MAX_ALLOWED_SIZE_DIFF_KM2)]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

############################################
#### 2.4. Invalid: REP_AREA >> GIS_AREA ####
############################################

def area_invalid_too_large_rep(wdpa_df, return_pid=False):
    '''
    Return True if REP_AREA is too large compared to GIS_AREA - based on thresholds specified below.
    Return list of WDPA_PIDs where REP_AREA is too large compared to GIS_AREA, if return_pid=True
    '''

    # Set maximum allowed absolute difference between GIS_AREA and REP_AREA (in km²)
    MAX_ALLOWED_SIZE_DIFF_KM2 = 50

    # Create two Series:
    # One to calculate the mean and stdev without outliers
    # One to use as index, to find WDPA_PIDs with a too large REP_AREA

    # Compare GIS_AREA to REP_AREA, replace outliers with NaN, then obtain mean and stdev
    # Settings
    calc =      (wdpa_df['REP_AREA'] + wdpa_df['GIS_AREA']) / wdpa_df['GIS_AREA']
    condition = [calc > 100,
                calc < 0]
    choice =    [np.nan,np.nan]

    # Produce Series without outliers
    relative_size_stats = pd.Series(
        np.select(condition, choice, default = calc))

    # Calculate the maximum and minimum allowed values for relative_size using mean and stdev
    max_rep = relative_size_stats.mean() + (2*relative_size_stats.std())

    # Series: compare REP_AREA to GIS_AREA
    relative_size = pd.Series((wdpa_df['REP_AREA'] + wdpa_df['GIS_AREA']) / wdpa_df['GIS_AREA'])

    # Find the rows with an incorrect REP_AREA
    invalid_wdpa_pid= wdpa_df[(relative_size > max_rep) & (abs(wdpa_df['REP_AREA']-wdpa_df['GIS_AREA']) > MAX_ALLOWED_SIZE_DIFF_KM2)]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

################################################
#### 2.5. Invalid: GIS_M_AREA >> REP_M_AREA ####
################################################

def area_invalid_too_large_gis_m(wdpa_df, return_pid=False):
    '''
    Return True if GIS_M_AREA is too large compared to REP_M_AREA - based on thresholds specified below.
    Return list of WDPA_PIDs where GIS_M_AREA is too large compared to REP_M_AREA, if return_pid=True
    '''

    # Set maximum allowed absolute difference between GIS_M_AREA and REP_M_AREA (in km²)
    MAX_ALLOWED_SIZE_DIFF_KM2 = 50

    # Create two Series:
    # One to calculate the mean and stdev without outliers
    # One to use as index, to find WDPA_PIDs with a too large GIS_M_AREA

    # Compare GIS_M_AREA to REP_M_AREA, replace outliers with NaN, then obtain mean and stdev
    # Settings
    calc =      (wdpa_df['REP_M_AREA'] + wdpa_df['GIS_M_AREA']) / wdpa_df['REP_M_AREA']
    condition = [calc > 100,
                calc < 0]
    choice =    [np.nan,np.nan]

    # Produce column without outliers
    relative_size_stats = pd.Series(
        np.select(condition, choice, default = calc))

    # Calculate the maximum and minimum allowed values for relative_size using mean and stdev
    max_gis = relative_size_stats.mean() + (2*relative_size_stats.std())

    # Series: compare REP_M_AREA to GIS_M_AREA
    relative_size = pd.Series((wdpa_df['REP_M_AREA'] + wdpa_df['GIS_M_AREA']) / wdpa_df['REP_M_AREA'])

    # Find the rows with an incorrect GIS_M_AREA
    invalid_wdpa_pid= wdpa_df[(relative_size > max_gis) & (abs(wdpa_df['GIS_M_AREA']-wdpa_df['REP_M_AREA']) > MAX_ALLOWED_SIZE_DIFF_KM2)]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

################################################
#### 2.6. Invalid: REP_M_AREA >> GIS_M_AREA ####
################################################

def area_invalid_too_large_rep_m(wdpa_df, return_pid=False):
    '''
    Return True if REP_M_AREA is too large compared to GIS_M_AREA - based on thresholds specified below.
    Return list of WDPA_PIDs where REP_M_AREA is too large compared to GIS_M_AREA, if return_pid=True
    '''

    # Set maximum allowed absolute difference between GIS_M_AREA and REP_M_AREA (in km²)
    MAX_ALLOWED_SIZE_DIFF_KM2 = 50

    # Create two Series:
    # One to calculate the mean and stdev without outliers
    # One to use as index, to find WDPA_PIDs with a too large REP_M_AREA

    # Compare GIS_M_AREA to REP_M_AREA, replace outliers with NaN, then obtain mean and stdev
    # Settings
    calc =      (wdpa_df['REP_M_AREA'] + wdpa_df['GIS_M_AREA']) / wdpa_df['GIS_M_AREA']
    condition = [calc > 100,
                calc < 0]
    choice =    [np.nan,np.nan]

    # Produce column without outliers
    relative_size_stats = pd.Series(
        np.select(condition, choice, default = calc))

    # Calculate the maximum and minimum allowed values for relative_size using mean and stdev
    max_rep = relative_size_stats.mean() + (2*relative_size_stats.std())

    # Series: compare REP_M_AREA to GIS_M_AREA
    relative_size = pd.Series((wdpa_df['REP_M_AREA'] + wdpa_df['GIS_M_AREA']) / wdpa_df['GIS_M_AREA'])

    # Find the rows with an incorrect REP_M_AREA
    invalid_wdpa_pid= wdpa_df[(relative_size > max_rep) & (abs(wdpa_df['REP_M_AREA']-wdpa_df['GIS_M_AREA']) > MAX_ALLOWED_SIZE_DIFF_KM2)]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

#######################################################
#### 2.7. Invalid: GIS_AREA <= 0.0001 km² (100 m²) ####
#######################################################

def area_invalid_gis_area(wdpa_df, return_pid=False):
    '''
    Return True if GIS_AREA is smaller than 0.0001 km²
    Return list of WDPA_PIDs where GIS_AREA is smaller than 0.0001 km², if return_pid=True
    '''

    # Arguments
    size_threshold = 0.0001
    field_gis_area = 'GIS_AREA'

    # Find invalid WDPA_PIDs
    invalid_wdpa_pid = wdpa_df[wdpa_df[field_gis_area] <= size_threshold]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

#######################################################
#### 2.8. Invalid: REP_AREA <= 0.0001 km² (100 m²) ####
#######################################################

def area_invalid_rep_area(wdpa_df, return_pid=False):
    '''
    Return True if REP_AREA is smaller than 0.0001 km²
    Return list of WDPA_PIDs where REP_AREA is smaller than 0.0001 km², if return_pid=True
    '''

    # Arguments
    size_threshold = 0.0001
    field_rep_area = 'REP_AREA'

    # Find invalid WDPA_PIDs
    invalid_wdpa_pid = wdpa_df[wdpa_df[field_rep_area] <= size_threshold]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

#################################################
#### 2.8.a Invalid: REP_AREA >= 500,000 km²  ####
#################################################

def area_invalid_big_rep_area(wdpa_df, return_pid=False):
    '''
    Return True if REP_AREA is larger than 500,000 km²
    Return list of WDPA_PIDs where REP_AREA is larger than 500,000 km², if return_pid=True
    '''

    # Arguments
    size_threshold = 500000
    field_rep_area = 'REP_AREA'

    # Find invalid WDPA_PIDs
    invalid_wdpa_pid = wdpa_df[wdpa_df[field_rep_area] >= size_threshold]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

############################################################
#### 2.9. Invalid: REP_M_AREA <= 0 when MARINE = 1 or 2 ####
############################################################

def area_invalid_rep_m_area_marine12(wdpa_df, return_pid=False):
    '''
    Return True if REP_M_AREA is smaller than or equal to 0 while MARINE = 1 or 2
    Return list of WDPA_PIDs where REP_M_AREA is invalid, if return_pid=True
    '''

    # Arguments
    field = 'REP_M_AREA'
    field_allowed_values = 0
    condition_field = 'MARINE'
    condition_crit = ['1','2']

    # Find invalid WDPA_PIDs
    invalid_wdpa_pid = wdpa_df[(wdpa_df[field] <= field_allowed_values) & (wdpa_df[condition_field].isin(condition_crit))]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

##########################################################
## 2.10. Invalid: GIS_M_AREA <= 0 when MARINE = 1 or 2 ###
##########################################################

def area_invalid_gis_m_area_marine12(wdpa_df, return_pid=False):
    '''
    Return True if GIS_M_AREA is smaller than or equal to 0 while MARINE = 1 or 2
    Return list of WDPA_PIDs where GIS_M_AREA is invalid, if return_pid=True
    '''

    # Arguments
    field = 'GIS_M_AREA'
    field_allowed_values = 0
    condition_field = 'MARINE'
    condition_crit = ['1','2']

    # Find invalid WDPA_PIDs
    invalid_wdpa_pid = wdpa_df[(wdpa_df[field] <= field_allowed_values) & (wdpa_df[condition_field].isin(condition_crit))]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

########################################################
## 2.11. Invalid: NO_TAKE, NO_TK_AREA and REP_M_AREA ####
########################################################

def invalid_no_take_no_tk_area_rep_m_area(wdpa_df, return_pid=False):
    '''
    Return True if NO_TAKE = 'All' while the REP_M_AREA is unequal to NO_TK_AREA
    Return list of WDPA_PIDs where NO_TAKE is invalid, if return_pid=True
    '''

    # Select rows with NO_TAKE = 'All'
    no_take_all = wdpa_df[wdpa_df['NO_TAKE']=='All']

    # Select rows where the REP_M_AREA is unequal to NO_TK_AREA
    invalid_wdpa_pid = no_take_all[no_take_all['REP_M_AREA'] != no_take_all['NO_TK_AREA']]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

############################################################################
## 2.12. Invalid: INT_CRIT & DESIG_ENG - non-Ramsar Site, non-WHS sites ####
############################################################################

def invalid_int_crit_desig_eng_other(wdpa_df, return_pid=False):
    '''
    Return True if DESIG_ENG is something else than Ramsar Site (...)' or 'World Heritage Site (...)'
    while INT_CRIT is unequal to 'Not Applicable'. Other-than Ramsar / WHS should not contain anything
    else than 'Not Applicable'.
    Return list of WDPA_PIDs where INT_CRIT is invalid, if return_pid is set True
    '''

    # Arguments
    field = 'DESIG_ENG'
    field_allowed_values = ['Ramsar Site, Wetland of International Importance',
                            'World Heritage Site (natural or mixed)']
    condition_field = 'INT_CRIT'
    condition_crit = ['Not Applicable']

    # Find invalid WDPA_PIDs
    invalid_wdpa_pid = wdpa_df[(~wdpa_df[field].isin(field_allowed_values)) & (~wdpa_df[condition_field].isin(condition_crit))]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

#########################################################################
#### 2.13. Invalid: DESIG_ENG & IUCN_CAT - non-UNESCO, non-WHS sites ####
#########################################################################

def invalid_desig_eng_iucn_cat_other(wdpa_df, return_pid=False):
    '''
    Return True if IUCN_CAT is unequal to the allowed values
    and DESIG_ENG is unequal to 'UNESCO-MAB (...)' or 'World Heritage Site (...)'
    Return list of WDPA_PIDs where IUCN_CAT is invalid, if return_pid is set True
    '''

    # Arguments
    field = 'IUCN_CAT'
    field_allowed_values = ['Ia',
                            'Ib',
                            'II',
                            'III',
                            'IV',
                            'V',
                            'VI',
                            'Not Reported',
                            'Not Assigned']
    condition_field = 'DESIG_ENG'
    condition_crit = ['UNESCO-MAB Biosphere Reserve',
                      'World Heritage Site (natural or mixed)']

    # Find invalid WDPA_PIDs
    invalid_wdpa_pid = wdpa_df[(~wdpa_df[field].isin(field_allowed_values)) & (~wdpa_df[condition_field].isin(condition_crit))]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

#########################################################
#### 3. Find inconsistent fields for the same WDPAID ####
#########################################################

#### Factory Function ####

def inconsistent_fields_same_wdpaid(wdpa_df,
                                        check_field,
                                        return_pid=False):
    '''
    Factory Function: this generic function is to be linked to
    the family of 'inconsistent' input functions stated below. These latter
    functions are to give information on which fields to check and pull
    from the DataFrame. This function is the foundation of the others.

    This function checks the WDPA for inconsistent values and
    returns a list of WDPA_PIDs that have invalid values for the specified field(s).

    Return True if inconsistent Fields are found for rows
    sharing the same WDPAID

    Return list of WDPA_PID where inconsistencies occur, if
    return_pid is set True

    ## Arguments ##
    check_field -- string of the field to check for inconsistency

    ## Example ##
    inconsistent_fields_same_wdpaid(
        wdpa_df=wdpa_df,
        check_field="DESIG_ENG",
        return_pid=True):
    '''

    if return_pid:
        # Group by WDPAID to find duplicate WDPAIDs and count the
        # number of unique values for the field in question
        wdpaid_groups = wdpa_df.groupby(['WDPAID'])[check_field].nunique()

        # Select all WDPAID duplicates groups with >1 unique value for
        # specified field ('check_attributtes') and use their index to
        # return the WDPA_PIDs
        return wdpa_df[wdpa_df['WDPAID'].isin(wdpaid_groups[wdpaid_groups > 1].index)]['WDPA_PID'].values

    # Sum the number of times a WDPAID has more than 1 value for a field
    return (wdpa_df.groupby('WDPAID')[check_field].nunique() > 1).sum() > 0

#### Input functions ####

#################################
#### 3.1. Inconsistent NAME #####
#################################

def inconsistent_name_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'NAME'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''

    check_field = 'NAME'

    # The command below loads the factory function
    # and adds the check_field and return_pid arguments in it
    # to evaluate the wdpa_df for these arguments
    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

#####################################
#### 3.2. Inconsistent ORIG_NAME ####
#####################################

def inconsistent_orig_name_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'ORIG_NAME'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''

    check_field = 'ORIG_NAME'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

#################################
#### 3.3. Inconsistent DESIG ####
#################################

def inconsistent_desig_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'DESIG'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''

    check_field = 'DESIG'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

#####################################
#### 3.4. Inconsistent DESIG_ENG ####
#####################################

def inconsistent_desig_eng_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'DESIG_ENG'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''

    check_field = 'DESIG_ENG'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

######################################
#### 3.5. Inconsistent DESIG_TYPE ####
######################################

def inconsistent_desig_type_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'DESIG_TYPE'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''

    check_field = 'DESIG_TYPE'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)


####################################
#### 3.6. Inconsistent INT_CRIT ####
####################################

def inconsistent_int_crit_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'INT_CRIT'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''

    check_field = 'INT_CRIT'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

###################################
#### 3.7. Inconsistent NO_TAKE ####
###################################

def inconsistent_no_take_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'NO_TAKE'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'NO_TAKE'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

##################################
#### 3.8. Inconsistent STATUS ####
##################################

def inconsistent_status_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'STATUS'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'STATUS'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

#####################################
#### 3.9. Inconsistent STATUS_YR ####
#####################################

def inconsistent_status_yr_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'STATUS_YR'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'STATUS_YR'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

#####################################
#### 3.10. Inconsistent GOV_TYPE ####
#####################################

def inconsistent_gov_type_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'GOV_TYPE'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'GOV_TYPE'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

#####################################
#### 3.11. Inconsistent OWN_TYPE ####
#####################################

def inconsistent_own_type_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'OWN_TYPE'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'OWN_TYPE'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

######################################
#### 3.12. Inconsistent MANG_AUTH ####
######################################

def inconsistent_mang_auth_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'MANG_AUTH'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''

    check_field = 'MANG_AUTH'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

######################################
#### 3.13. Inconsistent MANG_PLAN ####
######################################

def inconsistent_mang_plan_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'MANG_PLAN'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'MANG_PLAN'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

##################################
#### 3.14. Inconsistent VERIF ####
##################################

def inconsistent_verif_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'VERIF'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'VERIF'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

#######################################
#### 3.15. Inconsistent METADATAID ####
#######################################

def inconsistent_metadataid_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'METADATAID'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'METADATAID'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

####################################
#### 3.16. Inconsistent SUB_LOC ####
####################################

def inconsistent_sub_loc_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'SUB_LOC'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'SUB_LOC'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

#######################################
### 3.17. Inconsistent PARENT_ISO3 ####
#######################################

def inconsistent_parent_iso3_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'PARENT_ISO3'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'PARENT_ISO3'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

#################################
#### 3.18. Inconsistent ISO3 ####
#################################


def inconsistent_iso3_same_wdpaid(wdpa_df, return_pid=False):
    '''
    This function is to capture inconsistencies in the field 'ISO3'
    for records with the same WDPAID

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing field inconsistencies
    '''
    check_field = 'ISO3'

    return inconsistent_fields_same_wdpaid(wdpa_df, check_field, return_pid)

##########################################
#### 4. Find invalid values in fields ####
##########################################

#### Factory Function ####

def invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid=False):
    '''
    Factory Function: this generic function is to be linked to
    the family of 'invalid' input functions stated below. These latter
    functions are to give information on which fields to check and pull
    from the DataFrame. This function is the foundation of the others.

    This function checks the WDPA for invalid values and returns a list of WDPA_PIDs
    that have invalid values for the specified field(s).

    Return True if invalid values are found in specified fields.

    Return list of WDPA_PIDs with invalid fields, if return_pid is set True.

    ## Arguments ##

    field                -- a string specifying the field to be checked
    field_allowed_values -- a list of expected values in each field
    condition_field      -- a list with another field on which the evaluation of
                            invalid values depends; leave "" if no condition specified
    condition_crit       -- a list of values for which the condition_field
                            needs to be evaluated; leave [] if no condition specified

    ## Example ##
    invalid_value_in_field(
        wdpa_df,
        field="DESIG_ENG",
        field_allowed_values=["Ramsar Site, Wetland of International Importance",
                              "UNESCO-MAB Biosphere Reserve",
                              "World Heritage Site (natural or mixed)"],
        condition_field="DESIG_TYPE",
        condition_crit=["International"],
        return_pid=True):
    '''

    # if condition_field and condition_crit are specified
    if condition_field != '' and condition_crit != []:
        invalid_wdpa_pid = wdpa_df[(~wdpa_df[field].isin(field_allowed_values)) & (wdpa_df[condition_field].isin(condition_crit))]['WDPA_PID'].values

    # If condition_field and condition_crit are not specified
    else:
        invalid_wdpa_pid = wdpa_df[~wdpa_df[field].isin(field_allowed_values)]['WDPA_PID'].values

    if return_pid:
        # return list with invalid WDPA_PIDs
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

### Factory Function ####

def invalid_value_in_field_isnot(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid=False):
    '''
    Factory Function: this generic function is adapted from the above factory function and is
    to be linked to some of the family of 'invalid' input functions stated below. These latter
    functions are to give information on which fields to check and pull from the DataFrame.

    This function checks the WDPA for invalid values and returns a list of WDPA_PIDs
    that have invalid values for the specified field(s). It is used is condition_crit is anything
    but as specified string of values.

    Return True if invalid values are found in specified fields.

    Return list of WDPA_PIDs with invalid fields, if return_pid is set True.

    ## Arguments ##

    field                -- a string specifying the field to be checked
    field_allowed_values -- a list of expected values in each field
    condition_field      -- a list with another field on which the evaluation of
                            invalid values depends; leave "" if no condition specified
    condition_crit       -- a list of values to be excluded from the condition_field
    '''

    # if condition_field and condition_crit are specified
    if condition_field != '' and condition_crit != []:
        invalid_wdpa_pid = wdpa_df[(~wdpa_df[field].isin(field_allowed_values)) & (~wdpa_df[condition_field].isin(condition_crit))]['WDPA_PID'].values

    # If condition_field and condition_crit are not specified
    else:
        invalid_wdpa_pid = wdpa_df[~wdpa_df[field].isin(field_allowed_values)]['WDPA_PID'].values

    if return_pid:
        # return list with invalid WDPA_PIDs
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

#### Input functions ####

#############################
#### 4.1. Invalid PA_DEF ####
#############################

def invalid_pa_def(wdpa_df, return_pid=False):
    '''
    Return True if PA_DEF not 1
    Return list of WDPA_PIDs where PA_DEF is not 1, if return_pid is set True
    '''

    field = 'PA_DEF'
    field_allowed_values = ['1'] # WDPA datatype is string
    condition_field = ''
    condition_crit = []

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

################################################
#### 4.2. Invalid DESIG_ENG - international ####
################################################

def invalid_desig_eng_international(wdpa_df, return_pid=False):
    '''
    Return True if DESIG_ENG is invalid while DESIG_TYPE is 'International'
    Return list of WDPA_PIDs where DESIG_ENG is invalid, if return_pid is set True
    '''

    field = 'DESIG_ENG'
    field_allowed_values = ['Ramsar Site, Wetland of International Importance',
                            'UNESCO-MAB Biosphere Reserve',
                            'World Heritage Site (natural or mixed)']
    condition_field = 'DESIG_TYPE'
    condition_crit = ['International']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

#################################################
#### 4.3. Invalid DESIG_TYPE - international ####
#################################################

def invalid_desig_type_international(wdpa_df, return_pid=False):
    '''
    Return True if DESIG_TYPE is unequal to 'International', while DESIG_ENG is an allowed 'International' value
    Return list of WDPA_PIDs where DESIG_TYPE is invalid, if return_pid is set True
    '''

    field = 'DESIG_TYPE'
    field_allowed_values = ['International']
    condition_field = 'DESIG_ENG'
    condition_crit = ['Ramsar Site, Wetland of International Importance',
                      'UNESCO-MAB Biosphere Reserve',
                      'World Heritage Site (natural or mixed)']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)


###########################################
#### 4.4. Invalid DESIG_ENG - regional ####
###########################################

def invalid_desig_eng_regional(wdpa_df, return_pid=False):
    '''
    Return True if DESIG_ENG is invalid while DESIG_TYPE is 'Regional'
    Return list of WDPA_PIDs where DESIG_ENG is invalid, if return_pid is set True
    '''

    field = 'DESIG_ENG'
    field_allowed_values = ['Baltic Sea Protected Area (HELCOM)',
                            'Specially Protected Area (Cartagena Convention)',
                            'Marine Protected Area (CCAMLR)',
                            'Marine Protected Area (OSPAR)',
                            'Site of Community Importance (Habitats Directive)',
                            'Special Protection Area (Birds Directive)',
                            'Specially Protected Areas of Mediterranean Importance (Barcelona Convention)']
    condition_field = 'DESIG_TYPE'
    condition_crit = ['Regional']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

###########################################
#### 4.5. Invalid DESIG_TYPE - regional ###
###########################################

def invalid_desig_type_regional(wdpa_df, return_pid=False):
    '''
    Return True if DESIG_TYPE is unequal to 'Regional' while DESIG_ENG is an allowed 'Regional' value
    Return list of WDPA_PIDs where DESIG_TYPE is invalid, if return_pid is set True
    '''

    field = 'DESIG_TYPE'
    field_allowed_values = ['Regional']
    condition_field = 'DESIG_ENG'
    condition_crit = ['Baltic Sea Protected Area (HELCOM)',
                      'Specially Protected Area (Cartagena Convention)',
                      'Marine Protected Area (CCAMLR)',
                      'Marine Protected Area (OSPAR)',
                      'Site of Community Importance (Habitats Directive)',
                      'Special Protection Area (Birds Directive)',
                      'Specially Protected Areas of Mediterranean Importance (Barcelona Convention)']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)


#################################################################################
#### 4.6. Invalid INT_CRIT & DESIG_ENG  - Ramsar Site & World Heritage Sites ####
#################################################################################

def invalid_int_crit_desig_eng_ramsar_whs(wdpa_df, return_pid=False):
    '''
    Return True if INT_CRIT is unequal to the allowed values (>1000 possible values)
    and DESIG_ENG equals 'Ramsar Site (...)' or 'World Heritage Site (...)'
    Return list of WDPA_PIDs where INT_CRIT is invalid, if return_pid is set True
    '''

    # Function to create the possible INT_CRIT combination
    def generate_combinations():
        import itertools
        collection = []
        INT_CRIT_ELEMENTS = ['(i)','(ii)','(iii)','(iv)',
                             '(v)','(vi)','(vii)','(viii)',
                             '(ix)','(x)']
        for length_combi in range(1, len(INT_CRIT_ELEMENTS)+1): # for 1 - 10 elements
            for combi in itertools.combinations(INT_CRIT_ELEMENTS, length_combi): # generate combinations
                collection.append(''.join(combi)) # append to list, remove the '' in each combination
        return collection

    # Arguments
    field = 'INT_CRIT'
    field_allowed_values_extra = ['Not Reported']
    field_allowed_values =  generate_combinations() + field_allowed_values_extra
    condition_field = 'DESIG_ENG'
    condition_crit = ['Ramsar Site, Wetland of International Importance',
                      'World Heritage Site (natural or mixed)']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

#################################
#### 4.7. Invalid DESIG_TYPE ####
#################################

def invalid_desig_type(wdpa_df, return_pid=False):
    '''
    Return True if DESIG_TYPE is not "National", "Regional", "International" or "Not Applicable"
    Return list of WDPA_PIDs where DESIG_TYPE is invalid, if return_pid is set True
    '''

    field = 'DESIG_TYPE'
    field_allowed_values = ['National',
                            'Regional',
                            'International',
                            'Not Applicable']
    condition_field = ''
    condition_crit = []

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

###############################
#### 4.8. Invalid IUCN_CAT ####
###############################

def invalid_iucn_cat(wdpa_df, return_pid=False):
    '''
    Return True if IUCN_CAT is not equal to allowed values
    Return list of WDPA_PIDs where IUCN_CAT is invalid, if return_pid is set True
    '''

    field = 'IUCN_CAT'
    field_allowed_values = ['Ia', 'Ib', 'II', 'III',
                            'IV', 'V', 'VI',
                            'Not Reported',
                            'Not Applicable',
                            'Not Assigned']
    condition_field = ''
    condition_crit = []

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

#####################################################################
#### 4.9. Invalid IUCN_CAT - UNESCO-MAB and World Heritage Sites ####
#####################################################################

def invalid_iucn_cat_unesco_whs(wdpa_df, return_pid=False):
    '''
    Return True if IUCN_CAT is unqueal to 'Not Applicable'
    and DESIG_ENG is 'UNESCO-MAB (...)' or 'World Heritage Site (...)'
    Return list of WDPA_PIDs where IUCN_CAT is invalid, if return_pid is set True
    '''

    field = 'IUCN_CAT'
    field_allowed_values = ['Not Applicable']
    condition_field = 'DESIG_ENG'
    condition_crit = ['UNESCO-MAB Biosphere Reserve',
                      'World Heritage Site (natural or mixed)']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

##############################
#### 4.10. Invalid MARINE ####
##############################

def invalid_marine(wdpa_df, return_pid=False):
    '''
    Return True if MARINE is not in [0,1,2]
    Return list of WDPA_PIDs where MARINE is invalid, if return_pid is set True
    '''

    field = 'MARINE'
    field_allowed_values = ['0','1','2']
    condition_field = ''
    condition_crit = []

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

############################################
#### 4.11. Invalid NO_TAKE & MARINE = 0 ####
############################################

def invalid_no_take_marine0(wdpa_df, return_pid=False):
    '''
    Return True if NO_TAKE is not equal to 'Not Applicable' and MARINE = 0
    Return list of WDPA_PIDs where NO_TAKE is invalid, if return_pid is set True
    '''

    field = 'NO_TAKE'
    field_allowed_values = ['Not Applicable']
    condition_field = 'MARINE'
    condition_crit = ['0']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

################################################
#### 4.12. Invalid NO_TAKE & MARINE = [1,2] ####
################################################

def invalid_no_take_marine12(wdpa_df, return_pid=False):
    '''
    Return True if NO_TAKE is not in ['All', 'Part', 'None', 'Not Reported'] while MARINE = [1, 2]
    I.e. check whether coastal and marine sites (MARINE = [1, 2]) have an invalid NO_TAKE value.
    Return list of WDPA_PIDs where NO_TAKE is invalid, if return_pid is set True
    '''

    field = 'NO_TAKE'
    field_allowed_values = ['All', 'Part', 'None', 'Not Reported']
    condition_field = 'MARINE'
    condition_crit = ['1', '2']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

###########################################
#### 4.13. Invalid NO_TK_AREA & MARINE ####
###########################################

def invalid_no_tk_area_marine0(wdpa_df, return_pid=False):
    '''
    Return True if NO_TK_AREA is unequal to 0 while MARINE = 0
    Return list of WDPA_PIDs where NO_TAKE is invalid, if return_pid is set True
    '''

    field = 'NO_TK_AREA'
    field_allowed_values = [0]
    condition_field = 'MARINE'
    condition_crit = ['0']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

############################################
#### 4.14. Invalid NO_TK_AREA & NO_TAKE ####
############################################

def invalid_no_tk_area_no_take(wdpa_df, return_pid=False):
    '''
    Return True if NO_TK_AREA is unequal to 0 while NO_TAKE = 'Not Applicable'
    Return list of WDPA_PIDs where NO_TK_AREA is invalid, if return_pid is set True
    '''

    field = 'NO_TK_AREA'
    field_allowed_values = [0]
    condition_field = 'NO_TAKE'
    condition_crit = ['Not Applicable']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

##############################
#### 4.15. Invalid STATUS ####
##############################

def invalid_status(wdpa_df, return_pid=False):
    '''
    Return True if STATUS is unequal to any of the following allowed values:
    ["Proposed", "Designated", "Established"] for all sites except 2 designations (WH & Barcelona convention)
    Return list of WDPA_PIDs where STATUS is invalid, if return_pid is set True

    Note: "Inscribed" and "Adopted" are only valid for specific DESIG_ENG.
    '''

    field = 'STATUS'
    field_allowed_values = ['Proposed', 'Designated', 'Established']
    condition_field = 'DESIG_ENG'
    condition_crit != ['World Heritage Site (natural or mixed)',
			'Specially Protected Areas of Mediterranean Importance (Barcelona Convention)']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

########################################
#### 4.15.a Invalid STATUS WH Sites ####
########################################

def invalid_status_WH(wdpa_df, return_pid=False):
    '''
    Return True if STATUS is unequal to any of the following allowed values:
    ["Proposed", "Designated", "Established", "Inscribed"] and DESIG_ENG is unqual to 'World Heritage Site (natural or mixed)'
    Return list of WDPA_PIDs where STATUS is invalid, if return_pid is set True

    Note: Not sure if Designated and Established are allowed for WH sites.
    '''

    field = 'STATUS'
    field_allowed_values = ["Proposed", "Designated", "Established", "Inscribed"]
    condition_field = 'DESIG_ENG'
    condition_crit = ['World Heritage Site (natural or mixed)']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

####################################################
#### 4.15.b Invalid STATUS Barcelona Convention ####
####################################################

def invalid_status_Barca(wdpa_df, return_pid=False):
    '''
    Return True if STATUS is unequal to any of the following allowed values:
    ["Proposed", "Designated", "Established", "Adopted"] and DESIG_ENG is unqual to 'Specially Protected Areas of Mediterranean Importance (Barcelona Convention)'
    Return list of WDPA_PIDs where STATUS is invalid, if return_pid is set True

    Note: Not sure if Designated and Established are allowed for Barcelona Convention sites.
    '''

    field = 'STATUS'
    field_allowed_values = ["Proposed", "Designated", "Established", "Adopted"]
    condition_field = 'DESIG_ENG'
    condition_crit = ['Specially Protected Areas of Mediterranean Importance (Barcelona Convention)']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

#################################
#### 4.16. Invalid STATUS_YR ####
#################################

def invalid_status_yr(wdpa_df, return_pid=False):
    '''
    Return True if STATUS_YR is unequal to 0 or any year between 1750 and the current year
    Return list of WDPA_PIDs where STATUS_YR is invalid, if return_pid is set True
    '''

    field = 'STATUS_YR'
    year = datetime.date.today().year # obtain current year
    yearArray = [0] + np.arange(1750, year + 1, 1).tolist() # make a list of all years, from 0 to current year
    field_allowed_values = [str(x) for x in yearArray] # change all integers to strings
    condition_field = ''
    condition_crit = []

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

################################
#### 4.17. Invalid GOV_TYPE ####
################################

def invalid_gov_type(wdpa_df, return_pid=False):
    '''
    Return True if GOV_TYPE is invalid
    Return list of WDPA_PIDs where GOV_TYPE is invalid, if return_pid is set True
    '''

    field = 'GOV_TYPE'
    field_allowed_values = ['Federal or national ministry or agency',
                            'Sub-national ministry or agency',
                            'Government-delegated management',
                            'Transboundary governance',
                            'Collaborative governance',
                            'Joint governance',
                            'Individual landowners',
                            'Non-profit organisations',
                            'For-profit organisations',
                            'Indigenous peoples',
                            'Local communities',
                            'Not Reported']

    condition_field = ''
    condition_crit = []

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

################################
#### 4.18. Invalid OWN_TYPE ####
################################

def invalid_own_type(wdpa_df, return_pid=False):
    '''
    Return True if OWN_TYPE is invalid
    Return list of WDPA_PIDs where OWN_TYPE is invalid, if return_pid is set True
    '''

    field = 'OWN_TYPE'
    field_allowed_values = ['State',
                            'Communal',
                            'Individual landowners',
                            'For-profit organisations',
                            'Non-profit organisations',
                            'Joint ownership',
                            'Multiple ownership',
                            'Contested',
                            'Not Reported']
    condition_field = ''
    condition_crit = []

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

#############################
#### 4.19. Invalid VERIF ####
#############################

def invalid_verif(wdpa_df, return_pid=False):
    '''
    Return True if VERIF is invalid
    Return list of WDPA_PIDs where VERIF is invalid, if return_pid is set True
    '''

    field = 'VERIF'
    field_allowed_values = ['State Verified',
                            'Expert Verified',
                            'Not Reported']
    condition_field = ''
    condition_crit = []

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

###################################
#### 4.20. Invalid PARENT_ISO3 ####
###################################
def invalid_country_codes(wdpa_df, field, return_pid=False):

    def _correct_iso3(field):
        for each in field.split(';'):
            if each in iso3:
                pass
            else:
                return False

        return True

    invalid_wdpa_pid = wdpa_df[~wdpa_df[field].apply(_correct_iso3)]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    else:
        return len(invalid_wdpa_pid) > 0

def invalid_parent_iso3(wdpa_df, return_pid=False):

    return invalid_country_codes(wdpa_df, 'PARENT_ISO3', return_pid)

############################
#### 4.21. Invalid ISO3 ####
############################

def invalid_iso3(wdpa_df, return_pid=False):

    return invalid_country_codes(wdpa_df, 'ISO3', return_pid)

###########################################
#### 4.22. Invalid STATUS & DESIG_TYPE ####
###########################################

def invalid_status_desig_type(wdpa_df, return_pid=False):
    '''
    Return True if STATUS is unequal to 'Established', while DESIG_TYPE = 'Not Applicable'
    Return list of WDPA_PIDs for which the STATUS is invalid
    '''

    field = 'STATUS'
    field_allowed_values = ['Established']
    condition_field = 'DESIG_TYPE'
    condition_crit = ['Not Applicable']

    return invalid_value_in_field(wdpa_df, field, field_allowed_values, condition_field, condition_crit, return_pid)

###############################################################
#### 5. Area invalid size: GIS or Reported area is invalid ####
###############################################################

#### Factory Function ####

def area_invalid_size(wdpa_df, field_small_area, field_large_area, return_pid=False):
    '''
    Factory Function: this generic function is to be linked to
    the family of 'area' input functions stated below. These latter
    functions are to give information on which fields to check and pull
    from the DataFrame. This function is the foundation of the others.

    This function checks the WDPA for invalid areas and returns a list of WDPA_PIDs
    that have invalid values for the specified field(s).

    Return True if the size of the small_area is invalid compared to large_area

    Return list of WDPA_PIDs where small_area is invalid compared to large_area,
    if return_pid is set True

    ## Arguments ##
    field_small_area  -- string of the field to check for size - supposedly smaller
    field_large_area  -- string of the field to check for size - supposedly larger

    ## Example ##
    area_invalid_size(
        wdpa_df,
        field_small_area="GIS_M_AREA",
        field_large_area="GIS_AREA",
        return_pid=True):
    '''

    size_threshold = 1.0001 # due to the rounding of numbers, there are many false positives without a threshold.

    if field_small_area and field_large_area:
        invalid_wdpa_pid = wdpa_df[wdpa_df[field_small_area] >
                                 (size_threshold*wdpa_df[field_large_area])]['WDPA_PID'].values

    else:
        raise Exception('ERROR: field(s) to test is (are) not specified')

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

#### Input functions ####

######################################################
#### 5.1. Area invalid: NO_TK_AREA and REP_M_AREA ####
######################################################

def area_invalid_no_tk_area_rep_m_area(wdpa_df, return_pid=False):
    '''
    Return True if NO_TK_AREA is larger than REP_M_AREA
    Return list of WDPA_PIDs where NO_TK_AREA is larger than REP_M_AREA if return_pid=True
    '''

    field_small_area = 'NO_TK_AREA'
    field_large_area = 'REP_M_AREA'

    return area_invalid_size(wdpa_df, field_small_area, field_large_area, return_pid)

######################################################
#### 5.2. Area invalid: NO_TK_AREA and GIS_M_AREA ####
######################################################

def area_invalid_no_tk_area_gis_m_area(wdpa_df, return_pid=False):
    '''
    Return True if NO_TK_AREA is larger than GIS_M_AREA
    Return list of WDPA_PIDs where NO_TK_AREA is larger than GIS_M_AREA if return_pid=True
    '''

    field_small_area = 'NO_TK_AREA'
    field_large_area = 'GIS_M_AREA'

    return area_invalid_size(wdpa_df, field_small_area, field_large_area, return_pid)

####################################################
#### 5.3. Area invalid: GIS_M_AREA and GIS_AREA ####
####################################################

def area_invalid_gis_m_area_gis_area(wdpa_df, return_pid=False):
    '''
    Return True if GIS_M_AREA is larger than GIS_AREA
    Return list of WDPA_PIDs where GIS_M_AREA is larger than GIS_AREA, if return_pid=True
    '''

    field_small_area = 'GIS_M_AREA'
    field_large_area = 'GIS_AREA'

    return area_invalid_size(wdpa_df, field_small_area, field_large_area, return_pid)

####################################################
#### 5.4. Area invalid: REP_M_AREA and REP_AREA ####
####################################################

def area_invalid_rep_m_area_rep_area(wdpa_df, return_pid=False):
    '''
    Return True if REP_M_AREA is larger than REP_AREA
    Return list of WDPA_PIDs where REP_M_AREA is larger than REP_AREA, if return_pid=True
    '''

    field_small_area = 'REP_M_AREA'
    field_large_area = 'REP_AREA'

    return area_invalid_size(wdpa_df, field_small_area, field_large_area, return_pid)

#################################
#### 6. Forbidden characters ####
#################################

#### Factory Function ####

def forbidden_character(wdpa_df, check_field, return_pid=False):
    '''
    Factory Function: this generic function is to be linked to
    the family of 'forbidden character' input functions stated below. These latter
    functions are to give information on which fields to check and pull
    from the DataFrame. This function is the foundation of the others.

    This function checks the WDPA for forbidden characters and returns a list of WDPA_PIDs
    that have invalid values for the specified field(s).

    Return True if forbidden characters (specified below) are found in the DataFrame

    Return list of WDPA_PID where forbidden characters occur, if
    return_pid is set True

    ## Arguments ##
    check_field -- string of the field to check for forbidden characters

    ## Example ##
    forbidden_character(
        wdpa_df,
        check_field="DESIG_ENG",
        return_pid=True):
    '''

    # Import regular expression package and the forbidden characters
    forbidden_characters = ['<','>','?','*','\r','\n']
    forbidden_characters_esc = [re.escape(s) for s in forbidden_characters]

    pattern = '|'.join(forbidden_characters_esc)

    # Obtain the WDPA_PIDs with forbidden characters
    # remove those with nas
    wdpa_df = wdpa_df.dropna()
    invalid_wdpa_pid = wdpa_df[wdpa_df[check_field].str.contains(pattern, case=False)]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

#### Input functions ####

#########################################
#### 6.1. Forbidden character - NAME ####
#########################################

def forbidden_character_name(wdpa_df, return_pid=False):
    '''
    Capture forbidden characters in the field 'NAME'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing forbidden characters in field 'NAME'
    '''

    check_field = 'NAME'

    return forbidden_character(wdpa_df, check_field, return_pid)

##############################################
#### 6.2. Forbidden character - ORIG_NAME ####
##############################################

def forbidden_character_orig_name(wdpa_df, return_pid=False):
    '''
    Capture forbidden characters in the field 'ORIG_NAME'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing forbidden characters in field 'ORIG_NAME'
    '''

    check_field = 'ORIG_NAME'

    return forbidden_character(wdpa_df, check_field, return_pid)

##########################################
#### 6.3. Forbidden character - DESIG ####
##########################################

def forbidden_character_desig(wdpa_df, return_pid=False):
    '''
    Capture forbidden characters in the field 'DESIG'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing forbidden characters in field 'DESIG'
    '''

    check_field = 'DESIG'

    return forbidden_character(wdpa_df, check_field, return_pid)

##############################################
#### 6.4. Forbidden character - DESIG_ENG ####
##############################################

def forbidden_character_desig_eng(wdpa_df, return_pid=False):
    '''
    Capture forbidden characters in the field 'DESIG_ENG'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing forbidden characters in field 'DESIG_ENG'
    '''

    check_field = 'DESIG_ENG'

    return forbidden_character(wdpa_df, check_field, return_pid)

##############################################
#### 6.5. Forbidden character - MANG_AUTH ####
##############################################

def forbidden_character_mang_auth(wdpa_df, return_pid=False):
    '''
    Capture forbidden characters in the field 'MANG_AUTH'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing forbidden characters in field 'MANG_AUTH'
    '''

    check_field = 'MANG_AUTH'

    return forbidden_character(wdpa_df, check_field, return_pid)

##############################################
#### 6.6. Forbidden character - MANG_PLAN ####
##############################################

def forbidden_character_mang_plan(wdpa_df, return_pid=False):
    '''
    Capture forbidden characters in the field 'MANG_PLAN'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing forbidden characters in field 'MANG_PLAN'
    '''

    check_field = 'MANG_PLAN'

    return forbidden_character(wdpa_df, check_field, return_pid)

############################################
#### 6.7. Forbidden character - SUB_LOC ####
############################################

def forbidden_character_sub_loc(wdpa_df, return_pid=False):
    '''
    Capture forbidden characters in the field 'SUB_LOC'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing forbidden characters in field 'SUB_LOC'
    '''

    check_field = 'SUB_LOC'

    return forbidden_character(wdpa_df, check_field, return_pid)

########################
#### 7. NaN present ####
########################

#### Factory Function ####

def nan_present(wdpa_df, check_field, return_pid=False):
    '''
    Factory Function: this generic function is to be linked to
    the family of 'nan_present' input functions stated below. These latter
    functions are to give information on which fields to check and pull
    from the DataFrame. This function is the foundation of the others.

    This function checks the WDPA for NaN / NA / None values and returns
    a list of WDPA_PIDs that have invalid values for the specified field(s).

    Return True if NaN / NA values are found in the DataFrame

    Return list of WDPA_PID where forbidden characters occur, if
    return_pid is set True

    ## Arguments ##
    check_field -- string of field to be checked for NaN / NA values

    ## Example ##
    na_present(
        wdpa_df,
        check_field="DESIG_ENG",
        return_pid=True):
    '''

    invalid_wdpa_pid = wdpa_df[pd.isna(wdpa_df[check_field])]['WDPA_PID'].values

    if return_pid:
        return invalid_wdpa_pid

    return len(invalid_wdpa_pid) > 0

#### Input functions ####

#################################
#### 7.1. NaN present - NAME ####
#################################

def nan_present_name(wdpa_df, return_pid=False):
    '''
    Capture NaN / NA in the field 'NAME'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing NaN / NA in field 'NAME'
    '''

    check_field = 'NAME'

    return nan_present(wdpa_df, check_field, return_pid)

######################################
#### 7.2. NaN present - ORIG_NAME ####
######################################

def nan_present_orig_name(wdpa_df, return_pid=False):
    '''
    Capture NaN / NA in the field 'ORIG_NAME'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing NaN / NA in field 'ORIG_NAME'
    '''

    check_field = 'ORIG_NAME'

    return nan_present(wdpa_df, check_field, return_pid)

##################################
#### 7.3. NaN present - DESIG ####
##################################

def nan_present_desig(wdpa_df, return_pid=False):
    '''
    Capture NaN / NA in the field 'DESIG'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing NaN / NA in field 'DESIG'
    '''

    check_field = 'DESIG'

    return nan_present(wdpa_df, check_field, return_pid)

######################################
#### 7.4. NaN present - DESIG_ENG ####
######################################

def nan_present_desig_eng(wdpa_df, return_pid=False):
    '''
    Capture NaN / NA in the field 'DESIG_ENG'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing NaN / NA in field 'DESIG_ENG'
    '''

    check_field = 'DESIG_ENG'

    return nan_present(wdpa_df, check_field, return_pid)

######################################
#### 7.5. NaN present - MANG_AUTH ####
######################################

def nan_present_mang_auth(wdpa_df, return_pid=False):
    '''
    Capture NaN / NA in the field 'MANG_AUTH'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing NaN / NA in field 'MANG_AUTH'
    '''

    check_field = 'MANG_AUTH'

    return nan_present(wdpa_df, check_field, return_pid)

######################################
#### 7.6. NaN present - MANG_PLAN ####
######################################

def nan_present_mang_plan(wdpa_df, return_pid=False):
    '''
    Capture NaN / NA in the field 'MANG_PLAN'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing NaN / NA in field 'MANG_PLAN'
    '''

    check_field = 'MANG_PLAN'

    return nan_present(wdpa_df, check_field, return_pid)

####################################
#### 7.7. NaN present - SUB_LOC ####
####################################

def nan_present_sub_loc(wdpa_df, return_pid=False):
    '''
    Capture NaN / NA in the field 'SUB_LOC'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing NaN / NA in field 'SUB_LOC'
    '''

    check_field = 'SUB_LOC'

    return nan_present(wdpa_df, check_field, return_pid)

#######################################
#### 7.8. NaN present - METADATAID ####
#######################################

def nan_present_metadataid(wdpa_df, return_pid=False):
    '''
    Capture NaN / NA in the field 'METADATAID'

    Input: WDPA in pandas DataFrame
    Output: list with WDPA_PIDs containing NaN / NA in field 'METADATAID'
    '''

    check_field = 'METADATAID'



#################################################################
#### 8. METADATAID: WDPA and Source Table (on the Wish List) ####
#################################################################

#######################################################################
#### 8.1. Invalid: METADATAID present in WDPA, not in Source Table ####
#######################################################################

# def invalid_metadataid_not_in_source_table(wdpa_df, wdpa_source, return_pid=False):
#     '''
#     Return True if METADATAID is present in the WDPA but not in the Source Table
#     Return list of WDPA_PIDs for which the METADATAID is not present in the Source Table
#     '''

#     field = 'METADATAID'

    ########## OPTIONAL ##########
    #### Remove METADATAID = 840 (Russian sites that are restricted and not in Source Table)
    #condition_crit = [840]
    # Remove METADATAID = 840 from the WDPA
    #wdpa_df_no840 = wdpa_df[wdpa_df[field[0]] != condition_crit[0]]
    #invalid_wdpa_pid = wdpa_df_no840[~wdpa_df_no840[field[0]].isin(
    #                                  wdpa_source[field[0]].values)]['WDPA_PID'].values
    ##############################

    # Find invalid WDPA_PIDs
#     invalid_wdpa_pid = wdpa_df[~wdpa_df[field].isin(
#                                 wdpa_source[field].values)]['WDPA_PID'].values

#     if return_pid:
#         return invalid_wdpa_pid

#     return invalid_wdpa_pid > 0

#######################################################################
#### 8.2. Invalid: METADATAID present in Source Table, not in WDPA ####
#### Note: output is METADATAIDs.                                  ####
#######################################################################

# def invalid_metadataid_not_in_wdpa(wdpa_df, wdpa_point, wdpa_source, return_pid=False):
#     '''
#     Return True if METADATAID is present in the Source Table but not in the Source Table
#     Return list of METADATAIDs for which the METADATAID is not present in the Source Table
#     '''

#     field = ['METADATAID']

#     # Concatenate all METADATAIDs of the WDPA point and poly tables
#     field_allowed_values = np.concatenate((wdpa_df[field[0]].values,wdpa_point[field[0]].values),axis=0)

#     ########## OPTIONAL ##########
#     # Remove METADATA = 840 (Russian sites that are restricted and not in Source Table)
#     #metadataid_wdpa = np.concatenate((wdpa_df[field[0]].values,wdpa_point[field[0]].values),axis=0)
#     #field_allowed_values = np.delete(metadataid_wdpa, np.where(metadataid_wdpa == 840), axis=0)
#     #######################

#     # Find METADATAIDs in the Source Table that are not present in the WDPA
#     invalid_metadataid = wdpa_source[~wdpa_source[field[0]].isin(field_allowed_values)]['METADATAID'].values

#     if return_pid:
#         return invalid_metadataid

#     return len(invalid_metadataid) > 0

############################################################################################
#### Below is a dictionary that holds all checks' descriptive (as displayed in Excel)   ####
#### and script function names (as displayed in this script, qa.py).                    ####
#### These checks are subsequently called by the main functions, poly.py and point.py,  ####
#### to run all checks on the WDPA input feature class attribute table.                 ####
############################################################################################

# Checks to be run for both point and polygon data
core_checks = [
{'name': 'duplicate_wdpa_pid', 'func': duplicate_wdpa_pid},
{'name': 'tiny_rep_area', 'func': area_invalid_rep_area},
{'name': 'big_rep_area', 'func': area_invalid_big_rep_area},
{'name': 'zero_rep_m_area_marine12', 'func': area_invalid_rep_m_area_marine12},
{'name': 'ivd_rep_m_area_gt_rep_area', 'func': area_invalid_rep_m_area_rep_area},
{'name': 'ivd_no_tk_area_gt_rep_m_area', 'func': area_invalid_no_tk_area_rep_m_area},
{'name': 'ivd_no_tk_area_rep_m_area', 'func': invalid_no_take_no_tk_area_rep_m_area},
{'name': 'ivd_int_crit_desig_eng_other', 'func': invalid_int_crit_desig_eng_other},
{'name': 'ivd_desig_eng_iucn_cat_other', 'func': invalid_desig_eng_iucn_cat_other},
{'name': 'dif_name_same_id', 'func': inconsistent_name_same_wdpaid},
{'name': 'dif_orig_name_same_id', 'func': inconsistent_orig_name_same_wdpaid},
{'name': 'ivd_dif_desig_same_id', 'func': inconsistent_desig_same_wdpaid},
{'name': 'ivd_dif_desig_eng_same_id', 'func': inconsistent_desig_eng_same_wdpaid},
{'name': 'dif_desig_type_same_id', 'func': inconsistent_desig_type_same_wdpaid},
{'name': 'dif_int_crit_same_id', 'func': inconsistent_int_crit_same_wdpaid},
{'name': 'dif_no_take_same_id', 'func': inconsistent_no_take_same_wdpaid},
{'name': 'dif_status_same_id', 'func': inconsistent_status_same_wdpaid},
{'name': 'dif_status_yr_same_id', 'func': inconsistent_status_yr_same_wdpaid},
{'name': 'dif_gov_type_same_id', 'func': inconsistent_gov_type_same_wdpaid},
{'name': 'dif_own_type_same_id', 'func': inconsistent_own_type_same_wdpaid},
{'name': 'dif_mang_auth_same_id', 'func': inconsistent_mang_auth_same_wdpaid},
{'name': 'dif_mang_plan_same_id', 'func': inconsistent_mang_plan_same_wdpaid},
{'name': 'ivd_dif_verif_same_id', 'func': inconsistent_verif_same_wdpaid},
{'name': 'ivd_dif_metadataid_same_id', 'func': inconsistent_metadataid_same_wdpaid},
{'name': 'ivd_dif_sub_loc_same_id', 'func': inconsistent_sub_loc_same_wdpaid},
{'name': 'ivd_dif_parent_iso3_same_id', 'func': inconsistent_parent_iso3_same_wdpaid},
{'name': 'ivd_dif_iso3_same_id', 'func': inconsistent_iso3_same_wdpaid},
{'name': 'ivd_pa_def', 'func': invalid_pa_def},
{'name': 'ivd_desig_eng_international', 'func': invalid_desig_eng_international},
{'name': 'ivd_desig_type_international', 'func': invalid_desig_type_international},
{'name': 'ivd_desig_eng_regional', 'func': invalid_desig_eng_regional},
{'name': 'ivd_desig_type_regional', 'func': invalid_desig_type_regional},
{'name': 'ivd_int_crit', 'func': invalid_int_crit_desig_eng_ramsar_whs},
{'name': 'ivd_desig_type', 'func': invalid_desig_type},
{'name': 'ivd_iucn_cat', 'func': invalid_iucn_cat},
{'name': 'ivd_iucn_cat_unesco_whs', 'func': invalid_iucn_cat_unesco_whs},
{'name': 'ivd_marine', 'func': invalid_marine},
{'name': 'check_no_take_marine0', 'func': invalid_no_take_marine0},
{'name': 'ivd_no_take_marine12', 'func': invalid_no_take_marine12},
{'name': 'check_no_tk_area_marine0', 'func': invalid_no_tk_area_marine0},
{'name': 'ivd_no_tk_area_no_take', 'func': invalid_no_tk_area_no_take},
{'name': 'ivd_status', 'func': invalid_status},
{'name': 'ivd_status_WH', 'func': invalid_status_WH},
{'name': 'ivd_status_BarcelonaConv', 'func': invalid_status_Barca},
{'name': 'ivd_status_yr', 'func': invalid_status_yr},
{'name': 'ivd_gov_type', 'func': invalid_gov_type},
{'name': 'ivd_own_type', 'func': invalid_own_type},
{'name': 'ivd_verif', 'func': invalid_verif},
{'name': 'check_parent_iso3', 'func': invalid_parent_iso3},
{'name': 'check_iso3', 'func': invalid_iso3},
{'name': 'ivd_status_desig_type', 'func': invalid_status_desig_type},
{'name': 'ivd_character_name', 'func': forbidden_character_name},
{'name': 'ivd_character_orig_name', 'func': forbidden_character_orig_name},
{'name': 'ivd_character_desig', 'func': forbidden_character_desig},
{'name': 'ivd_character_desig_eng', 'func': forbidden_character_desig_eng},
{'name': 'ivd_character_mang_auth', 'func': forbidden_character_mang_auth},
{'name': 'ivd_character_mang_plan', 'func': forbidden_character_mang_plan},
{'name': 'ivd_character_sub_loc', 'func': forbidden_character_sub_loc},
{'name': 'ivd_nan_present_name', 'func': nan_present_name},
{'name': 'ivd_nan_present_orig_name', 'func': nan_present_orig_name},
{'name': 'ivd_nan_present_desig', 'func': nan_present_desig},
{'name': 'ivd_nan_present_desig_eng', 'func': nan_present_desig_eng},
{'name': 'ivd_nan_present_mang_auth', 'func': nan_present_mang_auth},
{'name': 'ivd_nan_present_mang_plan', 'func': nan_present_mang_plan},
{'name': 'ivd_nan_present_sub_loc', 'func': nan_present_sub_loc},
{'name': 'ivd_nan_present_metadataid', 'func': nan_present_metadataid}]

# Checks to be run for polygon data only (includes GIS_AREA and/or GIS_M_AREA)
area_checks = [
{'name': 'gis_area_gt_rep_area', 'func': area_invalid_too_large_gis},
{'name': 'rep_area_gt_gis_area', 'func': area_invalid_too_large_rep},
{'name': 'gis_m_area_gt_rep_m_area', 'func': area_invalid_too_large_gis_m},
{'name': 'rep_m_area_gt_gis_m_area', 'func': area_invalid_too_large_rep_m},
{'name': 'tiny_gis_area', 'func': area_invalid_gis_area},
{'name': 'no_tk_area_gt_gis_m_area', 'func': area_invalid_no_tk_area_gis_m_area},
{'name': 'ivd_gis_m_area_gt_gis_area', 'func': area_invalid_gis_m_area_gis_area},
{'name': 'zero_gis_m_area_marine12', 'func': area_invalid_gis_m_area_marine12},
{'name': 'ivd_marine_designation', 'func': area_invalid_marine},]

# Checks for polygons
poly_checks = core_checks + area_checks

# Checks for points (area checks excluded)
pt_checks = core_checks

#######################
#### END OF SCRIPT ####
#######################
