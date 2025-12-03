# RAMBO: a Quality Assurance Tool for the World Database on Protected Areas

## Description

RAMBO is a quality assurance tool developed for the World Database on Protected Areas, and can be run in ArcGIS Pro. The quality assurance (QA) scripts have been written in Python.

*Input*: Feature class attribute table that conforms to the official WDPA format. Currently, polygon and point feature class attribute tables are the only allowed inputs.

*Output*: Excel Workbook with identified errors.

**For in-depth usage and a description of Quality Assurance checks included, please check the [Wiki page](https://github.com/Yichuans/wdpa-qa/wiki) on this GitHub repository.**

---

## Installation requirements

- ArcGIS Pro `2.4` or later 

Required, but included in ArcGIS Pro `2.4`: 

- Python `3.6.8` or later (included in ArcGIS Pro `2.4`)
- Python packages required (versions stated below or later):
	- pandas `0.24.2`
	- numpy `1.16.2`
	- openpyxl `2.6.1`

Note: installing Anaconda is not required. Refrain from using any other Conda installation than the one that is installed by ArcGIS Pro by default.

## Quick start

1. Download RAMBO (the WDPA QA tool), from this GitHub repository.
2. Unzip the file in a folder of your choosing.
3. Open a (non-empty / empty) project in ArcGIS Pro.
4. On the ribbon, select Insert --> Toolbox --> Add Toolbox.
5. Go to the folder where you unzipped RAMBO, select the `.tbx` file (with red icon), and press 'OK'.
6. Open the Catalog pane --> Toolboxes --> The WDPA QA toolbox should now be visible.
7. Expand the toolbox, so that the embedded scripts become visible.
8. Right-click the script to run (e.g. for polygons or points, or with OECMS included), click Open, and specify the input table (feature class attribute table) to be checked, and the output directory.
9. Click Run, and click 'View Details' if you wish to see the progress.
10. The Excel output will be present in the previously specified output directory.
11. If you encounter errors, please refer to the Troubleshooting section in the Wiki.

## Notes

The scripts "Polygons" & "Points" can be run on either the WDPA, the OECM, or the combined datasets. However, a few tests will be missed.
The scripts "OECM_Polygons" and "OECM_Points" include additional tests that REQUIRE the additional fields just found in the OECM or the WDPA_WDOECM datasets. These will fail if tried on the stand alone WDPA version!

Please refrain from committing directly to the `master` branch. Instead, create a different branch containing edits and submit a pull request. 

```bash
git checkout -b {your branch} {base branch}
```

Run tests with

```bash
python -m unittest
```

##ideas

Please add new ideas, as well as problems, into the Sharepoint doc. Ask someone from the team if you don't have access to it

## Update 2025

In December 2025 the qa script was updated to allow the tool to be used with data in the combined WDPCA in the new schema and also align with checks in the Data Management Portal (DMP). 
The DMP applies the same checks when uploading data, but this toolbox is still used by people following data updates to ensure compliance. A summary of changes is below.

- Changed field names
	- WDPAID -> SITE_ID; WDPA_PID -> SITE_PID; PA_DEF -> SITE_TYPE; MARINE -> REALM; NAME -> NAME_ENG; ORIG_NAME -> NAME; PARENT_ISO3 -> PRNT_ISO3
- Changed coding associated with changed field names
	- SITE_TYPE changed to string with specified allowed values
	- REALM changed to string with specified allowed values
- Removed SUB_LOC field and associated checks
- Added new fields and associated checks for allowed values
	- GOVSUBTYPE; OWNSUBTYPE; OECM_ASMT; INLND_WTRS
- Various changes to allowed values for some fields
- Some now-irrelevant checks removed
- Added check for excessive vertices in polygons
- Changed WDPA to WDPCA in export.py

All PA and OECM checks have been combined to align with the WDPCA. The ArcGIS toolbox (.atbx) now has two scripts: Points and Polygons. A small change was made to these execution scripts
to fix a pathing error resulting from change to .atbx from .tbx format.

Known issues not resolved:
- Status check (ivd_status) will fail for sites with 'Not Reported'. This value is only allowed for old data and should still flag for new data.
- Update to new schema included adding ';' between values in INT_CRIT. This was accounted for, however if values are not in numerical order, the check will fail (ivd_int_crit).
- GOVSUBTYPE does not allow for ';' between values when multiple (ivd_govsubtype_shared / notshared)

## Credits

Original author: Stijn den Haan

Supervisor: Yichuan Shi

Further development: Claire Vincent (2020), Sara Pruckner, Kelsey Green (2025)
---
