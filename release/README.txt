CREATING THE TOOL

1) Copy updated 'wdpa' folder into this folder
2) Create new toolbox in ArcGIS and save it to this folder
3) Zip the above two files together into the portable toolbox. 

See Wiki for more details on how to create the tool. 


https://github.com/Yichuans/wdpa-qa/wiki#for-developers
Creating the ArcGIS toolbox
1. Open a blank project in ArcGIS Pro.
2. On the ribbon, select Insert --> Toolbox --> New Toolbox.
3. Specify the name of the toolbox.
4. Once created, the toolbox will be present in the Catalog, under Toolboxes.
5. Right-click the created toolbox --> New --> Script.
6. In tab 'General', specify the name and label, then point to a script containing the main functions (poly.py or point.py).
7. Ensure boxes 'Import script' and 'Store tool with relative path' are checked.
8. In tab 'Parameters', add two parameters to specify the input type and the output folder.
- Label: WDPA polygon, Data Type: Feature Class, Type: Required, Direction: Input.
- Label: Output folder, Data Type: Folder, Type: Required, Direction: Input.
9. Repeat steps 5-8 for another script if required (i.e. poly.py or point.py).
10. Finally, add the toolbox and the WDPA folder (containing the QA scripts) collectively in a single compressed file. Subsequently, wherever this folder is uncompressed, the toolbox will be in the same folder as the WDPA QA scripts. This makes RAMBO portable.