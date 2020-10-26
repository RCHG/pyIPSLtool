


.. pyanalysisIPSL documentation master file, created by
   sphinx-quickstart on Tue May 29 10:15:02 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

================
Structure 
================

Here is described the file structure of the library

File Tree
---------

Main files and folders::

     .
    ├── calibration_dust.py
    ├── dust_analysis.py
    ├── main_irene.py
    ├── regions_RIDLEY_study.py
    ├── regions_study.py
    ├── main.py
    ├── seasonal.py
    │ 
    ├── doc
    ├── etc
    │   ├── calibration_dust
    │   │   └── comparsion_calibration_LMDZORINCA_CRESCENDOWP6_PDv2.dat   
    │   ├── emis_alt.json
    │   ├── emissions.json
    │   ├── h2o_modes.json
    │   ├── loadings.json
    │   ├── mdw_modes.json
    │   ├── opt_depth.json
    │   ├── regions
    │   │   ├── ATmask_robinson.png
    │   │   ├── CORRECTIONmask_robinson.png
    │   │   ├── dataset_ridley_corrected_Jun2018_extended_for_calibrationIPSL.txt
    │   │   ├── dataset_ridley_corrected_Jun2018.txt
    │   │   ├── dataset_ridley_corrected.txt
    │   │   ├── dataset_ridley_original.txt
    │   │   ├── KIMmask_robinson.png
    │   │   ├── RHKZmask_PlateCarree.png
    │   │   └── RHKZmask_robinson.png
    │   └── rhv_files
    │       ├── rhv145x143.Sep2018v2.nc
    │       ├── rhv145x143.Sep2018v2_updated.nc
    │       ├── rhvEC_v2.nc
    │       └── rhvEC_v2_updated.nc
    ├── examples [17 entries exceeds filelimit, not opening dir]
    ├── lib
    │   ├── aerplots.py
    │   ├── aux_forplot.py
    │   ├── dic_definitions.py
    │   ├── __init__.py
    │   ├── loops.py
    │   ├── regions_old.py
    │   └── regions.py
    ├── non_used
    ├── old
    └── out [22 entries exceeds filelimit, not opening dir]

    28 directories, 94 files

Description
-----------

Main folders::

    * /lib               -> main modules 
    * /etc               -> specific information to be used by library
    * /non_used and /old -> old scripts
    * /doc               -> documentation
    * /examples          -> similar to /non_used??
    * /out               -> figures that are part of direct outputs and
                            temporary analysis.

Files::

    * calibration_dust.py  ----> IMPORTANT it creates correction files. 
    * dust_analysis.py
    * main_irene.py
    * regions_RIDLEY_study.py -> IMPORTANT it creates calibration tables.
    * regions_study.py
    * main.py
    * seasonal.py




