How to use: running in TGCC cluster
===================================


**PyPUMD** it is possible to use in other machines if the different 
libraries used are installed. The current version has been only checked
at the TGCC server, where you have to load the module of python3.
Note however that given that the TGCC 

TGCC machines
-------------

**COBALT**

The current version has been developed for cobalt users as curie doesn't
have the required libraries. Also by request of the R. Checa-Garcia
the TGCC team installed cartopy, xarray, netCDF4 and other libraries.

**CURIE**

Currently not all the required libraries are supported.

**IRENE**

It shares the same modules than cobalt so its works fine also in this
HPC sever. Like in other cases, it is needed to load the python module as

> module load python3 

before use any of the codes here described.



Cartopy-files
-------------

Normaly, cartopy first run download from several servers the mapping 
information. Given that TGCC users can not download files from extrnal
servers, they have to be installed manually by each users and cartopy
should be informed of that with the following piece of code,

.. code-block:: python

    from cartopy.io import Downloader

    # -- CARTOPY LOAD OF SHAPEFILES -----------------------------------------------------
    shaped='/ccc/cont003/home/dsm/checagar/.local/share/cartopy/shapefiles/'
    shapef=shaped+'natural_earth/physical/ne_110m_coastline.shp'
    Downloader.from_config(('shapefiles', 'natural_earth')).pre_downloaded_path_template=shapef
    Downloader.from_config(('shapefiles', 'natural_earth')).target_path_template=shapef


The above code should be present on any use of the library which aims to plots maps.


