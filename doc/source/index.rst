


.. COCS documentation master file, created by
   sphinx-quickstart on Tue May 29 10:15:02 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========================================================
Pack/Unpack IPSL imodel diagnostics
========================================================

**PyPUMD** (Python Pack/Unpack Model Diagnostics) is a library in 
python that helps the user to aspects like:

 -Create netcdf files for collaboration groups like:

   - CCMI: by recreating hierarchical folders
   - CMIP6: by testing metadata according to CMOR tables.
   - AeroCOM: by organizing files according to their requirements.

 - Provide automatic checks of the created netcdf files 


The suite is then organized in several logical steps:

1. Extract single netcdf for each variable from the 
   aggregated original climate files

2. Ascertain if this netcdf files are CF compilant if possible.

3. 

example by CCMI. Also it check that the metadata
is compatible with the CMOR tables, and give filenames structured
according with CMIP6/AeroCOM activities. So it provides acollection 
of modules and functions created tounpack diagnostics from netcdf 
files of climate model simulations, as well as perform test and 
evaluations of these diagnostics. It was created for IPSL-CM simulations, 
althought it is easy extend it to use with other climate models. 
The main objectives of this tool are:

- Create single diagnostics netcdf files of climate model outputs.
- Perform quick evaluations and test of these diagnostics.
- Perform simple operations with a given set of diagnostics to create
  derived ones.
- Prepare netcdf files for specifc projets: CMIP6, AeroCom etc, that
  respect the metadata and filenames requires structure.
- The tool automatize many common tasks with common and easy to use YAML
  setting files.

Dependencies
============


It uses several python libraries like xarray_, pandas_, yaml_, json_, tabulate_,
It has been selected the cartopy library over the basemap as it will have more
long-term support.

.. _xarray: http://xarray.pydata.org
.. _pandas: http://pandas.pydata.org
.. _yaml: http://pandas.pydata.org
.. _json: http://pandas.pydata.org
.. _netCDF4: http://www.unidata.ucar.edu/software/netcdf
.. _matplotlib: https://matplotlib.org/
.. _cartopy: https://scitools.org.uk/cartopy/docs/latest/

Additional it used the cmip6-cmor tables given as json files

Documentation
=============

**Getting Started**

* :doc:`Structure`
* :doc:`How-to-use`
* :doc:`api`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting Started

   Structure
   How-to-use
   api



Get in touch
============

- CRESCENDO project `Crescendo`_.
- LSCE-Mermaid `Mermaid`_.
- Author contact: Ramiro Checa-Garcia rcheca@lsce.ipsl.fr

.. _Crescendo: http://www.crescendoproject.eu/
.. _Mermaid: https://www.lsce.ipsl.fr/en/Phocea/Vie_des_labos/Ast/ast_groupe.php?id_groupe=94

License
=======

pyanalysisIPSL will be probably available under the open source `GNU3 License`__.
but maybe it will be specific of French Institutions Open Source

__ https://www.gnu.org/licenses/gpl-3.0.en.html

History
=======

* v0.7 Sept 2018. Contributors: Ramiro Checa-Garcia
* v0.6 July 2018. Contributors: Ramiro Checa-Garcia
* v0.5 June 2018. Contributors: Ramiro Checa-Garcia

See also
========

- Stephan Hoyer and Joe Hamman's `Journal of Open Research Software paper`_ describing the xarray project.
- Mignot et Bony, `IPSL model introductory paper`_ 

.. _Journal of Open Research Software paper: http://doi.org/10.5334/jors.148
.. _IPSL model introductory paper: Mignot, J. & Bony, S. Clim Dyn (2013) 40: 2089. https://doi.org/10.1007/s00382-013-1720-1



