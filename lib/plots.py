"""
##############################################################################################

  Project   : CRESCENDO/AEROCOM
  Filename  : check.py
  Author    : Ramiro Checa-Garcia
  email     : rcheca@lsce.ipsl.fr
  Purpose   : Specific type of checks for key variables

  Revision History ----------------------------------------------------------

  Date       Author     Ref    Revision

  2018-Oct   R.Checa           First version.
  2018-Nov   R.Checa           Added the global mean at level


  TODO LIST:


##############################################################################################
"""

import xarray as xr
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt              # Main wrapper to matplotlib library
import matplotlib.colors as mplcolors        # To handle with colors in matplotlib
from   matplotlib.ticker import MaxNLocator  # to change tickers for custom ones
from   matplotlib.ticker import LogLocator   # --
import matplotlib.ticker as mticker
from   matplotlib.ticker import StrMethodFormatter

import cartopy.crs as ccrs                   # To create maps on several projections
from   cartopy.util import add_cyclic_point  # To create cyclic coordinates for maps
from   cartopy.mpl.ticker import LongitudeFormatter
from   cartopy.mpl.ticker import LatitudeFormatter
import os.path
import math


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


# -- CARTOPY LOAD OF SHAPEFILES -----------------------------------------------------
from   cartopy.io import Downloader
shaped='/ccc/cont003/home/dsm/checagar/.local/share/cartopy/shapefiles/'
shapef=shaped+'natural_earth/physical/ne_110m_coastline.shp'
Downloader.from_config(('shapefiles', 'natural_earth')).pre_downloaded_path_template=shapef
Downloader.from_config(('shapefiles', 'natural_earth')).target_path_template=shapef




def myprint(mystr, finfo=None):
    print(mystr)
    if finfo!=None:
        finfo.write(mystr+'\n')
    return


def latlon_mean(datarray_input):
    '''
    INPUT: DataArray (lat, lon)
    RETURN: floa
    '''
    def lat_average(xarr, vlat):
        lat      = xarr[vlat]
        dg_to_rd = np.pi/180.
        cleaned_data = np.ma.masked_array(xarr.values, np.isnan(xarr.values))
        lataxis = xarr.get_axis_num(vlat)
        return np.ma.average(cleaned_data, axis=lataxis, weights=np.cos(lat*dg_to_rd))

    data_lat = datarray_input.mean(dim='lon')

    array_avg = lat_average(data_lat, 'lat')  # here we use the get_axis_num of xarray
    return array_avg

def singlemap(data2D, varname, pname, swap=False, extratitle='', style='discrete', **kwargs):
    """
    Simple function to plot a 2D map based on lat-lon dataset; it is important
    here that data2D has dimensions (lat,lon) in that order to perform correctly
    the average and cyclic added point.

    How to use (minimal):
        singlemap(data2D, 'myplot.png', **dic_forplot)

    The two style values are 'discrete' or anyother that sets a symetrical lognormal weight
    for the colorbar. For discrete colorbar the important values are vmin, vmax, custom_cmap,
    for the symetrical lognormal colorbar the linthresh and linscale together with vmin and
    vmas are the relevant values.

    Args:
        data2D      (data.xarray): with coordinates [lat, lon] data values to map
        pname       (string):      figure name to save.
        swap        (bool):        if data2D is [lon,lat] introduce swap=true. Default=False
        extratitle  (string):      extra info for title (named argument). Default=''
        style       (string):      select kind of colorbar. Default='discrete'
    kwargs:
        custom_cmap (cmap):        optional colormap for plot
        clev        (list):        optional list of color levels and ticks
        linthresh   (float):       linear thershold
        linscale    (float):       linear scale
        vmin        (float):       min value colormap
        vmax        (float):       max value colormap

    Returns:
        Nothing, it saves a figure with pname
    """

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    def get_omag(value, up=True):
        if (value == 0): 
            return value
        else:
            omag = int(math.floor(math.log10(abs(value))))
            scale=10**omag
            if up==True:
                newomag = scale*((value//scale)+1.0)
            else:
                newomag = scale*((value//scale))

            print(omag, value//scale, scale, value,  newomag)
        return newomag

    lat = data2D.coords['lat'].values  #lat is now a np.array
    lon = data2D.coords['lon'].values  #lon is now a np.array

    # Like in a pure numpy based function, lon and lat axis are still needed for
    # weighted means but the position on the data2D array could get from xarrays
    # Here we control which axis is lat and which is lon for legacy purposes with
    # respect to previous versions.

    # Default settings function ---------------------------------------------------

    lonaxis=1
    lataxis=0

    p_cmap='Spectral'                                     # plot cmap
    p_linthresh=1                                         # plot linthresh
    p_linscale=1                                          # plot linscale
    p_vmin=get_omag(np.amin(data2D[varname].values), up=False) # plot vmin
    p_vmax=get_omag(np.amax(data2D[varname].values), up=True ) # plot vmax

    # p_clev=[0.0001, 0.25, 0.5, 1, 2.5, 5, 10, 20, 30] # plot color levels

    if swap==True:
        lonaxis=0; lataxis=1

    if kwargs is not None:
        for key, value in kwargs.items():
            if key=='custom_cmap':
                p_cmap=kwargs['custom_cmap']
            if key=='linthresh':
                p_linthresh=kwargs['linthresh']
            if key=='linscale':
                p_linscale=kwargs['linscale']
            if key=='vmin':
                p_vmin=kwargs['vmin']
            if key=='vmax':
                p_vmax=kwargs['vmax']
            if key=='clev':
                p_clev=kwargs['clev']

    if p_vmin > 0 and p_vmin <= 0.001:
        n_vmin=0.0
    else:
        n_vmin=p_vmin

    bounds_figure = np.linspace(n_vmin, p_vmax, 11)

    norm_discrete = mplcolors.BoundaryNorm(boundaries=bounds_figure, ncolors=256)

    # xarrays does not offer weighted means

    data2Davg = latlon_mean(data2D[varname])

    data2Dcyclic, loncyclic = add_cyclic_point(data2D[varname].values, coord=lon, axis=lonaxis)

    normali = mplcolors.SymLogNorm(linthresh=p_linthresh,
                                   linscale=p_linscale,
                                   vmin=p_vmin,
                                   vmax=p_vmax)
    fig = plt.figure(figsize=(10, 5))                  # Figure size
    ax = fig.add_subplot(1, 1, 1,
                         projection=ccrs.Robinson())   # Our final projection 

    if style=='discrete':
        cb = ax.pcolormesh(loncyclic, lat, data2Dcyclic,
                           transform=ccrs.PlateCarree(),   # Dataset are in latlon
                           norm=norm_discrete,             # Normalization of cmap vs data
                           cmap=p_cmap)                    # cmap info.
    else:
        cb = ax.pcolormesh(loncyclic, lat, data2Dcyclic,
                           transform=ccrs.PlateCarree(),   # Dataset are in latlon
                           norm=normali,                   # Normalization of cmap vs data
                           cmap=p_cmap)                    # cmap info.

    ax.coastlines()
    ax.set_global()
    ax.gridlines(draw_labels=False, linewidth=1, color='gray', alpha=0.5, linestyle='--')

    cbar = plt.colorbar(cb) #,  format=mticker.FuncFormatter(fmt))
    if style=='discrete':
        cbar.set_ticks(bounds_figure)
        mylabels = ["{:3.3e}".format(a) for a in bounds_figure]
        cbar.set_ticklabels(mylabels)
        

    else:
        cbar.set_ticks(p_clev)
        cbar.set_ticklabels(p_clev)
        mylabels=["{:3.3e}".format(a) for a in p_clev]
        cbar.set_ticklabels(mylabels)


    #new_strunits = units_adjust(str(data2D[varname].attrs['units'])) #.replace('kg','mg'))
    #cbar.ax.set_ylabel(new_strunits)
    #cbar.ax.set_ylabel(str(data2D[varname].attrs['units']))

    cbar.ax.set_ylabel(varname+'  ['+str(kwargs['units'])+']')
    staval= (data2Davg,np.amin(data2D[varname].values),np.amax(data2D[varname].values))
    strval= ' Mean=%6.3f  Min=%6.3f  Max=%6.3f' % staval
    plt.text(0.5, -0.05, strval,
             horizontalalignment='center',
             verticalalignment='top',
             transform=ax.transAxes,  fontsize=10)
    
    if kwargs['title']=='At-Level': 
       plt.title(kwargs['title']+'-'+extratitle+' [Pa] '+varname)
    elif kwargs['title']=='Surface': 
       plt.title(kwargs['title']+' at surface '+varname)
    elif kwargs['title']=='Column': 
       plt.title(kwargs['title']+' '+varname)


    plt.tight_layout()
    plt.savefig(pname)
    plt.close()
    return


def zonal_mean(ncname, varname, finfo=None, settings=None, config=None, plotinfo=None):

    dir_plots=config['info_plots']
    data = xr.open_dataset(ncname).mean(dim='time')
    dataplt = data.mean(dim='lon')

    name = os.path.split(ncname)[-1].replace('.nc', "_ZONAL_MEAN.png")
    pname= os.path.join(dir_plots,name)

    plot_zm = dataplt[varname].plot(yincrease=False)
    plt.title(plotinfo['title']+' '+varname+'  ['+str(plotinfo['units'])+']')
    plt.savefig(pname)
    plt.close()

    return varname
 
def surface_map(ncname, varname, finfo=None, settings=None, config=None, plotinfo=None):
    """

    """

    if plotinfo['title']=='Column':
       set_level=None
    else:
       set_level = plotinfo['set_level']

    if set_level==None:

      dir_plots=config['info_plots']
      data = xr.open_dataset(ncname).mean(dim='time')
      name = os.path.split(ncname)[-1].replace('.nc', "_MAP.png")
      pname= os.path.join(dir_plots,name)

      # Dic plot examples ==============================================
      # dic_plot = {'vmin':250, 'vmax':350}

      if plotinfo==None:
          dic_plot = {}
      else:
          dic_plot = plotinfo

      singlemap(data, varname, pname, swap=False, extratitle='', style='discrete', **dic_plot)

    else:

      dir_plots=config['info_plots']
      data = xr.open_dataset(ncname).mean(dim='time').isel(pres=set_level)
      print(data)
      name = os.path.split(ncname)[-1].replace('.nc', "_atlevel"+str(set_level)+".png")
      pname= os.path.join(dir_plots,name)
      strtitle = str(data['pres'].values)
      # Dic plot examples ==============================================
      # dic_plot = {'vmin':250, 'vmax':350}

      if plotinfo==None:
          dic_plot = {}
      else:
          dic_plot = plotinfo

      singlemap(data, varname, pname, swap=False, extratitle=strtitle, style='discrete', **dic_plot)


    return varname, varname, varname







