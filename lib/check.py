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

def myprint(mystr, finfo=None):
    print(mystr)
    if finfo!=None:
        finfo.write(mystr+'\n')
    return

def total_load(ncname, varname, dicinfo,  area_fname='data/area_grid.nc', area_varid='area', finfo=None):
    """


    """

    data = xr.open_dataset(ncname)
    areadata = xr.open_dataset(area_fname)[area_varid]
    '''
    try:
        print(area_fname, area_varid)
        areadata = xr.open_dataset(area_fname)[area_varid]
    except:
        import iris
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == varname))
        cube = iris.load(fn_process, constraints=variable_constraint)[0]
        try:
           cube.coord('latitude').guess_bounds()
           cube.coord('longitude').guess_bounds()
        except ValueError:
           pass
    cube_area = iris.analysis.cartography.area_weights(cube)
    areadata=cube_area[0].data
    '''
    load_field = data[varname]*areadata
    load_month = load_field.sum(['lat','lon'])
    load_out = load_month.mean(dim='time')/1.e9

    myprint('       ... [checking values] [Tg]: '+ str(load_out.values), finfo=finfo)

    return varname, 'Tg', load_out.values



def _total_tendency(tendency_data, varname, area):
    """


    """

    tendency_field = tendency_data[varname]*area
    tendency_month = tendency_field.sum(['lat','lon'])
    seconds_month = np.array([k * 86400. for k in [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]])
    tendency_total= np.array([tendency_month[imonth]*seconds_month[imonth] for imonth in range(12)])

    acc_tot = np.round(tendency_total.sum()/1.e9,3)
    acc_mon = np.round(tendency_total/1.e9,1)

    return acc_tot, acc_mon

def tendency_mass(ncname, varname, area_fname, area_varid, finfo=None):

    data = xr.open_dataset(ncname)
    area = xr.open_dataset(area_fname)[area_varid]

    new = data.groupby('time.month')
    new1= new.mean(dim='time')
    new1.coords['month']=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    val_tendency, monthly = _total_tendency(new1, varname, area)

    myprint('       ... [checking values] [Tg yr-1]: '+ str(val_tendency), finfo=finfo)
    myprint('       ... [checking values]   monthly: '+ ', '.join([str(y) for y in monthly]), finfo=finfo)

    return varname, 'Tg yr-1', val_tendency

def global_mean(ncname, varname, dicinfo, finfo=None, level=None):

    dg_to_rd = np.pi/180.0001
    data = xr.open_dataset(ncname) # time, lat, lon
    if level==None:
        mymethod='global_mean'
        datavar = data[varname].values
        vlevel=''
    else:
        mymethod='global_mean_atlev'
        #print(data[varname])
        datavar = data[varname].isel(pres=level).values
        myprint('       evaluating at the pressure ... '+str(data['pres'].values[level])+
                    ' '+str(data['pres'].attrs['units']), finfo=finfo)
        vlevel = ' '+str(data['pres'].values[level])+' '+str(data['pres'].attrs['units'])


    lat = data.coords['lat'].values
    lon = data.coords['lon'].values
    time= data.coords['time'].values
    lat_nh = np.where(lat>0)
    lat_sh = np.where(lat<0)
    lataxis=0
    lonaxis=1
    ts_global=[]
    ts_nhemis=[]
    ts_shemis=[]
    for itime in range(0,len(time)):
        #print(datavar.shape)
        data2Dgb = np.squeeze(datavar[itime,:     ,:])
        data2Dnh = np.squeeze(datavar[itime,lat_nh,:])
        data2Dsh = np.squeeze(datavar[itime,lat_sh,:])
        data2Davggb = np.average(np.average(data2Dgb, axis=lonaxis),
                                 axis=lataxis, weights=np.cos(lat*dg_to_rd))

        data2Davgnh = np.average(np.average(data2Dnh, axis=lonaxis),
                                 axis=lataxis, weights=np.cos(lat[lat_nh]*dg_to_rd))

        data2Davgsh = np.average(np.average(data2Dsh, axis=lonaxis),
                                 axis=lataxis, weights=np.cos(lat[lat_sh]*dg_to_rd))

        ts_global.append(data2Davggb)
        ts_nhemis.append(data2Davgnh)
        ts_shemis.append(data2Davgsh)

    a_gl = np.array(ts_global) #; a_gl_avg = str(np.round(a_gl.mean(),4))
    a_nh = np.array(ts_nhemis) #; a_nh_avg = str(np.round(a_nh.mean(),4))
    a_sh = np.array(ts_shemis) # ; a_sh_avg = str(np.round(a_sh.mean(),4))
    #a_gl_min = "{0:0.3f}".format(a_gl.min())
    #a_gl_max = "{0:0.3f}".format(a_gl.max())
    #a_gl_std = "{0:0.3f}".format(a_gl.std())

    a_gl1 = "> 90S-90N {0:>8.5g} | NH  {1:>8.5g} | SH  {2:>8.5g} ".format(a_gl.mean(), a_nh.mean(), a_sh.mean())
    a_gl2 = "| min     {0:>8.5g} | max {1:>8.5g} | std {2:>8.5g} ".format(a_gl.min(), a_gl.max(), a_gl.std())
    myprint('       ... [checking values]: '+
            '['+str(data[varname].attrs['units']).rjust(20)+']'+a_gl1+a_gl2, finfo=finfo)

    if float(dicinfo['method'][mymethod]['factor'])!=1:
       na_gl =  a_gl*float(dicinfo['method'][mymethod]['factor'])
       na_nh =  a_nh*float(dicinfo['method'][mymethod]['factor'])
       na_sh =  a_sh*float(dicinfo['method'][mymethod]['factor'])

       na_gl1 = "> 90S-90N {0:>8.5g} | NH  {1:>8.5g} | SH  {2:>8.5g} ".format(na_gl.mean(), na_nh.mean(), na_sh.mean())
       na_gl2 = "| min     {0:>8.5g} | max {1:>8.5g} | std {2:>8.5g} ".format(na_gl.min(),  na_gl.max(),  na_gl.std() )

       myprint('       ... [checking values]: '+
               '['+str(dicinfo['method'][mymethod]['units']).rjust(20)+']'
               +na_gl1+na_gl2, finfo=finfo)

    else:
        na_gl= a_gl
    return varname, dicinfo['method'][mymethod]['units'], na_gl.mean(), vlevel


def global_mean_atlev(ncname, varname, dicinfo, level=0, finfo=None):

    return global_mean(ncname, varname, dicinfo, finfo=finfo, level=level)






