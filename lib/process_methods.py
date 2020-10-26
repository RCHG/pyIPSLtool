"""
##############################################################################################

  Project   : CRESCENDO/AEROCOM
  Filename  : liborganize.py
  Author    : Ramiro Checa-Garcia
  email     : rcheca@lsce.ipsl.fr
  Purpose   : Reorganize variables for LMDzINCAOR experiment runs.

  Revision History ----------------------------------------------------------

  Date       Author     Ref    Revision

  2018-Apr   R.Checa           First code
  2018-Jul   R.Checa           Implemented as module
  2018-Nov   R.Checa           Final first version.


  TODO LIST:

  0  Function add_str_attributes --> info should not hard coded. 
  1. Diagnostics on lev_var not in hard-code         (Improvement)
  2. Add reference values for key diagnostics        (Improvement)
  3. Add description to each function                (Documentation)
  4. Factors in new_var could be variables in netcdf (Improvement)

##############################################################################################
"""

import xarray as xr
import numpy  as np
import glob
from   netCDF4  import Dataset
import os
import json
from   pprint import pprint
import shutil
import yaml
import pandas as pd
import lib.check as check
import gc 
from lib.specific_IPSL import solve_time_IPSL, _get_var_checkfiles_IPSL
from lib.process_support import _get_var_checkdir, add_global_attributes
from lib.process_support import check_units, create_filename
from lib.process_support import directory_structure, find_tableID
from lib.prints import myprint


def get_var(mysettings, year, varname, freqsource, module, timeout,
            outname, expcase, study, extra='',index=1,finfo=None, model='IPSL'):
    """

    :param mysettings:
    :param year:
    :param varname:
    :param freqsource:
    :param module:
    :param timeout:
    :param outname:
    :param expcase:
    :param study:
    :param extra:
    :param index:
    :param finfo:
    :return:
    """
    """
    - Module= ATM, CHM, etc
    - name=vmro3 etc
    - freq=DA, HF, MO for ATM, for CHM all is DA
    """

    dirMO, dirDA, dirHF = directory_structure(study, mysettings, finfo=finfo)

    cwd  = mysettings['outdir']['path']
    subf = mysettings['outdir']['subfolder']
    storepath = mysettings['storepath']

    status = _get_var_checkdir(outname, freqsource, dirMO, dirDA, dirHF)

    if model=='IPSL':
       input_f, outfile  = _get_var_checkfiles_IPSL(outname, year, study, freqsource, module,
                                                    timeout, expcase, extra, storepath,
                                                    dirMO, dirDA, dirHF)
    else:
       myprint('No other models has been implemented'+'\n'+'... check your config file', finfo=finfo)

    
    str1=('    Extracting ('+str(index).rjust(4)+'): '+varname.ljust(20)+
          ' and save as '+outname.ljust(20)+' .... for year '+
          year + '  ... ['+freqsource+' '+module+']')

    myprint(str1, finfo=finfo)


    commandcdo =('cdo  --silent selvar,' + varname + ' ' + input_f + ' ' + outfile
                + ' &> delfiles/del' + varname + freqsource + module)
    os.system(commandcdo)

    if model=='IPSL':
      solve_time_IPSL(outfile, outfile, varname=varname, out=' &>> delfiles/del'+
                      varname+freqsource+module)

    if outname!=varname:
      commandcdo ='cdo --silent -O chname,'+varname+','+outname+' '+outfile+' '+outfile+'x'
      os.system(commandcdo)
      os.rename(outfile+'x', outfile)

    #commandncinfo = 'ncinfo -v '+varname+' '+outfile+' | grep shape'
    #print(outname, varname)

    return

def scale_variable(vardataset, d_scale_vars, outname, varname, outfile, finfo=None, model='IPSL'):
    """

    :param vardataset:
    :param d_scale_vars:
    :param outname:
    :param varname:
    :param outfile:
    :param finfo:
    :return:
    """

    var_initial =  d_scale_vars[outname]['var_initial']
    if 'scale_factor' in d_scale_vars[outname].keys():
        varscaled = vardataset[var_initial]*d_scale_vars[outname]['scale_factor']
    else:
        varscaled = vardataset[var_initial]

    # Check if global unit factor and new units is indicated
    if 'new_units' in d_scale_vars[outname].keys():
        units_factor = d_scale_vars[outname]['new_units']['units_factor']
        units_name   = d_scale_vars[outname]['new_units']['units_name']
        myprint('      We will change units from '+vardataset[varname].attrs['units']+
                '  to  '+units_name, finfo=finfo)
        myprint('      We use a unit conversion factor of '+str(units_factor), finfo=finfo) 
    else:
        units_factor = 1.0
        try:
            units_name   = vardataset[varname_sur[0]].attrs['units']
        except IndexError:
            units_name   = vardataset[varname_alt[0]].attrs['units']

    varfinal = varscaled*units_factor
    varfinal.name = outname
    varfinal.attrs['units']=units_name
    varfinal.to_netcdf(outfile, unlimited_dims='time_counter')
    varfinal.close()

    return
 

def new_variable_altsur(vardataset, d_new_vars, outname, varname, outfile, finfo=None, model='IPSL', extraname=''):
    """

    :param vardataset:
    :param d_new_vars:
    :param outname:
    :param varname:
    :param outfile:
    :param finfo:
    :return:
    """


    varname_sur = d_new_vars[outname]['var_surface']
    varname_alt = d_new_vars[outname]['var_altitud']
    factors_sur = [str(y) for y in d_new_vars[outname]['fac_surface']]
    factors_alt = [str(y) for y in d_new_vars[outname]['fac_altitud']]

    lunits =check_units(varname_sur+varname_alt, vardataset)
    if len(set(lunits))==1:
        myprint('       internal consistency in units is ... ok ... '+lunits[0], finfo=finfo)

    else:
        myprint(str(lunits)+str(varname_sur)+str(varname_alt), finfo=finfo)
        exit()

    # Aggregate values over vertical coordinates
    for ivar, var_alt in enumerate(varname_alt):
        if var_alt in vardataset.keys():
            if ivar==0:
                var_sum_alt = vardataset[var_alt].sum(dim='presnivs')
            else:
                var_sum_alt = (var_sum_alt +
                 vardataset[var_alt].sum(dim='presnivs')*d_new_vars[outname]['fac_altitud'][ivar])
        else:
            print(extraname)
            altdataset = xr.open_dataset(extraname[0])   # we assume that all variable are in same file.
            if ivar==0:
                var_sum_alt = altdataset[var_alt].sum(dim='presnivs')
            else:
                var_sum_alt = (var_sum_alt +
                 altdataset[var_alt].sum(dim='presnivs')*d_new_vars[outname]['fac_altitud'][ivar])

    # Aggregate values on surface
    for ivar, var_sur in enumerate(varname_sur):
        if 'presnivs' in vardataset[var_sur].dims:
            myprint('       selecting the pressure ... '+str(vardataset['presnivs'].values[0])+
                    ' '+str(vardataset['presnivs'].attrs['units']), finfo=finfo)
            if ivar==0:
                var_sum_sur  = vardataset[var_sur].isel(presnivs=0)
            else:
                var_sum_sur  = (var_sum_sur +
                                vardataset[var_sur].isel(presnivs=0)*d_new_vars[outname]['fac_surface'][ivar])
        else:
            if ivar==0:
                var_sum_sur  = vardataset[var_sur]
            else:
                var_sum_sur  = (var_sum_sur +
                                vardataset[var_sur]*d_new_vars[outname]['fac_surface'][ivar])

    # Check if global unit factor and new units is indicated
    if 'new_units' in d_new_vars[outname].keys():
        units_factor = d_new_vars[outname]['new_units']['units_factor']
        units_name   = d_new_vars[outname]['new_units']['units_name']
        myprint('      We will change units from '+
                vardataset[varname_sur[0]].attrs['units']+
                '  to  ' + units_name, finfo=finfo)
        myprint('      We use a unit conversion factor of '+str(units_factor), finfo=finfo) 
    else:
        units_factor = 1.0
        try:
            units_name   = vardataset[varname_sur[0]].attrs['units']
        except IndexError:
            units_name   = vardataset[varname_alt[0]].attrs['units']

    # Aggregate all together:
    if len(varname_sur)>=1 and len(varname_alt)>=1:
        varfinal = (var_sum_alt+var_sum_sur)*units_factor
    if len(varname_alt)==0:
        varfinal = var_sum_sur*units_factor
    if len(varname_sur)==0:
        varfinal = var_sum_alt*units_factor

    varfinal.name = outname
    varfinal.attrs['units']=units_name
    varfinal.to_netcdf(outfile, unlimited_dims='time_counter')
    varfinal.close()

    return


def new_var(mysettings, year, varname, freqsource, module, timeout, outname,
            expcase, study, extra='',index=1, finfo=None, model='IPSL'):
    """

    :param mysettings:
    :param year:
    :param varname:
    :param freqsource:
    :param module:
    :param timeout:
    :param outname:
    :param expcase:
    :param study:
    :param extra:
    :param index:
    :param finfo:
    :return:
    """

    """
    - Module= ATM, CHM, etc
    - name=vmro3 etc
    - freq=DA, HF, MO for ATM, for CHM all is DA
    """

    dirMO, dirDA, dirHF = directory_structure(study, mysettings)

    cwd  = mysettings['outdir']['path']
    subf = mysettings['outdir']['subfolder']
    storepath = mysettings['storepath']

    status = _get_var_checkdir(outname, freqsource, dirMO, dirDA, dirHF)

    if model=='IPSL':
       input_f, outfile  = _get_var_checkfiles_IPSL(outname, year, study, freqsource, module,
                                                    timeout, expcase, extra, storepath,
                                                    dirMO, dirDA, dirHF)
    else:
       myprint('No other models has been implemented'+'\n'+'... check your config file', finfo=finfo)

    vardataset = xr.open_dataset(input_f, chunks={'time_counter': 10})  # we assume that all variable are in same file.
    altname    = input_f.replace('emi.nc','species.nc'),
    f_new_vars = open("etc/new_variables.yaml", 'r')
    d_new_vars = yaml.load(f_new_vars)

    f_scale_vars = open("etc/scale_variables.yaml", 'r')
    d_scale_vars = yaml.load(f_scale_vars)

    if outname in d_scale_vars.keys():
        str1=(  '    Scaling   ('+str(index).rjust(4)+'): '+varname.ljust(20)+
                ' and save as '+outname.ljust(20)+
                ' .... for year '+ year + '  ... ['+freqsource+' '+module+']\n')
        myprint(str1, finfo=finfo)

        scale_variable(vardataset, d_scale_vars, outname, varname, outfile, finfo)

    if outname in d_new_vars.keys():
        varname_sur = d_new_vars[outname]['var_surface']
        varname_alt = d_new_vars[outname]['var_altitud']
        factors_sur = [str(y) for y in d_new_vars[outname]['fac_surface']]
        factors_alt = [str(y) for y in d_new_vars[outname]['fac_altitud']]

        str1=(  '    Creating   ('+str(index).rjust(4)+'): '+varname.ljust(20)+' and save as '
              +outname.ljust(20)+' .... for year '
              + year + '  ... ['+freqsource+' '+module+']\n')
        str2=(   '       using files ... '+', '.join(varname_sur)+
              '\n        ... with factors '+', '.join(factors_sur)+
              '\n       and files   ... '+', '.join(varname_alt)+
              '\n        ... with factors '+', '.join(factors_alt))

        myprint(str1+str2, finfo=finfo)

        new_variable_altsur(vardataset, d_new_vars, outname, varname, outfile, finfo, extraname=altname)

    if model=='IPSL':
      solve_time_IPSL(outfile, outfile, varname=varname, 
                      out=' &>> delfiles/del_'+varname+freqsource+module)

    return

def lev_var(mysettings, year, varname, freqsource, module, timeout, outname,
            expcase, study, extra='', lev='top', index=1, finfo=None, model='IPSL'):
    """

    :param mysettings:
    :param year:
    :param varname:
    :param freqsource:
    :param module:
    :param timeout:
    :param outname:
    :param expcase:
    :param study:
    :param extra:
    :param lev:
    :param index:
    :param finfo:
    :return:
    """

    """
    - Module= ATM, CHM, etc
    - name=vmro3 etc
    - freq=DA, HF, MO for ATM, for CHM all is DA
    """

    dirMO, dirDA, dirHF = directory_structure(study, mysettings, finfo=finfo)

    cwd  = mysettings['outdir']['path']
    subf = mysettings['outdir']['subfolder']
    storepath = mysettings['storepath']

    status = _get_var_checkdir(outname, freqsource, dirMO, dirDA, dirHF)

    if model=='IPSL':
       input_f, outfile  = _get_var_checkfiles_IPSL(outname, year, study, freqsource, module,
                                                    timeout, expcase, extra, storepath,
                                                    dirMO, dirDA, dirHF)
    else:
       myprint('No other models has been implemented'+'\n'+'... check your config file', finfo=finfo)

    vardataset = xr.open_dataset(input_f)  # we assume that all variable are in same file.

    if outname in ['rlut']:
        varname_lev = ['rlu']

    if outname in ['rlutcs']:
        varname_lev = ['rlucs']


    str1=('    Creating   ('+str(index).rjust(4)+'): '+varname.ljust(20)+
          ' and save as '+outname.ljust(20)+' .... for year '+
          year + '  ... ['+freqsource+' '+module+']')
    str2=('       using files ... '+','.join(varname_lev))
    myprint(str1+'\n'+str2, finfo=finfo)


    if lev=='top':
    #from  netCDF4  import Dataset
        presval = -1
    elif lev=='surf':
        presval = 0
    else:
        print('Problem in lev_var')
        exit()

    var_lev  = vardataset[varname_lev[0]].isel(presnivs=presval)
    myprint('       selecting the pressure ... '+str(vardataset['presnivs'].values[presval]), finfo=finfo)

    var_lev.name = outname
    var_lev.to_netcdf(outfile, unlimited_dims='time_counter')
    if model=='IPSL':
      solve_time_IPSL(outfile, outfile, varname=varname, 
                      out=' &>> delfiles/del'+varname+freqsource+module)

    return

