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
from  netCDF4  import Dataset
import os
import json
from pprint import pprint
import shutil
import yaml
import pandas as pd
import lib.check as check
import gc 


def _get_var_checkfiles_IPSL(outname, yr, study, freqsource, module, timeout, exp,
                             extra, storepath, dirMO, dirDA, dirHF):
        """
        Purpose: This function is specific of IPSL. For other models additional
                 functions has to be programmed.

                 It is indicated the typical form of the files found.
                 for each kind of 'module', therefore for each variable-id
                 the code found the netcdf files where the variable is stored.


        :param outname:
        :param yr:
        :param study:
        :param freqsource:
        :param module:
        :param timeout:
        :param exp:       
        :param extra:
        :param storepath:  where the raw netcdf files has to be found
        :param dirMO:
        :param dirDA:
        :param dirHF:

        :type outname: string
        :type yr: string
        :type study: string
        :type freqsource: string
        :type module: string
        :type timeout: string
        :type exp: string
        :type extra: string
        :type storepath:  string 
        :type dirMO: string
        :type dirDA: string
        :type dirHF: string

        :return:
        """

        #dic_struc = { 'ATM':{'MO':('_1M_histmth.nc','_mon'), 'DA':('_1M_histmth.nc','_mon'),
        #                     'HF':('_HF_histh1f.nc','_1hr')}

        f_tail   = '_IPSL_'+study+'_'+yr+'0101_'+yr+'1231'+'.nc'

        bpath = storepath+exp  # storepath is where the files are found

        if module=='ATM':
            basedir=bpath+'/ATM/Output/'+freqsource+'/'+exp+'_'+yr+'0101_'+yr+'1231'
            if freqsource=='MO':
                input_f=basedir +'_1M_histmth.nc'
                outfile = dirMO+outname+'/'+outname+'_mon'+f_tail
            if freqsource=='DA':
                input_f=basedir+'_1D_histday.nc'
                outfile = dirDA+outname+'/'+outname+'_day'+f_tail
            if freqsource=='HF':
                if timeout=='1h':
                    input_f = basedir+'_HF_histh1f.nc'
                    outfile    = dirHF+outname+'/'+outname+'_1hr'+f_tail
                if timeout=='6h':
                    input_f = basedir+'_HF_histh6f.nc'
                    outfile    = dirHF+outname+'/'+outname+'_6hr'+f_tail
                if timeout=='3h':
                    input_f = basedir+'_3H_remotesens.nc'
                    outfile    = dirHF+outname+'/'+outname+'_3hr'+f_tail
                if timeout=='xh':
                    input_f = basedir+'_HF_histhf.nc'
                    outfile    = dirHF+outname+'/'+outname+'_3hr'+f_tail


        if module=='CHM':
            if freqsource=='MO':
                if extra=='aerosols_from_inca.nc':
                    input_f=bpath+'/CHM/Output/MO/'+exp+'_'+yr+'0101_'+yr+'1231'+'_1M_'+extra
                else:
                    input_f=bpath+'/CHM/Output/MO/'+exp+'_'+yr+'0101_'+yr+'1231'+'_1M_inca_'+extra
                outfile = dirMO+outname+'/'+outname+'_mon'+f_tail

            if freqsource=='DA':
                input_f=bpath+'/CHM/Output/DA/'+exp+'_'+yr+'0101_'+yr+'1231'+'_1D_inca_'+extra
                outfile = dirDA+outname+'/'+outname+'_day'+f_tail

            if freqsource=='HF':
                input_f=bpath+'/CHM/Output/DA/'+exp+'_'+yr+'0101_'+yr+'1231'+'_1D_inca_'+extra
                if timeout=='1h':
                    outfile = dirHF+outname+'/'+outname+'_1h'+f_tail
                if timeout=='6h':
                    outfile = dirHF+outname+'/'+outname+'_6h'+f_tail

        if module=='SRF':
            input_f=bpath+'/SRF/Output/MO/'+exp+'_'+yr+'0101_'+yr+'1231'+'_1M_sechiba_history.nc'
            if freqsource=='MO':
                outfile = dirMO+outname+'/'+outname+'_mon'+f_tail


        #print(input_f)
        #print(outfile)
        return input_f, outfile

def solve_time_IPSL(outfile, outfilen, varname='', out=''):
    """

    :param outfile:
    :param outfilen:
    :param varname:
    :param out:
    :return:
    """


    commandnco1 = ('ncrename -O -d time_counter,time -v time_counter,time ' + outfile
                  + ' ' + outfilen + out)
    os.system(commandnco1)
    commandnco2 = 'ncrename -O -v .time_counter_bnds,time_bnds '+outfile+' '+outfilen+out
    os.system(commandnco2)
    commandnco3 = 'ncrename -O -d .presnivs,pres -v .presnivs,pres '+outfile+' '+outfilen+out
    os.system(commandnco3)
    commandnco4 = 'ncatted -O -a standard_name,pres,o,c,air_pressure '+outfile+' '+outfilen+out
    os.system(commandnco4)
    commandnco5 = 'ncatted -O -a bounds,time,o,c,time_bnds '+outfile+' '+outfilen+out
    os.system(commandnco5)
    commandnco6 = 'ncks -O -x -v time_centered '+outfile+' '+outfilen+out
    os.system(commandnco6)
    commandnco7 = 'ncks -O -x -v time_instant '+outfile+' '+outfilen+out
    os.system(commandnco7)
    '''
    print(commandnco1)
    print('\n')
    print(commandnco2)
    print('\n')
    print(commandnco3)

    print('\n')
    print(commandnco4)

    print('\n')
    print(commandnco5)

    print('\n')
    print(commandnco6)

    print('\n')
    print(commandnco7)
    '''
    if varname is not '':
        commandnco8 = 'ncatted -O -a coordinates,'+varname+',d,, '+outfile+' '+outfilen+out
        os.system(commandnco8)
        #print('\n')
        #print(commandnco8, varname)


    return


