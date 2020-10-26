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
new_var  emibc           MO           CHM          -              emibc        emi.nc
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
from lib.process_methods import lev_var, get_var, new_var, new_variable_altsur
from lib.specific_IPSL import solve_time_IPSL, _get_var_checkfiles_IPSL

def memory_usage_psutil():
    # return the memory usage in MB
    import psutil
    return psutil.virtual_memory()

def myprint(mystr, finfo=None):
    """myprint: just a python to manage the log prints

    :param mystr:
    :param finfo:
    """

    if finfo!=None:
        finfo.write(mystr+'\n')
    else:
        print(mystr)
    return

def find_tableID(varID, freq='mon', origin='ATM', dirtables=''):
    """find_tableID: This function look the varID in cmor-cmip6 tables 
                     recursively. When it finds the table it is using
                     it to get/check metatada.

    :param varID:
    :param freq:
    :param origin:
    :param dirtables:
    """

    table_out=False
    ltable=[]
    varIDwith = "'"+varID+"'"
    cmip6_tables = 'etc/cmip6-cmor-tables-master/Tables/'
    for fname in os.listdir(cmip6_tables):               # change directory as needed
        #print(fname)
        if os.path.isfile(cmip6_tables+fname):           # it's a file, not a directory entry
            with open(cmip6_tables+fname) as f:          # open file
                for line in f:                           # process line by line
                    if varIDwith in line:                # search for string
                        table_name = cmip6_tables+fname
                        ltable.append(table_name)
                        break

    varIDwith = '"'+varID+'"'
    for fname in os.listdir(cmip6_tables):               # change directory as needed
        #print(fname)
        if os.path.isfile(cmip6_tables+fname):           # it's a file, not a directory entry
            with open(cmip6_tables+fname) as f:          # open file
                for line in f:                           # process line by line
                    if varIDwith in line:                # search for string
                        table_name = cmip6_tables+fname
                        ltable.append(table_name)
                        break

    if len(ltable)==1:
        table_out = ltable[0]
    else:
        l2table=[]
        for table_name in ltable:
            if freq in table_name:
                l2table.append(table_name)
        if len(l2table)==1:
            table_out=l2table[0]
        else:
            for table_name in l2table:
                if 'CMIP6_AER' in table_name:
                    table_out = table_name
            if table_out==False:
                for table_name in l2table:
                    if 'CMIP6_A' in table_name:
                        table_out = table_name

    return table_out


def create_filename(mysettings, study, expID, varID, varfile, freq='mon', finfo=''):
    """ 
    Purpose: This function creates the new filename of the netcdf with the diagnostics 
             as derived from the current settings. In the case of CMIP6 netcdfs the 
             filename will have the structure:

             **varID_tableID_sourceID_experimentID_memberID_gridlevel_timerange**

             * varID  -------> is the variable name: varname
             * tableID ------> depends on variable AERday, Amon etc...
             * sourceID -----> IPSL-LMDZORINCAv6
             * experimentID -> CRESCWP3PD-amip
             * memberID -----> v5-r1i1p1f1  (v5 refers to our internal version)
             * gridlabel ----> gr           (refers to the original grid but we provide at preslev)
             * timerange ----> 20000101-20121231 for all cases.


    :param mysettings:
    :param study:
    :param expID:
    :param varID:   is the variable name: varname
    :param varfile:
    :param freq:
    :param finfo:

    :return:
        - save_name
        - info_var
    """

    table_name   = find_tableID(varID, freq=freq)
    if table_name == False:
        print('         ---> variable: ',varID, ' not found in CMIP6-cmor-tables\n')
        return False, False
    tableID      = table_name.split('/')[-1].split('.')[0].split('CMIP6')[1][1::]
    sourceID     = mysettings[study]['expCVs']['sourceID']
    experimentID = mysettings[study]['expID']
    memberID     = mysettings[study]['expCVs']['memberID']
    gridlabel    = mysettings[study]['expCVs']['gridlabel']
    timerange    = varfile.split('_')[-1]

    if mysettings[study]['expKIND'].upper() in ['CMIP6','CRESCENDO','JASPERKOK']:
        # for CMIP6 or CRESCENDO
        save_name = '_'.join([varID ,tableID, sourceID, experimentID, memberID,
                              gridlabel, timerange
                              ])

    if 'aerocom' in mysettings[study]['expKIND'].lower():
        # for aerocom-phaseIII
        # aerocom3_<ModelName>_<ExperimentName>_<VariableName>_<VerticalCoordinateType>
        #    _<Period>_<Frequency>.nc
        # expID is AP3-CTRL2016-PD or AP3-CTRL2016-PI or AP3-remotesensing etc..
        save_name = ('aerocom3_' + sourceID + '_' + mysettings[study]['expID'] + '_'
                     + varID + '_' + timerange.replace('.nc','') + '_' + freq + '.nc')

    with open(table_name, 'r') as cmip6tb:
        tbdata   = json.load(cmip6tb)
        info_var = tbdata['variable_entry'][varID]
    myprint('       ... final filename will be :'+ save_name, finfo=finfo)
    return save_name, info_var


def add_global_attributes(mysettings, ncname, newname, study, info_var, finfo=None):
    """
    This function should be improved

    :param mysettings:
    :param ncname:
    :param newname:
    :param study:
    :param info_var:
    :param finfo:
    :return:
    """

    dataset = xr.open_dataset(ncname)
    #, autoclose=True) autoclose not implemented in IRENE old xarray version
    # Here we opened a dataset 
    dataset.attrs['activity_id']   = mysettings[study]['expCVs']['activityID']
    dataset.attrs['institution_id']='IPSL-LSCE'
    dataset.attrs['source_id']     ='IPSL-LMDZORINCAv6'
    dataset.attrs['source_type']   ='AGCM-AER'
    dataset.attrs['realm']         ='aerosols'
    dataset.attrs['experiment_id'] = mysettings[study]['expID']
    dataset.attrs['packagedby']    ='pyIPSLpack (contact: rcheca@lsce.ipsl.fr)'
    dataset.attrs['contact']       = mysettings[study]['expCVs']['contactID']
    dataset.attrs['spinning-up']   = mysettings[study]['expCVs']['spin-up']

    for field in ['long_name', 'standard_name', 'modeling_realm', 'frequency',
                  'dimensions', 'comment']:
        dataset[info_var['out_name']].attrs[field]=info_var[field]

    if dataset[info_var['out_name']].attrs['units']!=info_var['units']:
        # Here it would be possible to create a file with units equivalences...
        lipsl = ['-', '-',         'kg/kg',   'kg/m2' , 'kg/m2/s',    'kg/(s*m2)' , 'kg/m3' , '-', 'm^-1', 'W/m2', 'kg/kg',   'm/s']
        l_cf  = ['%', 'mol mol-1', '1'    ,   'kg m-2', 'kg m-2 s-1', 'kg m-2 s-1', 'kg m-3', '1', 'm-1' , 'W m-2','kg kg-1', 'm s-1']

        for valipsl, valcf in zip(lipsl, l_cf):
            if dataset[info_var['out_name']].attrs['units']==valipsl and info_var['units']==valcf:
               dataset[info_var['out_name']].attrs['units']=info_var['units']
    if dataset[info_var['out_name']].attrs['units']!=info_var['units']:
        str1=('        non CMIP6 units '+ dataset[info_var['out_name']].attrs['units']+','+info_var['units'])
        str2=('        keeping non CMIP6')
        myprint(str1+'\n'+str2+'\n', finfo=finfo)

        unit_status=False
    else:
        str1=('      ... Units consistent with CF')
        myprint(str1, finfo=finfo)

        unit_status=True
    ltest = ['times_bnds', 'time_bnds', 'time_instant_bounds']
    for checkbnd in ltest:
        if 'bounds' in dataset['time'].attrs.keys():
            if dataset['time'].attrs['bounds']==checkbnd:
                if  checkbnd not in dataset.data_vars.keys():
                    del dataset['time'].attrs['bounds']
                    #print('removed?', dataset['time'].attrs)

    if 'AEROCOM' in newname.upper():
        if 'lat' in dataset.coords and 'lon' in dataset.coords and 'pres' in dataset.coords:
            finalname = newname.replace(info_var['out_name']+'_', info_var['out_name']+'_3D_')
        elif 'lat' in dataset.coords and 'lon' in dataset.coords and 'pres' in dataset.coords:
            finalname = newname.replace(info_var['out_name']+'_', info_var['out_name']+'_3D_')
        elif 'lat' in dataset.coords and 'lon' in dataset.coords:
            finalname = newname.replace(info_var['out_name']+'_', info_var['out_name']+'_2D_')
        else:
            finalname=newname
    else:
        finalname=newname
    #print('going to save', info_var['out_name'] )
    #if '3d' not in info_var['out_name']:
    if os.path.isfile(finalname):
        print('File ',finalname, ' is there')
        print(os.access(finalname, os.W_OK))

    if 'airmass' in finalname or 'ec550' in finalname:
        print('')
    else:
        dataset.to_netcdf(finalname, mode='w')
        dataset.close()
    #else:
    #    shutil.copy2(ncname,finalname)
    return unit_status, finalname

def directory_structure(study, mysettings, clean=False, create=False, finfo=None, model='IPSL'):
    """

    :param study:
    :param mysettings:
    :param clean:
    :param create:
    :param finfo:
    :return:
    """

    import shutil
    from os.path import join as pjoin
    cwd  = mysettings['outdir']['path']
    subf = mysettings['outdir']['subfolder']

    filename1 = pjoin(cwd, subf, study, 'monthly/')
    filename2 = pjoin(cwd, subf, study, 'daily/')
    filename3 = pjoin(cwd, subf, study, 'hourly/')

    if clean==True:
         for filename in [filename1, filename2, filename3]:
            myprint('    Cleaning of: '+ filename, finfo=finfo)
            try:
                listf = os.listdir(filename)
                if mysettings['safety']['overwrite']==False and len(listf)>0:
                    myprint(' Safety settings stop the program: out-tree is not empty: \n\n', finfo=finfo)
                    myprint(str(listf), finfo=finfo)
                    myprint('-------\n'*2, finfo=finfo)
                else:
                    shutil.rmtree(filename)
            except OSError:
                myprint('    ... ask to clean an empty directory ...', finfo=finfo)
                pass

    if create==True:
        for filename in [filename1, filename2, filename3]:
            try:
                os.makedirs(filename)
            except OSError:
                pass


    return filename1, filename2, filename3


def _get_var_checkdir(outname, freqsource, dirMO, dirDA, dirHF):
    """

    :param outname:
    :param freqsource:
    :param dirMO:
    :param dirDA:
    :param dirHF:
    :return:
    """

    if freqsource=='DA':
        try:
            os.stat(dirDA+outname)
        except:
            os.mkdir(dirDA+outname)
    if freqsource=='MO':
         try:
            os.stat(dirMO+outname)
         except:
            os.mkdir(dirMO+outname)
    if freqsource=='HF':
         try:
            os.stat(dirHF+outname)
         except:
            os.mkdir(dirHF+outname)

    return True

def check_units(varname_list, vars_dataset):
    """

    :param varname_list:
    :param vars_dataset:
    :return:
    """

    lunits = []
    for varn in varname_list:
        lunits.append(vars_dataset[varn].attrs['units'])

    return lunits


def read_process_file(mysettings, process_file, year, expcase, study, finfo=None, model='IPSL'):
    """
      Read the process-file and for each line in the process-file 
      perform the process of the category indicated: new_var, get_var, lev_var

    :param myconfig:
    :param mysettings:
    :param process_file:
    :param year:
    :param expcase:
    :param study:
    :param finfo:
    :return:
    """

    import pandas as pd
    pd.set_option('display.width', 120)
    df = pd.read_table(process_file, skiprows=0, header=0, delim_whitespace=True)
    two_times = df.duplicated(subset='varname-out')
    if True in two_times.values:
        myprint('\n\n\n ============= CHECK your diagnostics file '+process_file
                +'\n ============= it has two repeated diagnostics final varnames!!\n\n\n',
                finfo=finfo)
        exit()
    list_vars = []
    for index, row in df.iterrows():
        gc.collect()

        list_vars.append(row['varname-out'])
        if mysettings[study]['preprocess']==True:
            if row['method']=='new_var':
                new_var(mysettings, year, row['varname-source'], row['freq-source'],
                        row['module-IPSL'], row['freq-specific'], row['varname-out'],
                        expcase, study, extra=row['extra-info-source'], 
                        index=index, finfo=finfo, model=model)
            if row['method']=='get_var':
                get_var(mysettings, year, row['varname-source'], row['freq-source'],
                        row['module-IPSL'], row['freq-specific'], row['varname-out'],
                        expcase, study, extra=row['extra-info-source'], 
                        index=index, finfo=finfo, model=model)
            if row['method']=='lev_var':
                lev_var(mysettings, year, row['varname-source'], row['freq-source'],
                        row['module-IPSL'], row['freq-specific'], row['varname-out'],
                        expcase, study, extra=row['extra-info-source'], 
                        index=index, finfo=finfo, model=model)
    return list_vars

def process_files(myconfig, mysettings, year, expcase, study, proc_files, finfo=None, model='IPSL'):
      """ This is the initial function of the processing step. It just call the 
          read_process_file for each process-file (remember files ended in .process)
          the process file has 3 categories: get_var, new_var and lev_var


      :param myconfig:
      :param mysettings:
      :param year:
      :param expcase:
      :param study:
      :param proc_file:
      :param finfo:
      :return:
         total_processed_vars
      """
      total_processed_vars = []
      if isinstance(proc_files, str):
         lproc_file=[proc_files]
      else:
         lproc_file=proc_files

      for proc_file in lproc_file:
        if str(proc_file)!='none':
            myprint('    ... aprox. number variables: '+ str(sum(1 for line in open(proc_file))-1), finfo=finfo)
            myprint('    ... aprox. number variables: '+ str(sum(1 for line in open(proc_file))-1), finfo=None)
            processed_vars = read_process_file(mysettings, proc_file, year, expcase, study, finfo=finfo,
                                               model=model)
        else:
            myprint('    ... No processing ... [process-file==None]', finfo=finfo)
            myprint('    ... No processing ... [process-file==None]', finfo=None)
            processed_vars=[]
        total_processed_vars = total_processed_vars + processed_vars

      return total_processed_vars


