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
import lib.plots as plots

def memory_usage_psutil():
    # return the memory usage in MB
    import psutil
    return psutil.virtual_memory()


def myprint(mystr, finfo=None):
    print(mystr)
    if finfo!=None:
        finfo.write(mystr+'\n')
    return

def check_var_netcdf(vars_checks, l_non_processed, l_yes_processed, l_uni_processed, l_yes_checked, check_list, 
                     directory, direc, mysettings, study, expID, varid, fname_input, post_t, 
                     freq='mon', finfo='', myconfig=None):

    if post_t=='year' or 'AEROCOM' in expID:
        new_name = fname_input.replace('1_20','1_20')
    else:
        new_name = fname_input

    good_name, info_var = create_filename(mysettings, study, expID, varid, new_name, freq=freq, finfo=finfo)
    
    if good_name!=False:
        save_name = directory+direc+'/'+good_name
        unit_stat, myname = add_global_attributes(mysettings, new_name, save_name, study, info_var)
        l_yes_processed.append(direc)
        l_uni_processed.append(unit_stat)
        for checkl in check_list:
          if varid in checkl.keys():
            l_yes_checked.append(varid)
            for mymethod in  checkl[varid]['method'].keys():
                if mymethod=='tendency_mass':
                    a1, a2, a3 = check.tendency_mass(myname, varid, 'etc/area_grid.nc', 'area', finfo=finfo)
                    vars_checks.append({'variable': a1, 'method': mymethod,'value': a3, 'units':a2,
                                        'testing-file':myname.split('/')[-1]})

                elif mymethod=='global_mean':
                    a1, a2, a3, a4= check.global_mean(myname, varid, checkl[varid], finfo=finfo)
                    vars_checks.append({'variable': a1, 'method': mymethod, 'value': a3, 'units':a2,
                                        'testing-file':myname.split('/')[-1]})

                elif mymethod=='global_mean_atlev':
                    for level in checkl[varid]['method'][mymethod]['levels']:
                        a1, a2, a3 , a4= check.global_mean_atlev(myname, varid,  checkl[varid], level=level, finfo=finfo)
                        vars_checks.append({'variable': a1, 'method': 'glo mean at'+a4,
                                            'value': a3, 'units':a2, 'testing-file':myname.split('/')[-1]})

                elif mymethod=='total_load':
                    a1, a2, a3  = check.total_load(finalname, varid, checkl[varid], 
                                                   area_fname='data/area_grid.nc',  area_varid='area', finfo=finfo)
                    vars_checks.append({'variable': a1, 'method': mymethod,  'value': a3, 'units':a2,
                                        'testing-file':myname.split('/')[-1]})
                else:
                    l_non_processed.append(direc)

    else:
      l_non_processed.append(varid)
      for checkl in check_list:
         if varid in checkl.keys():
            l_yes_checked.append(varid)
            myname = new_name
            unit_stat = True    # take care of this. Dangerous.

            for mymethod in  checkl[varid]['method'].keys():
                if mymethod=='tendency_mass':
                    a1, a2, a3 = check.tendency_mass(myname, varid, 'etc/area_grid.nc', 'area', finfo=finfo)
                    vars_checks.append({'variable': a1, 'method': mymethod,'value': a3, 'units':a2,
                                        'testing-file':myname.split('/')[-1]})

                elif mymethod=='global_mean':
                    a1, a2, a3, a4= check.global_mean(myname, varid, checkl[varid], finfo=finfo)
                    vars_checks.append({'variable': a1, 'method': mymethod, 'value': a3, 'units':a2,
                                        'testing-file':myname.split('/')[-1]})

                elif mymethod=='global_mean_atlev':
                    for level in checkl[varid]['method'][mymethod]['levels']:
                        a1, a2, a3 , a4= check.global_mean_atlev(myname, varid,  checkl[varid], level=level, finfo=finfo)
                        vars_checks.append({'variable': a1, 'method': 'glo mean at'+a4,
                                            'value': a3, 'units':a2, 'testing-file':myname.split('/')[-1]})

                elif mymethod=='total_load':
                    a1, a2, a3  = check.total_load(finalname, varid, checkl[varid], 
                                                   area_fname='data/area_grid.nc',  area_varid='area', finfo=finfo)
                    vars_checks.append({'variable': a1, 'method': mymethod,  'value': a3, 'units':a2,
                                        'testing-file':myname.split('/')[-1]})
                else:
                    l_non_processed.append(direc)

    if mysettings[study]['plotchecks']==True:
       list_files_plots = glob.glob('etc/plots/plots_*.yaml')
       plotcheck_list=[]
       for file_check in list_files_plots:
           f_checkings= open(file_check, 'r')
           plotcheckf = yaml.load(f_checkings)
           plotcheck_list.append(plotcheckf)
           f_checkings.close()

       for plotcheck in plotcheck_list:
           if varid in plotcheck.keys():
              #l_yes_plotted.append(varid)
              myname = new_name
              unit_stat = True    # take care of this. Dangerous.

              for mymethod in plotcheck[varid]['method'].keys():
                 if 'surface_map' in mymethod:
                    a1, a2, a3 = plots.surface_map(myname, varid, finfo=finfo, 
                                                   settings=mysettings, config=myconfig, 
                                                   plotinfo=plotcheck[varid]['method'][mymethod])


                 elif mymethod=='zonal_mean':
                    a1 = plots.zonal_mean(myname, varid,  finfo=finfo, 
                                                   settings=mysettings, config=myconfig, 
                                                   plotinfo=plotcheck[varid]['method'][mymethod])
                 else:
                    #l_non_plotted.append(direc)
                    print("-")


    return vars_checks, l_non_processed, l_yes_processed, l_uni_processed, l_yes_checked

def post_processing (myconfig, 
                     mysettings, 
                     directory, 
                     list_variables, 
                     lyears, 
                     post_t, 
                     study, 
                     expID,
                     freq, 
                     finfo=None,
                     model='IPSL'):

    """
    This function is doing the post-processing of the extracted files

    Note that it is not the best function from the programming point of view as the
    results depends strongly on the files found in the directory.

    - It should be important to add something that test if the files on the directories
      are the files expected from the settings file

    """

    def _sanity_check_fnames(lfiles):
        """
        Ensure that we have a set of filenames reasonable
        """
        for indx, item in enumerate(lfiles):
             if '-accum' in item or '_accum' in item:
                 del lfiles[indx]
        for item in lfiles:
             if len(item)!=len(lfiles[0]):
                 str1=('\n     -- Problem not all filenames have the same length: unusual')
                 myprint(str1, finfo=finfo)
                 if mysettings['safety']['check_files']==True:
                     exit()

        return lfiles

    def _accumulate(lfiles, post_t):
         #### -------------- Generate ACCUM files -----------------------------------
         if len(lfiles)==1:
            iyr_ini = 0
            iyr_end = 0

            first_yr  = lfiles[iyr_ini].split('_')[-2]
            last_yr   = lfiles[iyr_end].split('_')[-1].split('.')[0]
            stringyr  = lyears[0]
            if len(stringyr)>4 or len(stringyr)<2:
                 str1=('\n     -- Problem with year in filename. Not usual. Check')
                 myprint(str1, finfo=finfo)
                 exit()

            new_file_t= lfiles[iyr_ini].replace(stringyr+'0101',first_yr)
            new_file  = new_file_t.replace('_'+stringyr+'1231','_'+last_yr+'-accum')
            varid = new_file.split('/')[-1].split('_')[0]
            #print('check this', new_file,varid)

         if len(lfiles)>1:
            iyr_ini =  1 ### Reason ?? -> we discard first file?? -> it is spining
            iyr_end = -1
            first_yr  = lfiles[iyr_ini].split('_')[-2]
            last_yr   = lfiles[iyr_end].split('_')[-1].split('.')[0]
            str_files = ' '.join(lfiles[iyr_ini::])
            stringyr  = lyears[1]
            new_file_t= lfiles[iyr_ini].replace(stringyr+'0101',first_yr)
            new_file  = new_file_t.replace('_'+stringyr+'1231','-'+last_yr+'-accum')
            varid = new_file.split('/')[-1].split('_')[0]

            # Concatenation of all files -------------------------------------------------
            if post_t=='all':
                # print('is concatenating ----')
                # this create the -accum file --------------------------------------------
                command_nco ='ncrcat -O '+str_files +' '+ new_file
                os.system(command_nco)

                varid = new_file.split('/')[-1].split('_')[0]
                if varid!=vname:
                    varid=vname

         return new_file, varid

    # This open a file with the list of test per variable. I call it cmip6 because I use
    # cmip6 names. The file name it can be stored in config.yaml

    list_files_checks = glob.glob('etc/checks/check_*.yaml')
    check_list=[]
    for file_check in list_files_checks:
        f_checkings= open(file_check, 'r')
        checkl = yaml.load(f_checkings)
        check_list.append(checkl)
        f_checkings.close()
        
    # We define here the columns of the info-checks file ------------
    dfcolumns=['variable', 'method', 'value', 'units', 'testing-file']

    vars_checks     = []
    l_non_processed = []
    l_yes_processed = []
    l_uni_processed = []
    l_yes_checked   = []


    write_process_vars=False

    for ivname, vname in enumerate(list_variables):
         str1=('\n    ... Post-Processing of files related with '+vname)
         myprint(str1, finfo=finfo)

         lfiles = sorted(glob.glob(directory+vname+'/'+vname+'*.nc'))

         ## we select only those files that begins with the variable name, as they were 
         ## generated on the pre-processing step.
         ## We also ensure that we are not including aggregations of several files
         ## and we remove those files with the '_accum' string on the name from the list

         lfiles          = _sanity_check_fnames(lfiles)

         if post_t == 'all':
            new_file, varid = _accumulate(lfiles, post_t)
            out = check_var_netcdf(vars_checks, l_non_processed, l_yes_processed, l_uni_processed,
                                   l_yes_checked, check_list,  directory, vname, mysettings, study, expID, varid, new_file, 
                                   post_t, freq=freq, finfo=finfo)
            vars_checks, l_non_processed, l_yes_processed, l_uni_processed, l_yes_checked = out

         if post_t=='year' or 'AEROCOM' in expID:
            for fname in lfiles:
                 new_file, varid = _accumulate([fname], post_t)
                 out = check_var_netcdf(vars_checks, l_non_processed, l_yes_processed, l_uni_processed,
                                   l_yes_checked, check_list, directory, vname, mysettings, study, expID, varid, fname, 
                                   post_t, freq=freq, finfo=finfo, myconfig=myconfig)
                 vars_checks, l_non_processed, l_yes_processed, l_uni_processed, l_yes_checked = out

    ## SUMMARIZE INFORMATION POST-PROCESSSING ==========================================================================

    non_cmip6_units =  [y for i,y in enumerate(l_yes_processed) if l_uni_processed[i]==False]
    yes_cmip6_units =  [y for i,y in enumerate(l_yes_processed) if l_uni_processed[i]==True]

    str1=('\n\n  Processed variables ...... ['+str(len(l_yes_processed))+'] '+ str(l_yes_processed)+ '\n')
    str2=('    without CMIP6 units .... '+ str([y for i,y in enumerate(l_yes_processed) if l_uni_processed[i]==False])+ '\n')
    str3=('    with    CMIP6 units .... '+ str([y for i,y in enumerate(l_yes_processed) if l_uni_processed[i]==True])+'\n')
    with open('vars_cf_cmip6.txt', 'w') as fcmip6:
        for y in yes_cmip6_units:
            fcmip6.write(y+'\n')
    str4=('  Non processed variables ... ['+str(len(l_non_processed))+'] '+ str(l_non_processed)+'\n')
    myprint(str1+str2+str3+str4, finfo=finfo)
    myprint('\n    Checked variables ... ['+str(len(l_yes_checked))+'] '+ str(l_yes_checked)+ '\n', finfo=finfo)
    l_non_checked = set(l_yes_processed)-set(l_yes_checked)
    myprint(  '  ... processed non checked ... ['+str(len(l_non_checked))+'] '+ str(l_non_checked)+ '\n', finfo=finfo)

    if write_process_vars==True:
       pd.set_option('display.width', 140)
       if freq=='hr':
          previ_process_file = mysettings[study]['process']['hrs']
       else:
          previ_process_file = mysettings[study]['process'][freq]
       cmip6_process_file = previ_process_file.replace('.process', '_cmip6.process')
       noncmip6_process_file = previ_process_file.replace('.process', '_noncmip6.process')

       df = pd.read_table(previ_process_file, skiprows=0, header=0, delim_whitespace=True)
       cmip6_df = df[df['varname-out'].isin(yes_cmip6_units)]
       cmip6_df.to_csv('tempf', sep=' ', index=False)
       os.system(('column -t tempf > '+ cmip6_process_file))
       non_df = df[df['varname-out'].isin(non_cmip6_units+l_non_processed)]
       non_df.to_csv('tempf', sep=' ', index=False)
       os.system(('column -t tempf > '+ noncmip6_process_file))
       os.remove('tempf')

    for varid in l_non_checked:
         vars_checks.append({'variable': varid, 'method': 'not-set-up', 'value': -9.99, 'units':'', 'testing-file':''})
    test_info = pd.DataFrame(vars_checks, columns=dfcolumns )

    return test_info

def find_tableID(varID, freq='mon', origin='ATM', dirtables=''):
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
    """ This function creates the new filename of the netcdf with the diagnostics 
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
        print('         ---> variable: * ',varID, ' not found in CMIP6-cmor-tables\n')
        return False, False

    tableID      = table_name.split('/')[-1].split('.')[0].split('CMIP6')[1][1::]
    sourceID     = mysettings[study]['expCVs']['sourceID']
    experimentID = mysettings[study]['expID']
    memberID     = mysettings[study]['expCVs']['memberID']
    gridlabel    = mysettings[study]['expCVs']['gridlabel']
    timerange    = varfile.split('_')[-1]


    #print(varID, 'here')

    if mysettings[study]['expKIND'].upper() in ['CMIP6','CRESCENDO','JASPERKOK','CLIMDO']:
        # for CMIP6 or CRESCENDO
        save_name = '_'.join([varID ,tableID, sourceID, experimentID, memberID,
                              gridlabel, timerange
                              ])
        print(save_name, timerange)

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
    myprint('       ... final filename will be: '+ save_name, finfo=finfo)
    return save_name, info_var


def add_global_attributes(mysettings, ncname, newname, study, info_var, finfo=None):
    """

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
        elif 'lat' in dataset.coords and 'lon' in dataset.coords and 'lev' in dataset.coords:
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

    if 'airmassiiii' in finalname or 'ec550iii' in finalname:
        print('')
    else:
        dataset.to_netcdf(finalname, mode='w')
        dataset.close()
    #else:
    #    shutil.copy2(ncname,finalname)
    return unit_status, finalname


 
