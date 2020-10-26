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

def memory_usage_psutil():
    # return the memory usage in MB
    import psutil
    return psutil.virtual_memory()


def myprint(mystr, finfo=None):
    print(mystr)
    if finfo!=None:
        finfo.write(mystr+'\n')
    return

def post_processing(mysettings, directory, list_variables, lyears, post_t, study, expID,
                    freq, finfo=None):
    '''
    This function is doing the post-processing of the extracted files

    Note that it is not the best function from the programming point of view as the
    results depends strongly on the files found in the directory.

    - It should be important to add something that test if the files on the directories
      are the files expected from the settings file

    '''

    gc.collect()  # to ensure that we are not collecting too much garbage.

    ### Add a new config file where the name of the files with check are added
    ### then we have to make a kind of loop for all files

    f_checkings= open('etc/check_cmip6_variables.yaml', 'r')
    cmip6check = yaml.load(f_checkings)

    # We define here the columns of the info-checks file ------------
    dfcolumns=['variable','method','value', 'units', 'testing-file']

    vars_checks     = []
    l_non_processed = []
    l_yes_processed = []
    l_uni_processed = []
    l_yes_checked   = []

    for idirec, direc in enumerate(list_variables):
         str1=('\n    ... Post-Processing of files related with '+direc)
         myprint(str1, finfo=finfo)
         lfiles = sorted(glob.glob(directory+direc+'/'+direc+'*.nc'))
         #print(lfiles)
         ## we select only those files that begins with the variable name, as they were 
         ## generated on the pre-processing step.
         ## We also ensure that we are not including aggregations of several files
         ## and we remove those files with the '_accum' string on the name.
         for indx, item in enumerate(lfiles):
             if '-accum' in item or '_accum' in item:
                 del lfiles[indx]
         for item in lfiles:
             if len(item)!=len(lfiles[0]):
                 str1=('\n     -- Problem not all filenames have the same length: unusual')
                 myprint(str1, finfo=finfo)

                 if mysettings['safety']['check_files']==True:
                     exit()

         if len(lfiles)==1:
            first_yr  = lfiles[0].split('_')[-2]
            last_yr   = lfiles[0].split('_')[-1].split('.')[0]
            stringyr  = lyears[0]
            new_file_t= lfiles[0].replace(stringyr+'0101',first_yr)
            new_file  = new_file_t.replace('_'+stringyr+'1231','-'+last_yr+'-accum')
            varid = new_file.split('/')[-1].split('_')[0]
            #print('check this', new_file,varid)

         if len(lfiles)>1:
            first_yr  = lfiles[1].split('_')[-2]
            last_yr   = lfiles[-1].split('_')[-1].split('.')[0]
            str_files = ' '.join(lfiles[1::])
            stringyr  = lyears[0]
            new_file_t= lfiles[0].replace(stringyr+'0101',first_yr)
            new_file  = new_file_t.replace('_'+stringyr+'1231','-'+last_yr+'-accum')
            #print(new_file)
            varid = new_file.split('/')[-1].split('_')[0]

            # Concatenation of all files -------------------------------------------------
            if post_t=='all':
                commandnco ='ncrcat -O '+str_files +' '+ new_file
                os.system(commandnco)
                varid = new_file.split('/')[-1].split('_')[0]
                if varid!=direc:
                    varid=direc
                good_name, info_var = create_filename(mysettings, study, expID, varid,
                                                      new_file, freq=freq, finfo=finfo)
                if good_name!=False:
                    save_name = directory+direc+'/'+good_name
                    unit_stat, finalname = add_global_attributes(mysettings, new_file, 
                                                                 save_name, study, info_var)
                    l_yes_processed.append(direc)
                    l_uni_processed.append(unit_stat)

                    if varid in cmip6check.keys():
                       l_yes_checked.append(varid)
                       for mymethod in  cmip6check[varid]['method'].keys():
                           if mymethod=='tendency_mass':
                                a1, a2, a3 = check.tendency_mass(finalname, varid, 
                                                      'etc/area_grid.nc', 'area', finfo=finfo)
                                vars_checks.append({'variable': a1,
                                                    'method': mymethod,
                                                    'value': a3, 'units':a2, 
                                                    'testing-file':finalname.split('/')[-1]})

                           if mymethod=='global_mean':
                                a1, a2, a3, a4= check.global_mean(finalname, varid,
                                                               cmip6check[varid], finfo=finfo)
                                vars_checks.append({'variable': a1,
                                                    'method': mymethod,
                                                    'value': a3, 'units':a2,
                                                    'testing-file':finalname.split('/')[-1]})

                           if mymethod=='global_mean_atlev':
                                for level in cmip6check[varid]['method'][mymethod]['levels']:
                                    a1, a2, a3 , a4= check.global_mean_atlev(finalname, varid, 
                                                  cmip6check[varid], level=level, finfo=finfo)
                                    vars_checks.append({'variable': a1,
                                                        'method': 'glo mean at'+a4,
                                                        'value': a3, 'units':a2,
                                                        'testing-file':finalname.split('/')[-1]})

                           if mymethod=='total_load':
                                a1, a2, a3  = check.total_load(finalname, varid,
                                                               cmip6check[varid], area_fname='data/area_grid.nc',
                                                               area_varid='area', finfo=finfo)
                                vars_checks.append({'variable': a1,
                                                    'method': mymethod,
                                                    'value': a3, 'units':a2,
                                                    'testing-file':finalname.split('/')[-1]})
                else:
                    l_non_processed.append(direc)

         if post_t=='year' or 'AEROCOM' in expID:
                #print(expID, '++++++++++++++++++++++++++++++++++++++++')
                for fname in lfiles:
                    good_name, info_var = create_filename(mysettings, study, expID, varid,
                                                          fname.replace('1_20','1-20'),
                                                          freq=freq, finfo=finfo)
                    if good_name!=False:
                        save_name = directory+direc+'/'+good_name
                        unit_stat, finalname = add_global_attributes(mysettings, fname, save_name, study, info_var)
                        l_yes_processed.append(direc)
                        l_uni_processed.append(unit_stat)
                        if varid in cmip6check.keys():
                            l_yes_checked.append(varid)
                            for mymethod in  cmip6check[varid]['method'].keys():
                               if mymethod=='tendency_mass':
                                    a1, a2, a3 = check.tendency_mass(finalname, varid, 
                                                          'etc/area_grid.nc', 'area', finfo=finfo)
                                    vars_checks.append({'variable': a1,
                                                        'method': mymethod,
                                                        'value': a3, 'units':a2, 
                                                        'testing-file':finalname.split('/')[-1]})

                               if mymethod=='global_mean':
                                    a1, a2, a3, a4= check.global_mean(finalname, varid,
                                                                   cmip6check[varid], finfo=finfo)
                                    vars_checks.append({'variable': a1,
                                                        'method': mymethod,
                                                        'value': a3, 'units':a2,
                                                        'testing-file':finalname.split('/')[-1]})

                               if mymethod=='global_mean_atlev':
                                    for level in cmip6check[varid]['method'][mymethod]['levels']:
                                        a1, a2, a3 , a4= check.global_mean_atlev(finalname, varid, 
                                                      cmip6check[varid], level=level, finfo=finfo)
                                        vars_checks.append({'variable': a1,
                                                            'method': 'glo mean at'+a4,
                                                            'value': a3, 'units':a2,
                                                            'testing-file':finalname.split('/')[-1]})

                    else:
                        l_non_processed.append(direc)

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

    for varid in l_non_checked:
        vars_checks.append({'variable': varid, 'method': 'not-set-up', 'value': -9.99, 'units':'', 'testing-file':''})
    test_info = pd.DataFrame(vars_checks, columns=dfcolumns )
    os.remove('tempf')
    return test_info

def find_tableID(varID, freq='mon', origin='ATM'):
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
    Filename will have the structure:
        varID_tableID_sourceID_experimentID_memberID_gridlevel_timerange

    varID  -------> is the variable name: varname
    tableID ------> depends on variable AERday, Amon etc...
    sourceID -----> IPSL-LMDZORINCAv6
    experimentID -> CRESCWP3PD-amip
    memberID -----> v5-r1i1p1f1  (v5 refers to our internal version)
    gridlabel ----> gr           (refers to the original grid but we provide at preslev)
    timerange ----> 20000101-20121231 for all cases.
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

def directory_structure(study, mysettings, clean=False, create=False, finfo=None):

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

def _get_var_checkfiles(outname, yr, study, freqsource, module, timeout, exp,
                        extra, storepath, dirMO, dirDA, dirHF):
        #dic_struc = { 'ATM':{'MO':('_1M_histmth.nc','_mon'), 'DA':('_1M_histmth.nc','_mon'),
        #                     'HF':('_HF_histh1f.nc','_1hr')}

        f_tail   = '_IPSL_'+study+'_'+yr+'0101_'+yr+'1231'+'.nc'
        bpath = storepath+exp
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

        if module=='CHM':
            if freqsource=='MO':
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

        return input_f, outfile

def solve_time(outfile, outfilen, varname='', out=''):

    commandnco = ('ncrename -O -d time_counter,time -v time_counter,time ' + outfile
                  + ' ' + outfilen + out)
    os.system(commandnco)
    commandnco = 'ncrename -O -v .time_counter_bnds,time_bnds '+outfile+' '+outfilen+out
    os.system(commandnco)
    commandnco = 'ncrename -O -d .presnivs,pres -v .presnivs,pres '+outfile+' '+outfilen+out
    os.system(commandnco)
    commandnco = 'ncatted -O -a standard_name,pres,o,c,air_pressure '+outfile+' '+outfilen+out
    os.system(commandnco)
    commandnco = 'ncatted -O -a bounds,time,o,c,time_bnds '+outfile+' '+outfilen+out
    os.system(commandnco)
    commandnco = 'ncks -O -x -v time_centered '+outfile+' '+outfilen+out
    os.system(commandnco)
    commandnco = 'ncks -O -x -v time_instant '+outfile+' '+outfilen+out
    os.system(commandnco)

    if varname is not '':
        commandnco = 'ncatted -O -a coordinates,'+varname+',d,, '+outfile+' '+outfilen+out
        os.system(commandnco)
        #print(commandnco, varname)


    return

def get_var(mysettings, year, varname, freqsource, module, timeout,
            outname, expcase, study, extra='',index=1,finfo=None):
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

    input_f, outfile= _get_var_checkfiles(outname, year, study, freqsource, module, timeout,
                                            expcase, extra, storepath, dirMO, dirDA, dirHF)
    
    str1=('    Extracting ('+str(index).rjust(4)+'): '+varname.ljust(20)+
          ' and save as '+outname.ljust(20)+' .... for year '+
          year + '  ... ['+freqsource+' '+module+']')

    myprint(str1, finfo=finfo)


    commandcdo =('cdo  --silent selvar,' + varname + ' ' + input_f + ' ' + outfile
                + ' &> delfiles/del' + varname + freqsource + module)
    #print(commandcdo)
    os.system(commandcdo)

    solve_time(outfile, outfile, varname=varname, out=' &>> delfiles/del'+
               varname+freqsource+module)

    if outname!=varname:
      commandcdo ='cdo --silent -O chname,'+varname+','+outname+' '+outfile+' '+outfile+'x'
      #print(commandcdo)
      os.system(commandcdo)
      os.rename(outfile+'x', outfile)

    return

def check_units(varname_list, vars_dataset):
    lunits = []
    for varn in varname_list:
        lunits.append(vars_dataset[varn].attrs['units'])

    return lunits


def scale_variable(vardataset, d_scale_vars, outname, varname, outfile, finfo=None):

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
 

def new_variable_altsur(vardataset, d_new_vars, outname, varname, outfile, finfo=None):

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
        if ivar==0:
            var_sum_alt = vardataset[var_alt].sum(dim='presnivs')
        else:
            var_sum_alt = (var_sum_alt +
             vardataset[var_alt].sum(dim='presnivs')*d_new_vars[outname]['fac_altitud'][ivar])

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
            expcase, study, extra='',index=1, finfo=None):

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

    input_f, outfile  = _get_var_checkfiles(outname, year, study, freqsource, module,
                                               timeout, expcase, extra, storepath,
                                               dirMO, dirDA, dirHF)
    vardataset = xr.open_dataset(input_f, chunks={'time_counter': 10})  # we assume that all variable are in same file.

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

        new_variable_altsur(vardataset, d_new_vars, outname, varname, outfile, finfo)

    solve_time(outfile, outfile, varname=varname, out=' &>> delfiles/del'
               +varname+freqsource+module)

    return

def lev_var(mysettings, year, varname, freqsource, module, timeout, outname,
            expcase, study, extra='', lev='top', index=1, finfo=None):

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

    input_f, outfile  = _get_var_checkfiles(outname, year, study, freqsource, module,
                                               timeout, expcase, extra, storepath,
                                               dirMO, dirDA, dirHF)

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
    solve_time(outfile, outfile, varname=varname, out=' &>> delfiles/del'+
               varname+freqsource+module)

    return

def read_process_file(mysettings, process_file, year, expcase, study, finfo=None):
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
                        expcase, study, extra=row['extra-info-source'], index=index, finfo=finfo)
            if row['method']=='get_var':
                get_var(mysettings, year, row['varname-source'], row['freq-source'],
                        row['module-IPSL'], row['freq-specific'], row['varname-out'],
                        expcase, study, extra=row['extra-info-source'], index=index, finfo=finfo)
            if row['method']=='lev_var':
                lev_var(mysettings, year, row['varname-source'], row['freq-source'],
                        row['module-IPSL'], row['freq-specific'], row['varname-out'],
                              expcase, study, extra=row['extra-info-source'], index=index, finfo=finfo)
    return list_vars

def process_files(mysettings, year, expcase, study, proc_file, finfo=None):
      if str(proc_file)!='none':
          myprint('    ... aprox. number variables: '+ str(sum(1 for line in open(proc_file))-1), finfo=finfo)
          processed_vars = read_process_file(mysettings, proc_file, year, expcase, study, finfo=finfo)
      else:
          myprint('    ... No processing ... [process-file==None]', finfo=finfo)

          processed_vars=[]

      return processed_vars


