"""
##############################################################################################

  Project   : CRESCENDO/AEROCOM
  Filename  : organize_yaml.py
  Author    : Ramiro Checa-Garcia
  email     : rcheca@lsce.ipsl.fr
  Purpose   : Reorganize variables for LMDzINCAOR experiment runs.

  Revision History ----------------------------------------------------------

  Date       Author     Ref    Revision

  2018-Apr   R.Checa           First version
  2018-Sep   R.Checa           Working as main with modules in lib
  2018-Oct   R.Checa           Implemented yaml
  2018-Nov   R.Checa           Final first version based on yaml

  TODO LIST:

##############################################################################################
"""

# External modules -----------------------------------------------------------------
import yaml                     # Manage the settings, new variables and check files
import platform                 # Just to print info of the computer used on calcul.
import datetime
import xarray as xr
import numpy  as np
import glob
import os
from   optparse import OptionParser

# Internal modules ----------------------------------------------------------------
from lib.organize import (directory_structure,
                          process_files, myprint)
from lib.postpro  import post_processing
from lib.tabulate import tabulate


#### Parsing function ------------------------------------------------------------------------
#

def opt_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-s", "--settings",
                      action="store",
                      dest="myfsettings",
                      default='settings.yaml',
                      help="Give the filename with the settings in YMAL format")
    parser.add_option("-c", "--config",
                      action="store",
                      dest="myfconfig",
                      default='settings/config.yaml',
                      help="Give the filename with the configuration in YMAL format")
    (options, args) = parser.parse_args()
    #print(options, args)
    #if len(args) != 1:
    #    parser.error("wrong number of arguments")

    return options, args

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

def manage_info(mysettings, study, screen=True):
    """manage_info

    :param mysettings:
    :param study:
    """
    datenow = datetime.datetime.now().strftime("%Y%m%d-%H%M%S_")           # datetime of calculation
    finfo   = open(myconfig['info_logs']+datenow+study+'.loginfo', 'w')    # we open a loginfo file.
    expcase = mysettings[study]['expIPSL']
    expID   = mysettings[study]['expID']
    lyears  = [str(y) for y in mysettings[study]['years']]
    pfiles  = mysettings[study]['process']
    post_t  = mysettings[study]['accum']

    # Show machine info 
    machine_info = platform.uname()._asdict()
    myprint('\n------------------------------------------------------------', finfo=finfo)
    myprint('--- Settings and info from: '+fname_set+ '\n',finfo=finfo)
    for key, value in machine_info.items():
        myprint('  '+key+':'+value, finfo=finfo)
    myprint('--- Configuration from:     '+fname_con+ '\n',finfo=finfo)

    # Show general study info from settings
    myprint('\n   Processing  ...............  to subdir: '+ study+
            ' ... with experiment name: '+ expID, finfo=finfo)
    myprint('    from '+myconfig['model_name']+'  experiment ..... '+ expcase, finfo=finfo)
    myprint('    for years ................ '+ str(lyears), finfo=finfo)
    myprint('    with processing files .... '+ str(pfiles)+ '\n\n', finfo=finfo)

    if screen==True:
      # Show machine info 
      machine_info = platform.uname()._asdict()
      myprint('\n------------------------------------------------------------', finfo=None)
      myprint('--- Settings and info from: '+fname_set+ '\n',finfo=None)
      for key, value in machine_info.items():
          myprint('  '+key+':'+value, finfo=finfo)
      myprint('--- Configuration from:     '+fname_con+ '\n',finfo=None)

      # Show general study info from settings
      myprint('\n   Processing  ...............  to subdir: '+ study+
              ' ... with experiment name: '+ expID, finfo=None)
      myprint('    from '+myconfig['model_name']+'  experiment ..... '+ expcase, finfo=None)
      myprint('    for years ................ '+ str(lyears), finfo=None)
      myprint('    with processing files .... '+ str(pfiles)+ '\n\n', finfo=None)


    return finfo, lyears, pfiles, post_t, expcase, datenow, expID



if __name__ == '__main__':

    #### GENERAL SETTINGS ------------------------------------------------------------------------
    #

    options, args = opt_parser()
    fname_set  = options.myfsettings
    fsettings  = open(fname_set, 'r')
    mysettings = yaml.load(fsettings)

    fname_con  = options.myfconfig
    fconfig    = open(fname_con, 'r')
    myconfig   = yaml.load(fconfig)

    work_dir = myconfig['workdir']

    list_studies = mysettings['studies']  # As much as possible here the program should be
                                          # functional (same result for same settings)

    # The program does a loop on the different studies requested on the settings

    for study in list_studies:
        finfo, lyears, pfiles, post_t, expcase, datenow, expID = manage_info(mysettings, study)

        # Checking file directories structure ---------------------------------------

        if mysettings['safety']['clean']==True:
            # This clean directory outputs
            directory_structure(study, mysettings, clean=True, finfo=finfo, 
                                model=myconfig['model_name'])
        # This create directory outputs if necessary and returns directory structure
        dmon, dday, dhrs = directory_structure(study, mysettings, create=True, 
                                               finfo=finfo, model=myconfig['model_name'])

        myprint('\n----- Created the directory structure ----- \n\n', finfo=None)
        # Pre-processing step -------------------------------------------------------
        if mysettings[study]['preprocess']==True:
           for year in lyears:
              myprint('\n----- '+year+' ------ day datasets '+study, finfo=finfo)
              myprint('\n----- '+year+' ------ day datasets '+study, finfo=None)
              daily_vars = process_files(myconfig,  mysettings, year, expcase, study, 
                                         pfiles['day'], finfo=finfo,  model=myconfig['model_name'])

              myprint('\n----- '+year+' ------ mon datasets '+study, finfo=finfo)
              myprint('\n----- '+year+' ------ mon datasets '+study, finfo=None)
              month_vars = process_files(myconfig, mysettings, year, expcase, study, 
                                        pfiles['mon'], finfo=finfo,  model=myconfig['model_name'])

              myprint('\n----- '+year+' ------ hrs datasets '+study, finfo=finfo)
              myprint('\n----- '+year+' ------ hrs datasets '+study, finfo=None)
              hours_vars = process_files(myconfig, mysettings, year, expcase, study, 
                                         pfiles['hrs'], finfo=finfo,  model=myconfig['model_name'])


        # Post-processing + checking values steps -----------------------------------
        if mysettings[study]['postprocess']==True:
            if len(month_vars)>=1 and month_vars!='none':
                test_vars_mon = post_processing(myconfig, mysettings, dmon, month_vars, lyears, post_t,
                                                study, expID, freq='mon', finfo=finfo,  model=myconfig['model_name'])
                with open('info_checks/'+datenow+study+'_mon_checks.txt', 'w') as ftable:
                    ftable.write(tabulate(test_vars_mon, tablefmt='simple'))

            if len(daily_vars)>=1 and daily_vars!='none':
                test_vars_day = post_processing(myconfig, mysettings, dday, daily_vars, lyears, post_t,
                                                study, expID, freq='day', finfo=finfo,  model=myconfig['model_name'])

                with open('info_checks/'+datenow+study+'_day_checks.txt', 'w') as ftable:
                    ftable.write(tabulate(test_vars_day, tablefmt='simple'))


            if len(hours_vars)>=1 and  hours_vars!='none':
                test_vars_hrs = post_processing(myconfig, mysettings, dhrs, hours_vars, lyears, post_t,
                                                study, expID, freq='hr' , finfo=finfo,  model=myconfig['model_name'])
                with open('info_checks/'+datenow+study+'_hrs_checks.txt', 'w') as ftable:
                    ftable.write(tabulate(test_vars_hrs, tablefmt='simple'))


        finfo.close()


