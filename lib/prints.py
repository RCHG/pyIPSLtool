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


def myprint(mystr, finfo=None):
    """myprint: just a python to manage the log prints

    :param mystr:
    :param finfo:
    """

    if finfo!=None:
        finfo.write(mystr+'\n')
    return



