#!/usr/bin/python

""" fetches a one day worthy of data, concatenate them and do boxcar filtering """

import sys
sys.path.append("/home/muhammad/Documents/Important/midlat_convection/dopsearch_py")
sys.path.append("/home/muhammad/Documents/Important/sd-data_update_indices")
from  dopsearch import prepare_file
from functions import scp_to_sddata, remove_tmp_dir, create_tmp_dir, get_password
import datetime as dt
from davitpy.pydarn.sdio.fetchUtils import fetch_local_files
import os
import multiprocessing as mp
#davitpy.rcParams['verbosity'] = "debug"

# input parameters

# create center datetimes 
""" NOTE: Do a year of data for a certain ftype at a time
Also, set the channel value correctly
"""
yr = 2014
ftype = "fitacf"
#ftype = "fitex"
sctr_time = dt.datetime(yr,1,1)
#sctr_time = dt.datetime(yr,9,13)
ectr_time = dt.datetime(yr+1,1,1)
#ectr_time = dt.datetime(yr,2,28)
num_days = (ectr_time - sctr_time).days + 1
dt_range = [dt.timedelta(days=i)+sctr_time for i in xrange(num_days)]

#rad_list = ["bks", "wal", "fhe", "fhw", "cve", "cvw"]
#rad_list = ["bks",  "wal", "fhe"]
#rad_list = ["fhw", "cve", "cvw"]
rad_list = ["cve", "cvw"]
#rad_list = ["tig", "unw", "bpk", "hok", "hkw"]
#rad_list = ["bpk", "hok", "hkw"]
#rad_list = ["hok", "hkw"]
#rad_list = ["ade", "adw"]
channel = None
#channel = '.'
scr = "local"
localdirfmt = "/sd-data/{year}/{ftype}/{radar}/"
fnamefmt = ['{date}.{hour}......{radar}.{channel}.{ftype}', '{date}.{hour}......{radar}.{ftype}']

# get sd-data password
#sddata_password = get_password()

for rad in rad_list:
    # create tmpdir to store the data
    tmpdir = create_tmp_dir(tmp_dir=rad+"_tmp")

    localdict = {"ftype" : ftype, "radar" : rad, "channel" : channel}
    for ctr_dtm in dt_range:

        print "start boxcarfiltering for data on ", ctr_dtm
        # concat 1 day of data, do boxcar filtering 
        ffname = prepare_file(ctr_dtm, localdirfmt, localdict, tmpdir,
                              fnamefmt, oneday_file_only=True)
        print "end boxcarfiltering for data on ", ctr_dtm
        print "generated ", ffname

        # remove the original not filtered data
        if ffname is not None:
            os.system("rm " + ffname[:-1])
        else:
            continue

    # scp to sd-data
#    sddata_file_path = "/sd-data/backup/muhammad/muhammad_data/filtered_concat3days/" + \
#                        str(yr) + "/" + ftype + "/" + rad + "/"
#    scp_to_sddata(tmpdir, sddata_file_path, sddata_password)

#    # delete tmpdir
#    remove_tmp_dir(tmpdir)  

