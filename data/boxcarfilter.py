# load modules
import sys
sys.path.append("../dopsearch_py")
sys.path.append("/home/muhammad/Documents/Important/sd-data_update_indices")
from  dopsearch import prepare_file
import datetime as dt
from davitpy.pydarn.sdio.fetchUtils import fetch_local_files
import os
import multiprocessing as mp
import logging
#davitpy.rcParams['verbosity'] = "debug"

def do_boxcarfiltering(sctr_day, ectr_day, rad_list, ftype, channel, n_jobs=None):

    """ fetches, concatenates and does boxcar median filtering for data
    in the time iterval between sctr_day, ectr_day for a given radar.
    Multiprocessing is implemented in the code,
    where a day of data from a radar is defined as a process.

    Parameters
    ----------
    sctr_day : datetime.datetime
        Start day
    ectr_day : datetime.datetime
        End day
    rad_list : list
        Three-letter radar codes
    ftype : str
        Data file type. e.g., "fitacf", "fitex"
    channel : str
        radar channel. e.g., "a", "b", "c", "d", or "." which is all
    n_jobs : int or None
        Number of jobs that run in parallel.
        Default to None, in which case all the CPUs but one will used.

    Return
    ------
    Nothing
        
    """

    
    # calculate number of days 
    sctr_day = dt.datetime(sctr_day.year, sctr_day.month, sctr_day.day)
    ectr_day = dt.datetime(ectr_day.year, ectr_day.month, ectr_day.day)
    num_days = (ectr_day - sctr_day).days + 1  # includes the end day
    dt_range = [dt.timedelta(days=i)+sctr_day for i in xrange(num_days)]


    # Define an output queue
    output = mp.Queue()

    # number of jobs to be run in parallel
    if not n_jobs:
        # get the number of CUPs
        n_jobs = mp.cpu_count() - 1

    # loop throughs the radars
    for rad in rad_list:

        # create tmpdir to store the data
        rad_tmp_dir = rad + "_tmp"
        os.system("mkdir -p " + rad_tmp_dir)
        rad_tmp_dir = os.getcwd() + "/"+ rad_tmp_dir + "/"
        print rad_tmp_dir + " is created"


        # cteate tmpdirs, one for each n_jobs
        tmp_dirs = []
        for j in range(n_jobs):

            # create tmpdirs to store the data, one tmpdir for each process
            tmp_dir = rad + "_tmp" + "_" + str(j)
            os.system("mkdir -p " + tmp_dir)
            tmp_dir = os.getcwd() + "/"+ tmp_dir + "/"
            print tmp_dir + " is created"
            tmp_dirs.append(tmp_dir)
        
        # iter through the days, one day at a time
        i = 0
        while 1:

            # run n_jobs in parallel
            dts_tmp = dt_range[i*n_jobs:(i+1)*n_jobs]
            i = i + 1
            if len(dts_tmp) == 0:
                break
            else:
                procs = []
                # send jobs in parallel
                for j, ctr_dtm in enumerate(dts_tmp):
                   p = mp.Process(target=worker, args=(rad, ctr_dtm, ftype, channel, tmp_dirs[j]))
                   procs.append(p)
                   p.start()

                # exit the completed processes
                for p in procs:
                    p.join()

        # move processed data from tmpdirs of all processes to a single tmpdir
        for j in range(n_jobs):
            os.system("mv " + tmp_dirs[j] + "* " + rad_tmp_dir)

        # remove tmpdirs
        os.system("rm -rf " + rad + "_tmp_*")
        print "tmpdirs have been deleted"

    return

def worker(rad, ctr_dtm, ftype, channel, tmp_dir):
    """ A worker function fetches a one day worthy of data, 
    concatenate them and do boxcar median filtering.
    
    Parameters
    ----------
    rad : str
	Three-letter radar code
    ctr_dtm : datetime.datetime
    ftype : str
        Data file type. e.g., "fitacf", "fitex"
    channel : str
        radar channel. e.g., "a", "b", "c", "d", or "." which is all.
    tmp_dir : str
        temprory folder to store the processed data

    Return
    ------
    Nothing
            
    """

    print "start boxcarfiltering for data on ", ctr_dtm
    # concat 1 day of data, do boxcar filtering 
    scr = "local"
    localdict = {"ftype" : ftype, "radar" : rad, "channel" : channel}
    localdirfmt = "/sd-data/{year}/{ftype}/{radar}/"
    fnamefmt = ['{date}.{hour}......{radar}.{channel}.{ftype}',
                '{date}.{hour}......{radar}.{ftype}']

    ffname = prepare_file(ctr_dtm, localdirfmt, localdict, tmp_dir,
                          fnamefmt, oneday_file_only=True)
    print "end boxcarfiltering for data on ", ctr_dtm
    print "generated ", ffname

    # remove the original not filtered data
    if ffname is not None:
        os.system("rm " + ffname[:-1])
    else:
        pass

    return

def main():
    """Execute the codes above"""

    # input parameters
    # create center datetimes 
    """ NOTE: Set the channel value correctly. """

    yrs = [2015, 2016]
    #yrs = [2016]
    for yr in yrs:
        ftype = "fitacf"
        #ftype = "fitex"
        sctr_day = dt.datetime(yr,1,1)
        ectr_day = dt.datetime(yr+1,1,1)

        #rad_list = ["wal", "fhe", "fhw", "cve", "cvw"]
        #rad_list = ["bks",  "wal", "fhe"]
        #rad_list = ["fhw", "cve", "cvw"]
        #rad_list = ["bks"]
        rad_list = ["tig", "unw"]
        #rad_list = ["bpk", "hok", "hkw"]
        #rad_list = ["hok", "hkw"]
        #rad_list = ["ade", "adw"]
        channel = None
        #channel = '.'

        do_boxcarfiltering(sctr_day, ectr_day, rad_list, ftype, channel, n_jobs=20)

    return

if __name__ == "__main__":
    main()


