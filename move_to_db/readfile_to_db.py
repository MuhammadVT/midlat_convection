'''
Created on Aug. 5, 2016

Muhamamd
'''

import pdb

class read_file_to_db(object):
    """ A class that holds the boxcar median filtered data after reading them from a file.
    It writes the data into a sqlite db using its move_to_db method."""

    def __init__(self, rad, ctr_date, ftype="fitacf", params=["velocity"], ffname=None):

        """ 
        Parameters
        ----------

        rad : str
            Three-letter code for a rad
        ctr_date : datetime.datetime
            a full day for which data are to be read. 
        ftype : str, default to "fitacf"
            SuperDARN file type. Valid inputs are "fitacf", "fitex"
        params : list
            NOTE: works for params=["velocity"] only
        ffname : string, default to None
            Full path of a file to be read. if ffname is not set to None, 
            ffname will be constructed.

        Returns
        -------
        read_file_to_db object 
            A dict of dicts in the form of {bmnum:dict}, which holds a 
            one day worth of data from all beams of a certain radar.
            
        """

        import sys
        sys.path.append("../")
        from dopsearch_py.dopsearch import read_file
        import datetime as dt
        
        # build the attributes
        self.rad = rad
        self.ftype = ftype
        self.ctr_date = dt.datetime(ctr_date.year, ctr_date.month, ctr_date.day)
        if ffname is None:
            # construct ffname (file full path)
            self.ffname = self._construct_filename()
        else:
            self.ffname = ffname

        # create stm (start time) and etm (end time) 
        stm = self.ctr_date
        # add two minute to the etm to read the last record from
        # the boxcar filtered concatenated fitacf(ex) files
        etm = self.ctr_date + dt.timedelta(days=1) + dt.timedelta(minutes=2)

        # read data from file 
        # Note: data_from_db argument has to be False
        self.data = read_file(self.ffname, rad, stm, etm, params,
                              ftype=self.ftype, data_from_db=False,
                              plotrti=False)

    def _construct_filename(self, basedir="../data/"):
        """ constructs filename with full path for a file of interest

        Parameters
        ----------
        basedir : str
            Relative path for data files

        Returns
        -------
        str
            Filename with its full path
        
        """

        import datetime as dt

        # create stm (start time) and etm (end time) 
        stm = self.ctr_date
        etm = self.ctr_date + dt.timedelta(days=1)

        # consturc a file name with full path
        ffname = stm.strftime("%Y%m%d.%H%M%S") + "." + \
                 etm.strftime("%Y%m%d.%H%M%S") + "." + \
                 self.rad + "." + self.ftype + "f"
        ffname = basedir + self.rad + "/" + ffname

        return ffname

    def move_to_db(self, conn):
        """ writes the data into sqlite db

        Parameters
        ----------
        conn : sqlite3.connect

        Returns
        -------
        Nothing
        """

        import json 
        import sqlite3

        cur = conn.cursor()

        # loop through all the beams
        for bmnum in self.data.keys():
            data_dict = self.data[bmnum]

            # create a table in sqlite db
            table_name = self.rad + "_bm" + str(bmnum)
            command = "CREATE TABLE IF NOT EXISTS {tb} (\
                      vel TEXT, rsep REAL, frang REAL, bmazm REAL,\
                      slist TEXT, gsflg TEXT,\
                      datetime TIMESTAMP PRIMARY KEY)".format(tb=table_name)
            cur.execute(command)

            # loop through each scan time, usually 2 minutes,
            # and write the data into table_name in the sqlite db
            for i in xrange(len(data_dict['datetime'])):
                command = "INSERT OR IGNORE INTO {tb} (vel, rsep, frang, bmazm,\
                            slist, gsflg, datetime) VALUES (?, ?, ?, ?, ?, ?, ?)"\
                            .format(tb=table_name)
                cur.execute(command, (json.dumps(data_dict["vel"][i]), data_dict["rsep"][i],
                            data_dict["frang"][i], data_dict["bmazm"][i],
                            json.dumps(data_dict["slist"][i]),
                            json.dumps(data_dict["gsflg"][i]),\
                            data_dict["datetime"][i]))

        # commit the change, once at one day of data points
        conn.commit()

def main():
    """ Call the functions above. Acts as an example code.
    Multiprocessing has been implemented to do parallel computing. The unit process is for """
    
    import datetime as dt
    import sqlite3
    from month_to_season import get_season_by_month
    import multiprocessing as mp

    # input parameters
    #season = "winter"
    #season = "summer"
    season = "equinox"

    # set the time interval for each season
    # NOTE: The code processes two years (full years specified by year_1 and year_2) of 
    #data at a time.
    year_0 = 2014; year_1 = 2015; year_2 = 2016; year_3 = 2017
    #year_0 = 2010; year_1 = 2011; year_2 = 2012; year_3 = 2013

    # run the code for the following radars in parallel
    #rad_list = ["bks", "wal", "fhe", "fhw", "cve", "cvw"]
    #rad_list = ["bks", "wal", "fhe"]
    #rad_list = ["fhw", "cve", "cvw"]
    rad_list = ["tig", "unw"]
    #rad_list = ["unw"]
    channel = None
    params=['velocity']
    ftype = "fitacf"
    #ftype = "fitex"
    ffname = None

    # make a db connection
#    season =  get_season_by_month(sdate.month)
    baseLocation = "../data/sqlite3/" + season + "/"

    # define the staring and ending dates for a given season
    # NOTE : winter (Nov - Feb), summer (May - Sep), equinox (Mar - Apr, Oct - Nov)
    if season == "winter":
        stms = [dt.datetime(year_0,12,31), dt.datetime(year_1,10,31), dt.datetime(year_2,10,31)]
        etms = [dt.datetime(year_1,3,1), dt.datetime(year_2,3,1), dt.datetime(year_3,1,1)]

    if season == "summer":
        stms = [dt.datetime(year_1,4,30), dt.datetime(year_2,4,30)]
        etms = [dt.datetime(year_1,9,1), dt.datetime(year_2,9,1)]

    if season == "equinox":
        stms = [dt.datetime(year_1,2,28), dt.datetime(year_1,8,31), dt.datetime(year_2,2,29), dt.datetime(year_2,8,31)]
        etms = [dt.datetime(year_1,5,1), dt.datetime(year_1,11,1), dt.datetime(year_2,5,1), dt.datetime(year_2,11,1)]

    # loop through the time intervals
    for dd in range(len(stms)):
        sdate = stms[dd]
        edate = etms[dd]
        num_days = (edate - sdate).days + 1
        dtm_range = [sdate + dt.timedelta(days=i) for i in xrange(num_days)]

        # loop through the radars
        for rad in rad_list:

            # make a db connection
            dbName = rad + "_" + ftype + ".sqlite"
            conn = sqlite3.connect(baseLocation + dbName)

            # loop through the days within a time interval:
            for ctr_date in dtm_range:

                # collect the data 
                t1 = dt.datetime.now()
                print "creating an object for " + rad + " for " + str(ctr_date)
                rf = read_file_to_db(rad, ctr_date, ftype=ftype, params=params, ffname=ffname)
                print "created an object for " + rad + " for " + str(ctr_date)
                if rf.data is not None:

                    # move data to db
                    rf.move_to_db(conn)
                    print ("object has been moved to db")

                t2 = dt.datetime.now()
                print ("creating and moving object to the db took " + str((t2-t1).total_seconds() / 60.)) + " mins\n"

            # close db connection
            conn.close()

    return


if __name__ == '__main__':
    main()

