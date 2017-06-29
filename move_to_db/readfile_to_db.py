'''
Created on Aug. 5, 2016

Muhamamd
'''

import pdb

class read_file_to_db(object):
    """ A class that holds the boxcar filtered data. """

    def __init__(self, rad, ctr_date, ftype="fitacf", params=["velocity"], ffname=None):

        """ 
        Parameters
        ----------

        rad : str
            Three-letter code for a rad
        ctr_date : datetime.datetime
            a full day for which data are to be read. 
        ftype : str
            SuperDARN file type. Valid inputs are "fitacf", "fitex"
        params : list
            NOTE: works for params=["velocity"] only
        ffname : string
            if ffname is not set to None, ffname will be constructed.

        Returns : read_file_to_db object 
            A dict of dicts in the form of {bmnum:dict} that holds a 
            one day worth of data of a certain radar.
            
        """

        import sys
        sys.path.append("../")
        from dopsearch_py.dopsearch import read_file
        import datetime as dt
        
        self.rad = rad
        self.ftype = ftype
        self.ctr_date = dt.datetime(ctr_date.year, ctr_date.month, ctr_date.day)
        if ffname is None:
            self.ffname = self._construct_filename()
        else:
            self.ffname = ffname

        # collect the data 
        stm = self.ctr_date

        # add two minute to the etm to read the last record from
        # the boxcar filtered concatenated fitacf(ex) files
        etm = self.ctr_date + dt.timedelta(days=1) + dt.timedelta(minutes=2)

        # read data from file. data_from_db has to be False
        self.data = read_file(self.ffname, rad, stm, etm, params,
                              ftype=self.ftype, data_from_db=False,
                              plotrti=False)

    def _construct_filename(self, basedir="../data/"):
        """ constructs filename of a file of interest """

        import datetime as dt
        # expend the time
        stm = self.ctr_date
        etm = self.ctr_date + dt.timedelta(days=1)

        ffname = stm.strftime("%Y%m%d.%H%M%S") + "." + \
                 etm.strftime("%Y%m%d.%H%M%S") + "." + \
                 self.rad + "." + self.ftype + "f"
        ffname = basedir + self.rad + "/" + ffname

        return ffname

    def move_to_db(self, conn):

        import json 
        import sqlite3

        cur = conn.cursor()
        for bmnum in self.data.keys():
            data_dict = self.data[bmnum]
            table_name = self.rad + "_bm" + str(bmnum)
            command = "CREATE TABLE IF NOT EXISTS {tb} (\
                      vel TEXT, rsep REAL, frang REAL, bmazm REAL,\
                      slist TEXT, gsflg TEXT,\
                      datetime TIMESTAMP PRIMARY KEY)".format(tb=table_name)

            cur.execute(command)
            for i in xrange(len(data_dict['datetime'])):
                command = "INSERT OR IGNORE INTO {tb} (vel, rsep, frang, bmazm,\
                            slist, gsflg, datetime) VALUES (?, ?, ?, ?, ?, ?, ?)"\
                            .format(tb=table_name)

                cur.execute(command, (json.dumps(data_dict["vel"][i]), data_dict["rsep"][i],
                            data_dict["frang"][i], data_dict["bmazm"][i],
                            json.dumps(data_dict["slist"][i]),
                            json.dumps(data_dict["gsflg"][i]),\
                            data_dict["datetime"][i]))


        # commit the change
        conn.commit()

def main():
    
    import datetime as dt
    import sqlite3
    from month_to_season import get_season_by_month

    # input parameters
    #sdate = dt.datetime(2008,9,16)
    #edate = dt.datetime(2008,9,18)      # includes the whole edate
    #sdate = dt.datetime(2010,1,14)
    #edate = dt.datetime(2010,1,16)      # includes the whole edate

    #sdate = dt.datetime(2010,12,31)
    #edate = dt.datetime(2011,3,1)
    #sdate = dt.datetime(2011,10,31)
    #edate = dt.datetime(2012,3,1)
    #sdate = dt.datetime(2012,10,31)
    #edate = dt.datetime(2013,1,1)

    #sdate = dt.datetime(2012,4,30)
    #edate = dt.datetime(2012,9,1)

    #sdate = dt.datetime(2011,2,28)
    #edate = dt.datetime(2011,5,1)
    #sdate = dt.datetime(2012,2,29)
    #edate = dt.datetime(2012,5,1)
    #sdate = dt.datetime(2012,8,31)
    #edate = dt.datetime(2012,11,1)

    #season = "winter"
    #season = "summer"
    season = "equinox"

    # set the time interval for each season
    if season == "winter":
        stms = [dt.datetime(2010,12,31), dt.datetime(2011,10,31), dt.datetime(2012,10,31)]
        etms = [dt.datetime(2011,3,1), dt.datetime(2012,3,1), dt.datetime(2013,1,1)]

    if season == "summer":
        stms = [dt.datetime(2011,4,30), dt.datetime(2012,4,30)]
        etms = [dt.datetime(2011,9,1), dt.datetime(2012,9,1)]

    if season == "equinox":
        stms = [dt.datetime(2011,2,28), dt.datetime(2011,8,31), dt.datetime(2012,2,29), dt.datetime(2012,8,31)]
        etms = [dt.datetime(2011,5,1), dt.datetime(2011,11,1), dt.datetime(2012,5,1), dt.datetime(2012,11,1)]

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

    # loop through time interval
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

            # loop through dates:
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

