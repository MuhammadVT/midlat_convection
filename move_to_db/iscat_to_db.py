'''
Created on Aug. 5, 2016

Muhamamd
'''

import pdb

class iscat(object):

    def __init__(self, ctr_date, localdict, 
               params=["velocity"], low_vel_iscat_event_only=True,
               search_allbeams=True, bmnum=7, no_gscat=True, ffname=None):

        """ A class that holds the iscat events that have been obtained by
            file prepareting, file reading, and searching for iscat events.
        
        ctr_date : datetime.datetime
            a full day for which iscat events are to be searched. 
        localdirc : dict
            holds radar related informations 
        params : list
            works for params=["velocity"] only
        search_allbeams : bool
            if set to true, iscat event searching will be performed on all the 
            beams, and ignores the bmnum argument. 
        bmnum : int
            bmnum argument only works in search_allbeams is set to False
        no_gscat : bool
            removes all the gscat
        ffname : string
            if ffname is not set to None, ffname will be constructed.

        Returns : iscat object 
            A dict of dicts in the form of {bmnum:dict}.
            if no_gscat==False, returns data with all its points'
            gsflg marked as 1 (gscat) except for iscat.
            if no_gscat==True, returns only the iscat (gsflg=0)
            
        """

        import sys
        sys.path.append("../")
        from dopsearch_py.dopsearch import iscat_event_searcher
        import datetime as dt
        
        self.rad = localdict["radar"]
        self.ftype = localdict["ftype"]
        self.channel = localdict["channel"]
        self.ctr_date = dt.datetime(ctr_date.year, ctr_date.month, ctr_date.day)
        self.low_vel_iscat_event_only = low_vel_iscat_event_only 
        self.search_allbeams = search_allbeams
        if self.search_allbeams:
            self.bmnum = None
        else:
            self.bmnum = bmnum 
        self.no_gscat = no_gscat
        if ffname is None:
            self.ffname = self._construct_filename()
        else:
            self.ffname = ffname

        # collect the iscat events 
        self.events = iscat_event_searcher(ctr_date, localdict,
                           params=params, low_vel_iscat_event_only=low_vel_iscat_event_only,
                           search_allbeams=search_allbeams, bmnum=bmnum,
                           no_gscat=no_gscat, ffname=self.ffname)

        # remove None type events. If there is no iscat points in an event for a day
        # for a specific beam then the output is {bm:None}
        if self.events is not None:
            self.events = self._remove_None()

    def _remove_None(self):
        for bn in self.events.keys():
            if self.events[bn] is None:
                self.events.pop(bn)
        if self.events == {}:
            self.events = None
        return self.events

    def _construct_filename(self, basedir="../data/"):
        """ constructs filename of a file of interest """

        import datetime as dt
        # expend the time to three days
        stm = self.ctr_date - dt.timedelta(days=1)
        etm = self.ctr_date + dt.timedelta(days=2)

        ffname = stm.strftime("%Y%m%d.%H%M%S") + "." + \
                 etm.strftime("%Y%m%d.%H%M%S") + "." + \
                 self.rad + "." + self.ftype + "f"
        ffname = basedir + self.rad + "/" + ffname

        return ffname

    def join_list_as_str(self, sep=","):
        """ joins list entry to a string with a certain seperator 
        This is useful to enter list entriy as a whole into a db

        
        Returns : nothing
            change each list entry as a joined strings seperated by sep

        """
        
        # join a list entriy as string seperated by sep
        for bn in self.events.keys():

            # find the variable names whose values are to be joined
            kys_tmp = self.events[bn].keys()
            kys = []
            for ky in kys_tmp:
                if isinstance(self.events[bn][ky][0], list):
                    kys.append(ky)

            for ky in kys:
                self.events[bn][ky] = [sep.join([str(round(itm,2)) for itm in x]) for x in self.events[bn][ky]]
        return

    def move_to_db(self, conn, column_map):

        from dbtools.db.dao.DBAccessObject import Table as T, TableDescriptor as TD, ColumnTypes as CT
        from dbtools.db.connection.Connector import Executor

        for bmnum in self.events.keys():
            data_dict = self.events[bmnum]
            table_name = self.rad + "_bm" + str(bmnum)
            td = TD(table_name, column_map)

            e = Executor(conn = conn)
            e._set_command(td._build_create_command())
            e._execute_command()
            for i in xrange(len(data_dict['datetime'])):
                col_value_map = {kywrd:data_dict[kywrd][i] for kywrd in column_map.keys()}
                #col_value_map = {"datetime":data_dict['datetime'][i], "vel":data_dict['vel'][i],
                #                 "slist":data_dict["slist"][i], "bmazm":data_dict["bmazm"][i]} 
                e._set_command(td._build_insert_command(col_value_map))
                e._execute_insert_command(col_value_map.values())
        conn.connection.commit()

def worker(rads, season, baseLocation):
    
    import datetime as dt
    import sys
    sys.path.append("../")
    from dbtools.db.dao.DBAccessObject import Table as T, TableDescriptor as TD, ColumnTypes as CT
    from dbtools.db.connection.Connector import Connectors

    # input parameters
    channel = None
    params=['velocity']
    ftype = "fitacf"
    #ftype = "fitex"
    #low_vel_iscat_event_only=False
    low_vel_iscat_event_only=True
    search_allbeams=True
    no_gscat=True

    # columns names and their types
    # data_keys = ['vel', 'datetime', 'slist', 'rsep', 'nrang', 'frang', 'gsflg', 'bmazm']
    column_map = {"datetime":CT.DATETIME.value, "vel":CT.TEXT.value, "slist":CT.TEXT.value,
                  "rsep":CT.FLOAT.value, "frang":CT.FLOAT.value,
                  "bmazm":CT.FLOAT.value}

    # set the time interval for each season
    if season == "summer":
        stms = [dt.datetime(2011,5,1), dt.datetime(2012,5,1)]
        etms = [dt.datetime(2011,8,31), dt.datetime(2012,8,31)]

    if season == "winter":
        stms = [dt.datetime(2011,1,1), dt.datetime(2011,11,1), dt.datetime(2012,11,1)]
        etms = [dt.datetime(2011,2,28), dt.datetime(2012,2,29), dt.datetime(2012,12,31)]

    if season == "equinox":
        stms = [dt.datetime(2011,3,1), dt.datetime(2011,9,1), dt.datetime(2012,3,1), dt.datetime(2012,9,1)]
        etms = [dt.datetime(2011,4,30), dt.datetime(2011,10,31), dt.datetime(2012,4,30), dt.datetime(2012,10,31)]

    # loop through time interval
    for dd in range(len(stms)):
        sdate = stms[dd]
        edate = etms[dd]
                
        num_days = (edate - sdate).days + 1
        dtm_range = [sdate + dt.timedelta(days=i) for i in xrange(num_days)]
       
        # loop through the radars
        for rad in rads:
            localdict = {"ftype" : ftype, "radar" : rad, "channel" : channel}

            # make a db connection
            dbName = "dopsearch_" + rad + "_" + ftype + ".sqlite"
            conn = Connectors(baseLocation = baseLocation, dbName = dbName,
                              isInMem = False, isAutoSave = False)

            # loop through dates:
            for ctr_date in dtm_range:

                # collect the iscat events 
                t1 = dt.datetime.now()
                print "creating an iscat object from " + rad + " for " + str(ctr_date)
                print "searching all beams of " + rad
                iscat_events = iscat(ctr_date, localdict,
                                     params=params, low_vel_iscat_event_only=low_vel_iscat_event_only,
                                     search_allbeams=search_allbeams, no_gscat=no_gscat, ffname=None)
        

                if iscat_events.events is not None:
                    # join a list entriy as string seperated by sep
                    iscat_events.join_list_as_str(sep=",")

                    # move iscat events to db
                    #t1 = dt.datetime.now()
                    iscat_events.move_to_db(conn, column_map)
                    #t2 = dt.datetime.now()
                    #print ("move_to_db takes " + str((t2-t1).total_seconds() / 60.)) + " mins"
                    print ("iscat has been moved to db")
                else:
                    print "iscat_events.events is None"

                t2 = dt.datetime.now()
                print ("Finishing an iscat object took " + str((t2-t1).total_seconds() / 60.)) + " mins\n"

            # close db connection
            conn._close_connection()

    return

if __name__ == '__main__':
    import multiprocessing

    seasons = ["winter", "summer", "equinox"]
    rads_list = [["bks", "wal", "fhe"], ["fhw", "cve", "cvw"]]
    #seasons = ["equinox"]
    #rads_list = [["fhw", "cve", "cvw"]]
    
    jobs = []
    for season in seasons:
        #baseLocation="../data/sqlite3/" + season + "/dopsearch_data/" + "all_iscat_event" + "/"
        baseLocation="../data/sqlite3/" + season + "/dopsearch_data/" + "low_vel_iscat_event_only" + "/"
        for i in range(len(rads_list)):
            #worker(rads_list[i], season, baseLocation)
            p = multiprocessing.Process(target=worker, args=(rads_list[i], season, baseLocation))
            jobs.append(p)
            p.start()
