'''
Created on Aug. 5, 2016

Muhamamd
'''

import pdb

class iscatold(object):

    def __init__(self, stm, etm, localdict, 
                 params=["velocity"], search_allbeams=True,
                 bmnum=7, no_gscat=True, ffname=None):

        """ A class that holds the iscat points that have been classified using 
        the traditional method.
        
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
        ffname : string
            if ffname is not set to None, ffname will be constructed.

        Returns : iscat object 
            A dict of dicts in the form of {bmnum:dict}.
            
        """

        import datetime as dt
        
        self.rad = localdict["radar"]
        self.ftype = localdict["ftype"]
        self.channel = localdict["channel"]
        self.params = params
        self.stm = dt.datetime(stm.year, stm.month, stm.day)
        self.etm = dt.datetime(etm.year, etm.month, etm.day)
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
        self.events = self._create_events()

        # remove None type events. If there is no iscat points in an event for a day
        # for a specific beam then the output is {bm:None}
        self.events = self._remove_None()

    def _construct_filename(self, basedir="../data/"):
        """ constructs filename of a file of interest """

        import datetime as dt
        # expend the time to three days

        ffname = self.stm.strftime("%Y%m%d.%H%M%S") + "." + \
                 self.etm.strftime("%Y%m%d.%H%M%S") + "." + \
                 self.rad + "." + self.ftype + "f"
        ffname = basedir + self.rad + "/" + ffname

        return ffname


    def _create_events(self):
        import sys
        sys.path.append("../")
        from dopsearch_py.dopsearch import read_file

        allbeams = read_file(self.ffname, self.rad, self.stm, self.etm,self.params,
                            ftype=self.ftype)

        kys = ['datetime', "vel", "slist"]
        data_dict = {}
        # select only the keywords of interest 
        for bn in allbeams.keys():
            if allbeams[bn]["vel"] != []:
                data_dict[bn] = {ky:[] for ky in kys}
                data_dict[bn]['datetime'] = allbeams[bn]['datetime']
                data_dict[bn]["vel"] = [[] for tm in allbeams[bn]['datetime']]
                data_dict[bn]["slist"] = [[] for tm in allbeams[bn]["slist"]]
            else:
                allbeams.pop(bn)
        for bn in allbeams.keys():
            for i, lst in enumerate(allbeams[bn]["gsflg"]):
                if lst is None:
                    data_dict[bn]["vel"][i] = None
                    data_dict[bn]["slist"][i] = None
                else:
                    for j, itm in enumerate(lst):
                        if (allbeams[bn]["gsflg"][i][j]==0) and (allbeams[bn]["slist"][i][j]>=7):
                            data_dict[bn]["vel"][i].append(allbeams[bn]["vel"][i][j])
                            data_dict[bn]["slist"][i].append(allbeams[bn]["slist"][i][j])
                    if data_dict[bn]["vel"][i]==[]:
                        data_dict[bn]["vel"][i] = None
                        data_dict[bn]["slist"][i] = None

        return data_dict

        
    def _remove_None(self):
        for bn in self.events.keys():
            num_iter = len(self.events[bn]['datetime'])
            self.events[bn]['datetime'] = [x for i,x in enumerate(self.events[bn]['datetime']) if self.events[bn]["vel"][i] is not None]
            self.events[bn]["vel"] = [x for x in self.events[bn]["vel"] if x is not None]
            self.events[bn]["slist"] = [x for x in self.events[bn]["slist"] if x is not None]
            if self.events[bn]['datetime'] == []:
                self.events.pop(bn)
        if self.events == {}:
            self.events = None

        return self.events


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
                col_value_map = {"datetime":data_dict['datetime'][i], "vel":data_dict['vel'][i],
                                 "slist":data_dict['slist'][i]}
                e._set_command(td._build_insert_command(col_value_map))
                e._execute_insert_command(col_value_map.values())
        conn.connection.commit()

def main():
    
    import datetime as dt
    import sys
    sys.path.append("../")
    from dbtools.db.dao.DBAccessObject import Table as T, TableDescriptor as TD, ColumnTypes as CT
    from dbtools.db.connection.Connector import Connectors

    # input parameters
    #stm = dt.datetime(2008,10,1)
    #etm = dt.datetime(2008,12,31)
    stm = dt.datetime(2011,10,1)
    etm = dt.datetime(2011,12,31)
    num_days = (etm - stm).days
    dtm_range = [stm + dt.timedelta(days=i) for i in range(0, num_days, 3)]

    #rad_list = ["bks", "wal", "fhe", "fhw", "cve", "cvw"]
    #rad_list = ["bks"]
    rad_list = ["wal"]
    channel = None
    params=['velocity']
    ftype = "fitacf"
    #ftype = "fitex"
    #ffname = "../data/20100114.000000.20100117.000000.bks.fitacff"

    # columns names and their types
    # data_keys = ['vel', 'datetime', 'slist', 'rsep', 'nrang', 'frang', 'gsflg', 'bmazm']
    column_map = {"datetime":CT.DATETIME.value, "vel":CT.TEXT.value, "slist":CT.TEXT.value}

    # make a db connection
    baseLocation = "../data/sqlite3/"
    
    # loop through the radars
    for rad in rad_list:
        localdict = {"ftype" : ftype, "radar" : rad, "channel" : channel}

        # make a db connection
        dbName = "iscatold_" + rad + "_" + ftype + ".sqlite"
        conn = Connectors(baseLocation = baseLocation, dbName = dbName,
                          isInMem = False, isAutoSave = False)


        # loop through dates:
        for stm in dtm_range:

            # collect the iscatold events 
            t1 = dt.datetime.now()
            print "creating an iscatold object for " + str(stm) + "\n"

            iscatold_events = iscatold(stm, stm+dt.timedelta(days=3), localdict, params=params,
                                    search_allbeams=True, no_gscat=True, ffname=None)
    
            t2 = dt.datetime.now()
            print ("creating an iscatold object takes " + str((t2-t1).total_seconds() / 60.)) + " mins"

            if iscatold_events is not None:
                # join a list entriy as string seperated by sep
                iscatold_events.join_list_as_str(sep=",")

                # move iscatold events to db
                t1 = dt.datetime.now()
                iscatold_events.move_to_db(conn, column_map)
                t2 = dt.datetime.now()
                print ("move_to_db takes " + str((t2-t1).total_seconds() / 60.)) + " mins"

    # close db connection
    conn._close_connection()

    return


if __name__ == '__main__':
    #pass
    main()

