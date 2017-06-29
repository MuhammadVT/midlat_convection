class gmi_imf_to_db(object):
    """ collects geomagnetic indices into db """


    def __init__(self, stm, etm, dbName=None, baseLocation="../data/sqlite3/"):

        self.stm = stm
        self.etm = etm
        if dbName is None:
            self.dbName = "gmi_imf.sqlite"
        else:
            self.dbName = dbName
        self.baseLocation = baseLocation

    def _create_dbconn(self):
        import sqlite3

        # make a db connection
        conn = sqlite3.connect(self.baseLocation + self.dbName,
                               detect_types = sqlite3.PARSE_DECLTYPES)
        return conn

    def imf_to_db(self):

        import datetime as dt
        from davitpy.gme.ind import readOmni
        import numpy as np

        data_dict = {'datetime':[], 'Bx':[], 'By':[], 'Bz':[]}
        # read the data we want in GSM coords
        omni_list = readOmni(sTime=self.stm, eTime=self.etm, res=1)
        data_dict['datetime'] = [omni_list[i].time for i in range(len(omni_list))]
        data_dict['Bx'] = [omni_list[i].bx for i in range(len(omni_list))]
        data_dict['By'] = [omni_list[i].bym for i in range(len(omni_list))]
        data_dict['Bz'] = [omni_list[i].bzm for i in range(len(omni_list))]
        # clock angle
        #data_dict['theta'] = np.degrees(np.arctan2(data_dict['By'], data_dict['Bz'])) % 360
        data_dict['theta'] = [round(np.degrees(np.arctan2(data_dict['By'][i], data_dict['Bz'][i])) % 360, 2) \
                              if (data_dict['By'][i] is not None and data_dict['Bz'][i] is not None) \
                              else None \
                              for i in range(len(data_dict['By']))]

        
        # make db connection
        self.conn = self._create_dbconn()
        # create table
        table_name = "IMF"
        colname_type = "datetime TIMESTAMP PRIMARY KEY, Bx REAL, By REAL, Bz REAL, theta REAL"
        command = "CREATE TABLE IF NOT EXISTS {tb} ({colname_type})".format(tb=table_name,\
                                                                           colname_type=colname_type)
        self.conn.cursor().execute(command)
        
        # populate the table
        row_num = len(data_dict['datetime'])
        columns = "datetime, Bx, By, Bz, theta"
        for i in xrange(row_num):
            dtm = data_dict['datetime'][i]
            Bx = data_dict['Bx'][i]
            By = data_dict['By'][i]
            Bz = data_dict['Bz'][i]
            theta = data_dict['theta'][i]
            if (Bx is not None) and (By is not None) and (Bz is not None):
                command = "INSERT OR IGNORE INTO {tb}({columns}) VALUES (?, ?, ?, ?, ?)".\
                          format(tb=table_name, columns=columns)
                self.conn.cursor().execute(command, (dtm, Bx, By, Bz, theta))
        self.conn.commit()

        # close db connection
        self.conn.close()

    def f107_to_db(self):

        import datetime as dt
        import pandas as pd

        # read the data from txt file 
        date_parser = lambda x: pd.datetime.strptime(x, '%Y %j %H')
        fname = "../data/sqlite3/gmi_imf/F107.txt"
        df = pd.read_csv(fname, index_col=0, header=None, names=["Year", "DOY", "Hour", "F107"],
                         skipinitialspace=True, delim_whitespace=True,
                         parse_dates={'datetime': [0,1,2]}, date_parser=date_parser, na_values=999.9)

        data_dict = {}
        data_dict['datetime'] = df.index.to_pydatetime().tolist() 
        data_dict['F107'] = df.F107.tolist() 
        
        # make db connection
        self.conn = self._create_dbconn()
        # create table
        table_name = "F107"
        colname_type = "datetime TIMESTAMP PRIMARY KEY, F107 REAL"
        command = "CREATE TABLE IF NOT EXISTS {tb} ({colname_type})".format(tb=table_name,\
                                                                           colname_type=colname_type)
        self.conn.cursor().execute(command)
        
        # populate the table
        row_num = len(data_dict['datetime'])
        columns = "datetime, F107"
        for i in xrange(row_num):
            dtm = data_dict['datetime'][i]
            F107 = round(data_dict['F107'][i],1)
            command = "INSERT OR IGNORE INTO {tb}({columns}) VALUES (?, ?)".\
                      format(tb=table_name, columns=columns)
            self.conn.cursor().execute(command, (dtm, F107))
        self.conn.commit()

        # close db connection
        self.conn.close()


    def kp_to_db(self, kp_lim=None):

        from davitpy import gme
        import datetime as dt

        self.conn = self._create_dbconn()
        # Kp data
        data_dict = {'datetime':[], 'kp':[]}
        Kp_list = gme.ind.readKp(sTime=self.stm,eTime=self.etm)
        day_num = (self.etm-self.stm).days
        for n in xrange(day_num):
            kp_tmp = Kp_list[n].kp
            time_tmp = Kp_list[n].time
            if kp_tmp is not None:
                for l in range(len(kp_tmp)):
                    if len(kp_tmp) < 8:
                        print str(len(kp_tmp)) + " values for a day, should have 8 values"
                    if len(kp_tmp[l])== 2:
                        if kp_tmp[l][1] == '+':
                            data_dict['kp'].append(int(kp_tmp[l][0])+0.3)
                        elif kp_tmp[l][1] == '-':
                            data_dict['kp'].append(int(kp_tmp[l][0])-0.3)
                    else:
                        data_dict['kp'].append(int(kp_tmp[l][0]))
                    data_dict['datetime'].append(time_tmp + dt.timedelta(hours=3*l))

        # move to db
        if (data_dict['datetime'])!=[]:

            # create table
            table_name = "kp"
            colname_type = "datetime TIMESTAMP PRIMARY KEY, kp REAL"
            command = "CREATE TABLE IF NOT EXISTS {tb} ({colname_type})".format(tb=table_name,\
                                                                           colname_type=colname_type)
            self.conn.cursor().execute(command)
            
            # populate the table
            row_num = len(data_dict['datetime'])
            columns = "datetime, kp"
            for i in xrange(row_num):
                dtm = data_dict['datetime'][i]
                k = data_dict['kp'][i]
                if kp_lim is not None:
                    if k < kp_lim[0] or k >= kp_lim[1]:
                        store_kp = False
                    else:
                        store_kp = True
                else:
                    store_kp = True
                if store_kp:
                    command = "INSERT OR IGNORE INTO {tb}({columns}) VALUES (?, ?)".\
                              format(tb=table_name, columns=columns)
                    self.conn.cursor().execute(command, (dtm, k))

                    
                #col_value_map = {"datetime":data_dict['datetime'][i], "vel":data_dict['vel'][i],
                #                 "slist":data_dict["slist"][i], "bmazm":data_dict["bmazm"][i]} 
            self.conn.commit()

        # close db connection
        self.conn.close()



    def symh_to_db(self, symh_lim=None):

        from davitpy import gme
        import datetime as dt

        # make a db connection
        self.conn = self._create_dbconn()

        # read SYMH data
        data_dict = {'datetime':[], 'symh':[]}
        sym_list = gme.ind.symasy.readSymAsy(sTime=self.stm,eTime=self.etm)
        for i in xrange(len(sym_list)):
            data_dict['symh'].append(sym_list[i].symh)
            data_dict['datetime'].append(sym_list[i].time)

        # move to db
        if (data_dict['datetime'])!=[]:

            # create table
            table_name = "symh"
            colname_type = "datetime TIMESTAMP PRIMARY KEY, symh REAL"
            command = "CREATE TABLE IF NOT EXISTS {tb} ({colname_type})".format(tb=table_name,\
                                                                           colname_type=colname_type)
            self.conn.cursor().execute(command)
            
            # populate the table
            row_num = len(data_dict['datetime'])
            columns = "datetime, symh"
            for i in xrange(row_num):
                dtm = data_dict['datetime'][i]
                k = data_dict['symh'][i]
                if symh_lim is not None:
                    if k < symh_lim[0] or k >= symh_lim[1]:
                        store_symh = False
                    else:
                        store_symh = True
                else:
                    store_symh = True
                if store_symh:
                    command = "INSERT OR IGNORE INTO {tb}({columns}) VALUES (?, ?)".\
                              format(tb=table_name, columns=columns)
                    self.conn.cursor().execute(command, (dtm, k))

                    
                #col_value_map = {"datetime":data_dict['datetime'][i], "vel":data_dict['vel'][i],
                #                 "slist":data_dict["slist"][i], "bmazm":data_dict["bmazm"][i]} 
            self.conn.commit()

        # close db connection
        self.conn.close()

    def ae_to_db(self, ae_lim=None):

        from davitpy import gme
        import datetime as dt

        # make a db connection
        self.conn = self._create_dbconn()

        # read AE data
        data_dict = {'datetime':[], 'ae':[]}
        # reads AE data with 1-min resolution
        ae_list = gme.ind.ae.readAe(sTime=self.stm, eTime=self.etm, res=1)
        for i in xrange(len(ae_list)):
            data_dict['ae'].append(ae_list[i].ae)
            data_dict['datetime'].append(ae_list[i].time)

        # move to db
        if (data_dict['datetime'])!=[]:

            # create table
            table_name = "ae"
            colname_type = "datetime TIMESTAMP PRIMARY KEY, ae REAL"
            command = "CREATE TABLE IF NOT EXISTS {tb} ({colname_type})".format(tb=table_name,\
                                                                           colname_type=colname_type)
            self.conn.cursor().execute(command)
            
            # populate the table
            row_num = len(data_dict['datetime'])
            columns = "datetime, ae"
            for i in xrange(row_num):
                dtm = data_dict['datetime'][i]
                k = data_dict['ae'][i]
                if ae_lim is not None:
                    if k < ae_lim[0] or k >= ae_lim[1]:
                        store_ae = False
                    else:
                        store_ae = True
                else:
                    store_ae = True
                if store_ae:
                    command = "INSERT OR IGNORE INTO {tb}({columns}) VALUES (?, ?)".\
                              format(tb=table_name, columns=columns)
                    self.conn.cursor().execute(command, (dtm, k))

            self.conn.commit()

        # close db connection
        self.conn.close()



def main():
    import datetime as dt
    stm = dt.datetime(2010, 12, 31)
    etm = dt.datetime(2013, 1, 2)
    #etm = dt.datetime(2011, 1, 2)
    dbName = None
    baseLocation = "../data/sqlite3/gmi_imf/"

    kp_lim = None
    symh_lim = None
    dst_lim = None
    ae_lim = None

    # create an object
    gmi = gmi_imf_to_db(stm, etm, dbName=dbName, baseLocation=baseLocation)

#    # store IMF into db
#    print "storing IMF to db"
#    gmi.imf_to_db()
#    print "imf is done"

#    # store F107 into db
#    print "storing F107 to db"
#    gmi.f107_to_db()
#    print "F107 is done"


#    # store kp into db
#    print "storing kp to db"
#    gmi.kp_to_db(kp_lim=kp_lim)
#    print "kp is done"

#    # store symh into db
#    print "storing symh to db"
#    gmi.symh_to_db(symh_lim=symh_lim)
#    print "symh is done"

    # store symh into db
    print "storing AE to db"
    gmi.ae_to_db(ae_lim=ae_lim)
    print "AE is done"

if __name__ == "__main__":
    main()
