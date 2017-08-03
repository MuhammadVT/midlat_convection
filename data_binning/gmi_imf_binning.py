class gmi_imf_binning(object):
    """ bins geomagnetic indices and IMF data into various bins"""


    def __init__(self, stm, etm, indbName=None, baseLocation="../data/sqlite3/gmi_imf/"):

        self.stm = stm
        self.etm = etm
        if indbName is None:
            self.indbName = "gmi_imf.sqlite"
        else:
            self.indbName = indbName
        self.baseLocation = baseLocation

    def _create_dbconn(self):
        import sqlite3

        # make a db connection
        conn = sqlite3.connect(self.baseLocation + self.indbName,
                               detect_types = sqlite3.PARSE_DECLTYPES)
        return conn

    def bin_imf(self, outdbName, bin_by_clock_angle=True):

        import sqlite3
        import datetime as dt
        import numpy as np

        # make a db connection
        conn = sqlite3.connect(self.baseLocation + outdbName,
                               detect_types = sqlite3.PARSE_DECLTYPES)
        cur = conn.cursor()
        
        # attach the sorce db
        conn.execute("ATTACH DATABASE '{db}' as 'source'".format(db=self.baseLocation+self.indbName)) 

        # create table
        if bin_by_clock_angle:
            #bins = [[65, 115], [245, 295]]
            #bins = [[60, 120], [240, 300]]
            #bins = [[-25, 25], [155, 205]]
            #bins = [[-30, 30], [150, 210]]

            bins = [[-60, 60], [120, 240]]

            before_mins = 20
            after_mins = 10
            bvec_max = 0.9  # max lengh of the  bias vector
            del_tm = 10
            num_lim = before_mins + del_tm + after_mins - 5    # minimum number of points  

            # create tables: one table for each bin
            for bn in bins:
                # create a table to store each del_tm quite intervals of IMF
                table_name = "b" + str((bn[0]%360)) + "_b" + str(bn[1]%360) 
                colname_type = "datetime TIMESTAMP PRIMARY KEY, avg_clock_angle REAL"
                #colname_type = "datetime TIMESTAMP PRIMARY KEY, Bx REAL, By REAL, Bz REAL, avg_clock_angle REAL"
                command = "CREATE TABLE IF NOT EXISTS {tb} ({colname_type})".format(tb=table_name,\
                                                                                   colname_type=colname_type)
                cur.execute(command)

                # create a table to store each minute of IMF that belongs to the quite intervals
                table_name_2 = table_name + "_all"
                colname_type = "datetime TIMESTAMP PRIMARY KEY, Bx REAL, By REAL, Bz REAL, clock_angle REAL"
                command = "CREATE TABLE IF NOT EXISTS {tb} ({colname_type})".format(tb=table_name_2,\
                                                                                   colname_type=colname_type)
                cur.execute(command)

            # retreive the data
            iter_num = int((self.etm - self.stm.replace(minute=5)).total_seconds() / 60. / del_tm) + 1
            for i in xrange(iter_num):
                tm_target = self.stm + dt.timedelta(minutes=del_tm*i+5)
                tm_before = tm_target - dt.timedelta(minutes=(before_mins+5))
                tm_after = tm_target + dt.timedelta(minutes=(after_mins+5))
                command = "SELECT * FROM {tb2} WHERE (datetime >= ? AND datetime < ?)".\
                          format(tb2="source.IMF")
                cur.execute(command, (tm_before, tm_after))
                rows = cur.fetchall()
                if len(rows) < num_lim:
                    continue
                else:
                    # normolize the (By, Bz) vectors
                    bt = [np.array([x[2], x[3]]) / np.linalg.norm(np.array([x[2], x[3]])) for x in rows]
                    # find the bias vector
                    bias_vec = np.array([np.mean([x[0] for x in bt]), np.mean([x[1] for x in bt])])
                    #print "bias_vec = ", round(np.linalg.norm(bias_vec),2)
                    bias_vec_norm = np.linalg.norm(bias_vec)
                    if bias_vec_norm < bvec_max:
                        continue
                    clk_angle = np.degrees(np.arctan2(bias_vec[0], bias_vec[1]))
                    for bn in bins:
                        if np.sign(bn[0]) < 0 or np.sign(bn[1]) < 0:
                            clk_angle = round(clk_angle, 2)
                        else:
                            clk_angle = round(clk_angle % 360, 2)
                        table_name = "b" + str((bn[0]%360)) + "_b" + str(bn[1]%360) 
                        table_name_2 = table_name + "_all"
                        if clk_angle>=bn[0] and clk_angle<=bn[1]: 
                            clk_angle = round(clk_angle % 360, 2)
                            # populate the table that stores the quite time intervals
                            command = "INSERT OR IGNORE INTO {tb} (datetime, avg_clock_angle) VALUES (?, ?)"\
                            .format(tb=table_name)
                            cur.execute(command, (tm_target, clk_angle))

                            # populate the table that stores quite time IMF for each minute
                            for rw in rows:
                                command = "INSERT OR IGNORE INTO {tb} (datetime, Bx, By, Bz, clock_angle) VALUES (?, ?, ?, ?, ?)"\
                                .format(tb=table_name_2)
                                cur.execute(command, (rw[0], rw[1], rw[2], rw[3], rw[4]))
                        else:
                            continue
        conn.commit()
        
        # detach the source db
        cur.execute("DETACH DATABASE 'source'")

        # close db connection
        conn.close()

    def bin_f107(self, outdbName, bins=[[0,100], [100, 150], [150, 500]]):

        import datetime as dt
        import sqlite3

        # make a db connection
#        self.conn = self._create_dbconn()
        conn = sqlite3.connect(self.baseLocation + outdbName,
                               detect_types = sqlite3.PARSE_DECLTYPES)
        
        # attach the sorce db
        conn.execute("ATTACH DATABASE '{db}' as 'source'".format(db=self.baseLocation+self.indbName)) 

        # create table
        for bn in bins:
            table_name = "b" + str(bn[0]) + "_b" + str(bn[1]) 
            colname_type = "datetime TIMESTAMP PRIMARY KEY, F107 REAL"
            command = "CREATE TABLE IF NOT EXISTS {tb} ({colname_type})".format(tb=table_name,\
                                                                               colname_type=colname_type)
            conn.cursor().execute(command)
            
            # populate the table
            command = "INSERT OR IGNORE INTO {tb} SELECT * FROM {tb2} WHERE\
                      (datetime >= ? AND datetime < ?) AND (F107 >= {bn_low} AND F107 < {bn_high})".\
                      format(tb=table_name, tb2="source.F107", bn_low=bn[0], bn_high=bn[1])
            conn.cursor().execute(command, (self.stm, self.etm))
            conn.commit()
        
        # detach the source db
        conn.cursor().execute("DETACH DATABASE 'source'")

        # close db connection
        #self.conn.close()
        conn.close()

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


def main():
    import datetime as dt
    stm = dt.datetime(2011, 1, 1)
    etm = dt.datetime(2013, 1, 1)
    indbName = None
    baseLocation = "../data/sqlite3/gmi_imf/"

    kp_lim = None
    symh_lim = None
    dst_lim = None

    # create an object
    gmi = gmi_imf_binning(stm, etm, indbName=indbName, baseLocation=baseLocation)

#    # bin F107 
#    print "binning F107 to db"
#    #bins=[[0,105], [105, 130], [130, 500]]
#    #bins=[[0,105], [105, 125], [125, 500]]
#    #bins=[[0,110], [110, 500]]
#    #bins=[[0,120], [120, 500]]
#    bins=[[140,500]]
#    gmi.bin_f107("binned_f107.sqlite", bins=bins)
#    print "done"

    # bin IMF clock angle 
    print "binning IMF"
    gmi.bin_imf("binned_imf.sqlite", bin_by_clock_angle=True)
    print "done"


if __name__ == "__main__":
    main()
