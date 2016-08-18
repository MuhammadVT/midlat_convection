#!/usr/bin/python

def gmi_based_filter(rad, ftype="fitacf", isKp_based=True, kp_lim=[0,3],
                     isSymH_based=False, symh_min=-50,
                     dopsearch_dbName=None, gmi_dbName=None,
                     baseLocation="./"):
    """ Selects the quiet time data based on gmag indicies.

    Note
    ----
        SymH base filter has not been implemented yet
    
    """

    import sys
    sys.path.append("../../")
    from dbtools.db.connection.Connector import Connectors, Executor
    import datetime as dt
    import sqlite3

    # make db connection for dopsearch
    if dopsearch_dbName is None:
        dopsearch_dbName = "dopsearch_" + rad + "_" + ftype + ".sqlite"
    conn_dop = Connectors(baseLocation = baseLocation, dbName = dopsearch_dbName,
                      isInMem = False, isAutoSave = False)
    executor_dop = Executor(conn = conn_dop)

    # make db connection for gmi
    if gmi_dbName is None:
        gmi_dbName = "gmi.sqlite"
    conn_gmi = Connectors(baseLocation = baseLocation, dbName = gmi_dbName,
                      isInMem = False, isAutoSave = False)
    executor_gmi = Executor(conn = conn_gmi)

#    if isSymH_based:
#        tbl_sym = "symh"
#        executor_gmi._set_command("SELECT datetime FROM {tb} WHERE symh >= {sym_min}"\
#                                   .format(tb=tbl_sym, symh_min=symh_min))
#        executor_gmi._execute_command()
#        dtms = executor_gmi.conn.cursor.fetchall()


    # get all the tables names in dopsearch db
    executor_dop._set_command("SELECT name FROM sqlite_master WHERE type = 'table'")
    executor_dop._execute_command()
    tbl_names = executor_dop.conn.cursor.fetchall()
    tbl_names = [x[0] for x in tbl_names]

    # select the time intervals base on kp_lim
    if isKp_based:
        tbl_nm = "kp"
        command = "SELECT datetime FROM {tb} WHERE kp >= {kp_min} AND\
                  kp < {kp_max}".format(tb=tbl_nm, kp_min=kp_lim[0],kp_max=kp_lim[1])

        executor_gmi._set_command(command)
        executor_gmi._execute_command()
        dtms = executor_gmi.conn.cursor.fetchall()
        dtms = [x[0] for x in dtms]

    # filter dopsearch db based on kp
    # do it one table at a time
    for tbl in tbl_names:

        # store rowids_dop in a memory db
        tmp_con = sqlite3.connect("file::memory:?cache=shared")
        with tmp_con:
            
            # attache the in-memory database into dopsearch db
            conn_dop.cursor.execute("ATTACH DATABASE 'file::memory:?cache=shared' AS inmem")

            tmp_cur = tmp_con.cursor()    
            # drop of the tmp table if it exists
            tmp_cur.execute("DROP TABLE IF EXISTS tmp")

            # create a tmp table inside the memory db
            tmp_cur.execute("CREATE TABLE tmp (rid INTEGER PRIMARY KEY)")

            # loop through every 3 hour interval of kp
            for dtm in dtms:
                sdtm = dtm
                edtm = sdtm + dt.timedelta(hours=3)

                # select row ids of needed records
                command = "SELECT rowid FROM {tb} WHERE DATETIME(datetime) >= ?\
                          AND DATETIME(datetime) < ?".format(tb=tbl)
                conn_dop.cursor.execute(command, (sdtm, edtm))
                rowids_tmp = executor_dop.conn.cursor.fetchall()
                
                if len(rowids_tmp) > 0:
                    rowids_tmp = [x[0] for x in rowids_tmp]
                    for rwid in rowids_tmp:

                        # populate the tmp table in the memory db
                        tmp_cur.execute("INSERT OR IGNORE INTO tmp (rid) VALUES (?)", (rwid,))
            
            ## commit the change
            #tmp_con.connection.commit()

            # delete unwanted records from dopsearch db
            command = "DELETE FROM {tb} WHERE rowid NOT IN (SELECT rid FROM inmem.tmp)"\
                      .format(tb=tbl)
            conn_dop.cursor.execute(command)

            # commit the change
            conn_dop.connection.commit()
            print "finished filtering " + tbl + " of " + rad
            
            # detach teh im-memory db from dopsearch db
            conn_dop.cursor.execute("DETACH DATABASE 'inmem'")
            
    # close the connection
    conn_dop.connection.close()

def main():


    # input parameters
    #rad_list = ["bks", "wal", "fhe", "fhw", "cve", "cvw"]
    rad_list = ["bks"]
    #rad_list = ["wal"]
    ftype = "fitacf"
    dopsearch_dbName= None
    gmi_dbName = None
    baseLocation="./"

    for rad in rad_list:
        gmi_based_filter(rad, ftype=ftype, isKp_based=True, kp_lim=[0,3],
                         dopsearch_dbName=dopsearch_dbName, gmi_dbName=gmi_dbName,
                         baseLocation=baseLocation)
    return

if __name__ == "__main__":
    main()

