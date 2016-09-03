def combine_radars(ftype="fitacf", baseLocation="../data/sqlite3/",
                   dbName=None):
    
    """ combines all the ten-min median filtered radar data into one master table.
    The results are stored in a different db file.
    """
    import sqlite3
    import datetime as dt
    import numpy as np 

    # make db connection for dopsearch
    if dbName is None:
        dbName = "ten_min_median_" + ftype + ".sqlite"
    conn_ten = sqlite3.connect(baseLocation + dbName)
    cur_ten = conn_ten.cursor()

    # create a db that combines all the radar data into one master table
    conn = sqlite3.connect(baseLocation + "master_" + ftype + ".sqlite")
    cur = conn.cursor()

    # attach ten_min_median db into master db
    cur.execute("ATTACH DATABASE '{db}' as 'rads'".format(db=baseLocation+dbName))

    # create a table that combines all the radar data into one master table
    T1 = "master"
    cur.execute("CREATE TABLE IF NOT EXISTS {tb}\
                (vel REAL, glatc REAL, glonc REAL, gazmc INTEGER)"\
                .format(tb=T1))

    # get all the table names
    cur_ten.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
    tbl_names = cur_ten.fetchall()
    tbl_names = [x[0] for x in tbl_names]

    for tbl in tbl_names:
        
        # copy data from ten_min_median db into master db
        command = "INSERT INTO {tb1} (vel, glatc, glonc, gazmc)\
                   SELECT vel, glatc, glonc, gazmc FROM {tb2}"\
                   .format(tb1=T1, tb2="rads."+tbl)
        cur.execute(command)

    # commit the change
    conn.commit()

    # detach ten_min_median db from master db
    cur.execute("DETACH DATABASE 'rads'")

    # close db connections
    conn_ten.close()
    conn.close()

    return

def master_summary(ftype="fitacf", baseLocation="../data/sqlite3/",
                   dbName=None):
    
    """ stores the summay of the master table into a different table in the same database.
    """
    import sqlite3
    import datetime as dt
    import numpy as np 

    # make db connection for dopsearch
    if dbName is None:
        dbName = "master_" + ftype + ".sqlite"
    conn_master = sqlite3.connect(baseLocation + dbName)
    cur_master = conn_master.cursor()

    # create a db that summarizes the data in the master table
    conn = sqlite3.connect(baseLocation + "master_" + ftype + ".sqlite")
    cur = conn.cursor()
    
    # create a summary table
    T1 = "master"
    T2 = "master_summary"
    T3 = "master_copy"
    cur.execute("CREATE TABLE IF NOT EXISTS {tb}\
                (vel TEXT, median_vel REAL, vel_count INTEGER, glatc REAL, glonc REAL, gazmc INTEGER,\
                 PRIMARY KEY (glatc, glonc, gazmc))".format(tb=T2))

    command = "INSERT OR IGNORE INTO {tb2} (vel_count, glatc, glonc, gazmc)\
               SELECT COUNT(vel), glatc, glonc, gazmc FROM {tb1}\
               GROUP BY glatc, glonc, gazmc".format(tb1=T1, tb2=T2)
    cur.execute(command)

    # commit the change
    conn.commit()
    
    # create a copy of the master table
    command = "CREATE TABLE IF NOT EXISTS {tb3}\
               AS SELECT * FROM {tb1}".format(tb1=T1, tb3=T3)
    cur.execute(command)

    # commit the change
    conn.commit()

    # select the velocity data grouping by lat-lon-azm bins
    command = "SELECT rowid, glatc, glonc, gazmc FROM {tb2} ORDER BY vel_count DESC".format(tb2=T2)
    cur.execute(command)
    rws = cur.fetchall()

    for ii, rw in enumerate(rws):
        rwid, lat, lon, az = rw
        command = "SELECT vel, rowid FROM {tb3}\
                   WHERE glatc={lat}\
                    AND glonc={lon}\
                    AND gazmc={az}".format(tb3=T3, lat=lat, lon=lon, az=az)
        cur.execute(command)
        vel_rowid = cur.fetchall()
        bin_vel = [x[0] for x in vel_rowid]
        rowids = [str(x[1]) for x in vel_rowid]

        # delete alread inquired rows
        command = "DELETE FROM {tb} WHERE rowid in ({rowids})"\
                    .format(tb=T3, rowids=",".join(rowids)) 
        cur.execute(command)

        # take the median value
        median_vel = round(np.median(bin_vel),2)
        
        # convert list into commas seperated string
        bin_vel = ",".join([str(x) for x in bin_vel])

        # populate the table 
        command = "UPDATE {tb} SET vel='{bin_vel}', median_vel='{median_vel}'\
                  WHERE rowid=={rwid}".format(tb=T2, bin_vel=bin_vel,
                                              median_vel=median_vel, rwid=rwid)
        cur.execute(command)
        print ii


    # drop master_copy table
    cur.execute("DROP TABLE {tb}".format(tb=T3))

    # commit the change
    conn.commit()

    # close db connection
    conn.close()

    return

def main():

    import datetime as dt

    # input parameters
    ftype = "fitacf"
    #ftype = "fitex"

    season = "winter"
    baseLocation="../data/sqlite3/" + season + "/"
            
    # combine the all the radars' data into one master table 
    combine_radars(ftype=ftype, baseLocation=baseLocation, dbName=None)
    print "created the master table"

    # build a summary table of the master table
    master_summary(ftype=ftype, baseLocation=baseLocation, dbName=None)
    print "created the master_summary table"


    return
if __name__ == "__main__":
    main()
