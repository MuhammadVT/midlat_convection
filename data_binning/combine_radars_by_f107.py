def combine_radars(f107_bin, ftype="fitacf", inbaseLocation="./",
                   outbaseLocation="./", indbName=None, outdbName=None):
    
    """ combines all the ten-min median filtered radar data into a number of master tables based
    on the number of F107 bins. That is, one master table for each F107 bin.
    The results are stored in different db files.
    
    parameters
    ----------
    f107_bin : list
        the lower and upper limit of certain F10.7 bin
    """
    import sqlite3
    import datetime as dt
    import numpy as np 

    # make db connection for dopsearch
    if indbName is None:
        indbName = "ten_min_median_" + ftype + ".sqlite"
    conn_ten = sqlite3.connect(inbaseLocation + indbName)
    cur_ten = conn_ten.cursor()

    # create a db that combines all the radar data into one master table
    if outdbName is None:
        outdbName = "f107_" + str(f107_bin[0]) + "_to_" + str(f107_bin[1]) +\
                    "_" + ftype + ".sqlite"
    conn = sqlite3.connect(outbaseLocation + outdbName)
    cur = conn.cursor()

    # attach ten_min_median db into master db
    cur.execute("ATTACH DATABASE '{db}' as 'rads'".format(db=inbaseLocation+indbName))
    cur.execute("ATTACH DATABASE '{db}' as 'binned_f107'".format(db="../data/sqlite3/gmi_imf/binned_f107.sqlite"))

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
                   SELECT vel, glatc, glonc, gazmc FROM {tb2}\
                   WHERE DATE(datetime) IN (SELECT DATE(datetime) FROM {tb3})\
                   ".format(tb1=T1, tb2="rads."+tbl,\
                   tb3="binned_f107."+"b"+str(f107_bin[0]) + "_b" + str(f107_bin[1]))
        cur.execute(command)

    # commit the change
    conn.commit()

    # detach ten_min_median db from master db
    cur.execute("DETACH DATABASE 'rads'")
    cur.execute("DETACH DATABASE 'binned_f107'")

    # close db connections
    conn_ten.close()
    conn.close()

    return

def master_summary(f107_bin, ftype="fitacf", baseLocation="./", dbName=None):
    
    """ stores the summay of the master table into a different table in the same database.
    """
    import sqlite3
    import datetime as dt
    import numpy as np 

    # make db connection for dopsearch
    if dbName is None:
        dbName = "f107_" + str(f107_bin[0]) + "_to_" + str(f107_bin[1]) + "_" + ftype + ".sqlite"

    # create a db that summarizes the data in the master table
    conn = sqlite3.connect(baseLocation + dbName)
    cur = conn.cursor()
    
    # create a summary table
    T1 = "master"
    T2 = "master_summary"
    cur.execute("CREATE TABLE IF NOT EXISTS {tb}\
                (vel TEXT, median_vel REAL, vel_count INTEGER, glatc REAL, glonc REAL, gazmc INTEGER,\
                 PRIMARY KEY (glatc, glonc, gazmc))".format(tb=T2))

    command = "INSERT OR IGNORE INTO {tb2} (vel, vel_count, glatc, glonc, gazmc)\
               SELECT group_concat(vel), COUNT(vel), glatc, glonc, gazmc FROM {tb1}\
               GROUP BY glatc, glonc, gazmc".format(tb1=T1, tb2=T2)
    cur.execute(command)

    # commit the change
    conn.commit()
    
    # select the velocity data grouping by lat-lon-azm bins
    command = "SELECT rowid, vel FROM {tb2}".format(tb2=T2)
    cur.execute(command)
    rws = cur.fetchall()

    for ii, rw in enumerate(rws):
        rwid, vel_txt = rw
        bin_vel = np.array([float(x) for x in vel_txt.split(",")])

        # take the median value
        median_vel = round(np.median(bin_vel),2)
        
        # populate the table 
        command = "UPDATE {tb} SET median_vel={median_vel}\
                  WHERE rowid=={rwid}".format(tb=T2, median_vel=median_vel, rwid=rwid)
        cur.execute(command)

    # commit the change
    conn.commit()

    # close db connection
    conn.close()

    return

def worker(f107_bins, inbaseLocation, outbaseLocation):

    import datetime as dt

    # input parameters
    ftype = "fitacf"
    #ftype = "fitex"

    for f107_bin in f107_bins: 
        # combine the all the radars' data into one master table 
        combine_radars(f107_bin, ftype=ftype, inbaseLocation=inbaseLocation,
                       outbaseLocation=outbaseLocation, indbName=None, outdbName=None)
        print "created the master table for bin ", f107_bin

        # build a summary table of the master table
        master_summary(f107_bin, ftype=ftype, baseLocation=outbaseLocation, dbName=None)
        print "created the master_summary table"


    return
if __name__ == "__main__":
    import multiprocessing

    #seasons = ["winter", "summer", "equinox"]
    seasons = ["winter"]
    #f107_bins = [[0, 100], [100, 175], [175, 500]]
    #f107_bins = [[0, 105], [105, 125], [125, 500]]
    #f107_bins = [[0, 120], [120, 500]]
    #f107_bins = [[0, 110], [110, 500]]
    #f107_bins = [[140, 500]]
    f107_bins = [[0, 95], [105, 130], [140, 500]]
    jobs = []
    for season in seasons:
        inbaseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_mlt" + "/"
        outbaseLocation=inbaseLocation + "binned_by_f107" + "/"
        p = multiprocessing.Process(target=worker, args=(f107_bins, inbaseLocation, outbaseLocation))
        jobs.append(p)
        p.start()

