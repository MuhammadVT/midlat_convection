def combine_radars(imf_bin, ftype="fitacf", inbaseLocation="./",
                   outbaseLocation="./", indbName=None, outdbName=None):
    
    """ combines all the ten-min median filtered radar data into a number of master tables based
    on the number of imf clock angle bins. That is, one master table for each imf clock angle bin.
    The results are stored in different db files.
    
    parameters
    ----------
    imf_bin : list
        the lower and upper limit of certain imf clock angle bin
    """
    import sqlite3
    import datetime as dt
    import numpy as np 

    # make db connection for dopsearch
    if indbName is None:
        indbName = "ten_min_median_" + ftype + ".sqlite"
    conn_ten = sqlite3.connect(inbaseLocation + indbName, detect_types=sqlite3.PARSE_DECLTYPES)
    cur_ten = conn_ten.cursor()

    # create a db that combines all the radar data into one master table
    if outdbName is None:
        outdbName = "imf_clock_angle_" + str(imf_bin[0]) + "_to_" + str(imf_bin[1]) +\
                    "_" + ftype + ".sqlite"
    conn = sqlite3.connect(outbaseLocation + outdbName, detect_types=sqlite3.PARSE_DECLTYPES)
    cur = conn.cursor()

    # attach ten_min_median db into master db
    cur.execute("ATTACH DATABASE '{db}' as 'rads'".format(db=inbaseLocation+indbName))
    cur.execute("ATTACH DATABASE '{db}' as 'binned_imf'".format(db="../data/sqlite3/gmi_imf/binned_imf.sqlite"))

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
                   WHERE datetime IN (SELECT datetime FROM {tb3})\
                   ".format(tb1=T1, tb2="rads."+tbl,\
                   tb3="binned_imf."+"b"+str(imf_bin[0]) + "_b" + str(imf_bin[1]))
        cur.execute(command)

    # commit the change
    conn.commit()

    # detach ten_min_median db from master db
    cur.execute("DETACH DATABASE 'rads'")
    cur.execute("DETACH DATABASE 'binned_imf'")

    # close db connections
    conn_ten.close()
    conn.close()

    return

def master_summary(imf_bin, ftype="fitacf", baseLocation="./", dbName=None):
    
    """ stores the summay of the master table into a different table in the same database.
    """
    import sqlite3
    import datetime as dt
    import numpy as np 

    # make db connection for dopsearch
    if dbName is None:
        dbName = "imf_clock_angle_" + str(imf_bin[0]) + "_to_" + str(imf_bin[1]) + "_" + ftype + ".sqlite"

    # create a db that summarizes the data in the master table
    conn = sqlite3.connect(baseLocation + dbName, detect_types=sqlite3.PARSE_DECLTYPES)
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

def worker(imf_bins, inbaseLocation, outbaseLocation):

    import datetime as dt

    # input parameters
    ftype = "fitacf"
    #ftype = "fitex"

    for imf_bin in imf_bins: 
        # combine the all the radars' data into one master table 
        combine_radars(imf_bin, ftype=ftype, inbaseLocation=inbaseLocation,
                       outbaseLocation=outbaseLocation, indbName=None, outdbName=None)
        print "created the master table for bin ", imf_bin

        # build a summary table of the master table
        master_summary(imf_bin, ftype=ftype, baseLocation=outbaseLocation, dbName=None)
        print "created the master_summary table"


    return
if __name__ == "__main__":
    import multiprocessing

    seasons = ["winter", "summer", "equinox"]
    #seasons = ["winter"]
    #imf_bins = [[65, 115], [245, 295]]
    #imf_bins = [[335, 25], [155, 205]]
    #imf_bins = [[330, 30], [150, 210]]
    #imf_bins = [[60, 120], [240, 300]]

    #imf_bins = [[315, 45], [135, 225]]
    imf_bins = [[300, 60], [120, 240]]


    jobs = []
    for season in seasons:
        inbaseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_mlt" + "/"
        outbaseLocation=inbaseLocation + "binned_by_imf_clock_angle" + "/"
        p = multiprocessing.Process(target=worker, args=(imf_bins, inbaseLocation, outbaseLocation))
        jobs.append(p)
        p.start()
    
    for p in jobs:
        p.join()

