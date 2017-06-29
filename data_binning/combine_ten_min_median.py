def combine_ten_min_median(rad, ftype="fitacf", baseLocation="../data/sqlite3/",
                                dbName=None):
    import sqlite3
    import numpy as np 

    # make db connection
    if dbName is None:
        dbName = "ten_min_median_" + rad + "_" + ftype + ".sqlite"
#    conn_rad = sqlite3.connect(baseLocation + dbName)
#    cur_rad = conn_rad.cursor()

    # create a db that stores each radar data as a table
    conn = sqlite3.connect(baseLocation + "ten_min_median_" + ftype + ".sqlite")
    cur = conn.cursor()

    # create a table that combines all the beams of a radar
    cur.execute("CREATE TABLE IF NOT EXISTS {tb}\
                (vel REAL, glatc REAL, glonc REAL, gazmc REAL, datetime TIMESTAMP)"\
                .format(tb=rad))

    # attach ten_min_median_rad db into ten_min_median db
    cur.execute("ATTACH DATABASE '{db}' as 'rad'".format(db=baseLocation+dbName))

    # move db file 
    cur.execute("INSERT INTO {tb1} SELECT * FROM {tb2}".format(tb1=rad, tb2="rad."+rad))

    # detach the tmpdb from dopsearch db
    cur.execute("DETACH DATABASE 'rad'")

    conn.close()


def worker(baseLocation):
    # input parameters
    rad_list = ["bks", "wal", "fhe", "fhw", "cve", "cvw"]
    ftype = "fitacf"
            
    for rad in rad_list:

        # take ten minutes median values
        print "moving ten_min_median of " + rad + " into ten_min_median db"
        combine_ten_min_median(rad, ftype=ftype, 
                               baseLocation=baseLocation, dbName=None)
        print "finished moving"

    return
if __name__ == "__main__":
    import multiprocessing

    seasons = ["winter", "summer", "equinox"]
    jobs = []
    for season in seasons:
        #baseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_mlt" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_geo" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_2/" + "data_in_geo" + "/"
        baseLocation="../data/sqlite3/" + season + "/kp_l_1/" + "data_in_mlt" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_1/" + "data_in_geo" + "/"
        p = multiprocessing.Process(target=worker, args=(baseLocation,))
        jobs.append(p)
        p.start()

