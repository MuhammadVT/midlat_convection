def ten_min_median(rad, stm, etm, ftype="fitacf", 
                   baseLocation="../data/sqlite3/", dopsearch_dbName=None):
    
    """ bins ten minutes of data from a single radar by choosing 
    the median vector in each azimuth bin within each grid cell. 
    The results are stored in a different db file that combines all the beams's 
    data into a single table named by the radar.
    """
    import sqlite3
    import datetime as dt
    import numpy as np 

    # make db connection for dopsearch
    if dopsearch_dbName is None:
        dopsearch_dbName = "dopsearch_" + rad + "_" + ftype + ".sqlite"
    conn_dop = sqlite3.connect(baseLocation + dopsearch_dbName)
    cur_dop = conn_dop.cursor()

    # create a db that stores each radar data as a table
    conn = sqlite3.connect(baseLocation + "ten_min_median_" + rad + "_" + ftype + ".sqlite")
    cur = conn.cursor()

    # create a table that combines all the beams of a radar
    cur.execute("CREATE TABLE IF NOT EXISTS {tb}\
                (vel REAL, glatc REAL, glonc REAL, gazmc REAL, datetime TIMESTAMP)"\
                .format(tb=rad))

    # get all the table names
    cur_dop.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
    tbl_names = cur_dop.fetchall()
    tbl_names = [x[0] for x in tbl_names]

    
    # lengh of the time interval during which the data is median filtered 
    len_tm = 10    # minutes

    # initial starting and ending time of the time interval given by len_tm
    sdtm = stm
    edtm = sdtm + dt.timedelta(minutes=len_tm)

    # initilize parameters
    # bin_vel stores the velocity data as {glatc-glonc-gazmc: [velocites]}  
    while edtm <= etm:
        bin_vel = {}
        for tbl in tbl_names:
            # select column variables from dopsearch for a ten minute interval
            command = "SELECT vel, glatc, glonc, gazmc FROM {tb} WHERE datetime >= ?\
                      AND datetime < ?".format(tb=tbl)
            cur_dop.execute(command, (sdtm, edtm))
            rows_tmp = cur_dop.fetchall()

            if rows_tmp:
                #rows_tmp = [x[0] for x in rows_tmp]
                for row in rows_tmp:
                    vel, lat, lon, az = row
                    
                    if None not in row:
                        # convert from string to float
                        vel = [float(x) for x in vel.split(",")]
                        lat = [float(x) for x in lat.split(",")]
                        lon = [float(x) for x in lon.split(",")]
                        az = [float(x) for x in az.split(",")]

                        for i in range(len(vel)):
                            try:
                                bin_vel[(lat[i],lon[i],az[i])].append(vel[i])
                            except KeyError:
                                bin_vel[(lat[i],lon[i],az[i])] = [vel[i]]
                    else:
                        continue
            
            else:
                continue

        if bin_vel:

            # take the mid point of sdtm and edtm
            mid_tm = sdtm + dt.timedelta(minutes=len_tm/2.)
                
            # populate the rad table 
            for ky in bin_vel.keys(): 
                # take the median value
                bin_vel[ky] = round(np.median(bin_vel[ky]),2)

                cur.execute("INSERT INTO {tb} (vel, glatc, glonc, gazmc, datetime)\
                            VALUES (?, ?, ?, ?, ?)".format(tb=rad),\
                            (bin_vel[ky], ky[0], ky[1], ky[2], mid_tm))

        print "finished taking the median velocity between " + str(sdtm)\
                + " and " + str(edtm) + " of " + rad
        # update starting and ending time of the time interval given by len_tm
        sdtm = edtm
        edtm = sdtm + dt.timedelta(minutes=len_tm)

    # commit the change
    conn.commit()

    # close db connections
    conn.close()
    conn_dop.close()

    return

def worker(rads, season, baseLocation):

    import datetime as dt

    # input parameters
    ftype = "fitacf"

    # set the time interval for each season
    if season == "summer":
        stms = [dt.datetime(2011,5,1), dt.datetime(2012,5,1)]
        etms = [dt.datetime(2011,9,1), dt.datetime(2012,9,1)]

    if season == "winter":
        stms = [dt.datetime(2011,1,1), dt.datetime(2011,11,1), dt.datetime(2012,11,1)]
        etms = [dt.datetime(2011,3,1), dt.datetime(2012,3,1), dt.datetime(2013,1,1)]

    if season == "equinox":
        stms = [dt.datetime(2011,3,1), dt.datetime(2011,9,1), dt.datetime(2012,3,1), dt.datetime(2012,9,1)]
        etms = [dt.datetime(2011,5,1), dt.datetime(2011,11,1), dt.datetime(2012,5,1), dt.datetime(2012,11,1)]

    # loop through time interval
    for dd in range(len(stms)):
        stm = stms[dd]
        etm = etms[dd]
                
        # loop through radars
        for rad in rads:

            # take ten minutes median values
            print "start working on table " + rad + " for interval between " + str(stm) + " and " + str(etm)
            ten_min_median(rad, stm, etm, ftype=ftype, 
                           baseLocation=baseLocation, dopsearch_dbName=None)
            print "finish taking ten mimute median filtering on " + rad + " for interval between " + str(stm) + " and " + str(etm)

    return
if __name__ == "__main__":
    import multiprocessing

    seasons = ["winter", "summer", "equinox"]
    rads_list = [["bks", "wal", "fhe"], ["fhw", "cve", "cvw"]]
    #seasons = ["equinox"]
    #rads_list = [["cvw"]]
    
    jobs = []
    for season in seasons:
        #baseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_mlt" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_3/" + "data_in_geo" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_2/" + "data_in_geo" + "/"
        #baseLocation="../data/sqlite3/" + season + "/kp_l_1/" + "data_in_geo" + "/"
        baseLocation="../data/sqlite3/" + season + "/kp_l_1/" + "data_in_mlt" + "/"
        for i in range(len(rads_list)):
            p = multiprocessing.Process(target=worker, args=(rads_list[i], season, baseLocation))
            jobs.append(p)
            p.start()
